"""Load Ensembl BioMart exports into the database.

The loader validates as it goes and reports what it rejected. Silently skipping
malformed rows is the failure mode worth guarding against here: roughly 10% of
the GO export carries an accession with no domain, and a loader that drops those
without saying so leaves a database that looks complete and is not.
"""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from .db import Connection

# BioMart's GO domain values, mapped to the two-letter namespace the schema uses.
GO_NAMESPACES = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}

# Insert in batches rather than one row at a time; 97k transcript_exon rows is
# enough for the round-trip cost to dominate otherwise.
BATCH_SIZE = 5_000

TABLES_IN_DEPENDENCY_ORDER = (
    "gene_go",
    "transcript_exon",
    "go_term",
    "exon",
    "transcript",
    "gene",
)


@dataclass
class LoadReport:
    """What the loader inserted, and what it refused."""

    inserted: Counter = field(default_factory=Counter)
    rejected: Counter = field(default_factory=Counter)
    examples: dict[str, str] = field(default_factory=dict)

    def reject(self, reason: str, detail: str) -> None:
        self.rejected[reason] += 1
        self.examples.setdefault(reason, detail)

    @property
    def total_rejected(self) -> int:
        return sum(self.rejected.values())

    def as_dict(self) -> dict[str, object]:
        return {
            "inserted": dict(self.inserted),
            "rejected": dict(self.rejected),
            "rejection_examples": self.examples,
            "total_inserted": sum(self.inserted.values()),
            "total_rejected": self.total_rejected,
        }

    def render(self) -> str:
        lines = ["Inserted:"]
        lines += [
            f"  {table:18s} {count:>8,}" for table, count in sorted(self.inserted.items())
        ]
        if self.rejected:
            lines.append("Rejected:")
            for reason, count in sorted(self.rejected.items(), key=lambda kv: -kv[1]):
                lines.append(f"  {reason:40s} {count:>8,}   e.g. {self.examples[reason]}")
        else:
            lines.append("Rejected: nothing")
        return "\n".join(lines)


def _rows(handle: Iterable[str]) -> Iterator[dict[str, str]]:
    # QUOTE_NONE because BioMart emits true TSV, not quoted CSV. GO definitions
    # are written as "text" [source] and a reader that honours quoting can
    # swallow the delimiters that follow one.
    yield from csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)


def open_table(path: Path) -> Iterator[dict[str, str]]:
    """Read a TSV, transparently handling gzip."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield from _rows(handle)
    else:
        with path.open(newline="", encoding="utf-8") as handle:
            yield from _rows(handle)


def _batched(rows: list[tuple], size: int = BATCH_SIZE) -> Iterator[list[tuple]]:
    for start in range(0, len(rows), size):
        yield rows[start : start + size]


def _insert(conn: Connection, sql: str, rows: list[tuple]) -> int:
    total = 0
    for batch in _batched(rows):
        total += max(conn.executemany(sql, batch), 0) or len(batch)
    conn.commit()
    return len(rows)


def create_schema(conn: Connection, sql_dir: Path, with_indexes: bool = True) -> None:
    """Drop everything and rebuild from the DDL files."""
    for table in TABLES_IN_DEPENDENCY_ORDER:
        conn.execute(f"DROP TABLE IF EXISTS {table}")
    conn.commit()

    conn.executescript((sql_dir / "schema.sql").read_text(encoding="utf-8"))
    if with_indexes:
        conn.executescript((sql_dir / "indexes.sql").read_text(encoding="utf-8"))
    conn.commit()


def load_genes(conn: Connection, path: Path, report: LoadReport) -> set[str]:
    """Load the gene table; returns the accepted gene identifiers."""
    rows: list[tuple] = []
    accepted: set[str] = set()

    for record in open_table(path):
        gene_id = record["Gene stable ID"].strip()
        strand = record["Strand"].strip()

        if strand not in {"1", "-1"}:
            report.reject("gene: strand not 1 or -1", f"{gene_id} strand={strand!r}")
            continue
        if gene_id in accepted:
            report.reject("gene: duplicate gene_id", gene_id)
            continue

        rows.append(
            (
                gene_id,
                record["Gene name"].strip() or None,
                record["Chromosome/scaffold name"].strip(),
                int(record["Gene start (bp)"]),
                int(record["Gene end (bp)"]),
                int(strand),
                record["Gene type"].strip(),
            )
        )
        accepted.add(gene_id)

    report.inserted["gene"] = _insert(
        conn,
        "INSERT INTO gene (gene_id, gene_name, chrom, gene_start, gene_end, strand, biotype)"
        " VALUES (?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    return accepted


def load_transcripts_and_exons(
    conn: Connection, path: Path, genes: set[str], report: LoadReport
) -> None:
    """Load transcripts, exons and the junction between them from one export.

    The export is one row per exon-in-transcript, so transcripts and exons both
    repeat and have to be de-duplicated on the way in.
    """
    transcripts: dict[str, tuple] = {}
    exons: dict[str, tuple] = {}
    pairs: dict[tuple[str, str], tuple] = {}

    for record in open_table(path):
        gene_id = record["Gene stable ID"].strip()
        transcript_id = record["Transcript stable ID"].strip()
        exon_id = record["Exon stable ID"].strip()

        # A transcript whose gene was rejected would violate the foreign key.
        if gene_id not in genes:
            report.reject(
                "transcript: gene_id not in gene table", f"{transcript_id} -> {gene_id}"
            )
            continue

        if transcript_id not in transcripts:
            length = record.get("Transcript length (including UTRs and CDS)", "").strip()
            transcripts[transcript_id] = (
                transcript_id,
                gene_id,
                int(record["Transcript start (bp)"]),
                int(record["Transcript end (bp)"]),
                int(length) if length.isdigit() else None,
            )

        if exon_id not in exons:
            exons[exon_id] = (
                exon_id,
                int(record["Exon region start (bp)"]),
                int(record["Exon region end (bp)"]),
            )

        key = (transcript_id, exon_id)
        if key in pairs:
            report.reject("transcript_exon: duplicate pair", f"{transcript_id}/{exon_id}")
            continue
        pairs[key] = (transcript_id, exon_id, int(record["Exon rank in transcript"]))

    report.inserted["transcript"] = _insert(
        conn,
        "INSERT INTO transcript (transcript_id, gene_id, tx_start, tx_end, tx_length)"
        " VALUES (?, ?, ?, ?, ?)",
        list(transcripts.values()),
    )
    report.inserted["exon"] = _insert(
        conn,
        "INSERT INTO exon (exon_id, exon_start, exon_end) VALUES (?, ?, ?)",
        list(exons.values()),
    )
    report.inserted["transcript_exon"] = _insert(
        conn,
        "INSERT INTO transcript_exon (transcript_id, exon_id, exon_rank) VALUES (?, ?, ?)",
        list(pairs.values()),
    )


def load_go(conn: Connection, path: Path, genes: set[str], report: LoadReport) -> None:
    """Load GO terms and their gene assignments."""
    terms: dict[str, tuple] = {}
    assignments: set[tuple[str, str]] = set()

    for record in open_table(path):
        gene_id = record["Gene stable ID"].strip()
        go_id = record["GO term accession"].strip()
        domain = record["GO domain"].strip()

        if not go_id:
            report.reject("go: no accession", f"gene {gene_id}")
            continue
        if domain not in GO_NAMESPACES:
            report.reject(
                "go: unmapped or missing domain", f"{go_id} domain={domain or 'empty'!r}"
            )
            continue
        if gene_id not in genes:
            report.reject("go: gene_id not in gene table", f"{go_id} -> {gene_id}")
            continue

        terms.setdefault(
            go_id,
            (go_id, record["GO term name"].strip() or go_id, GO_NAMESPACES[domain]),
        )
        assignments.add((gene_id, go_id))

    report.inserted["go_term"] = _insert(
        conn,
        "INSERT INTO go_term (go_id, go_name, go_namespace) VALUES (?, ?, ?)",
        list(terms.values()),
    )
    report.inserted["gene_go"] = _insert(
        conn,
        "INSERT INTO gene_go (gene_id, go_id) VALUES (?, ?)",
        sorted(assignments),
    )


def build(
    conn: Connection,
    data_dir: Path,
    sql_dir: Path,
    with_indexes: bool = True,
) -> LoadReport:
    """Create the schema and load every table."""
    report = LoadReport()
    create_schema(conn, sql_dir, with_indexes=with_indexes)

    genes = load_genes(conn, _find(data_dir, "genes"), report)
    load_transcripts_and_exons(conn, _find(data_dir, "transcripts_exons"), genes, report)
    load_go(conn, _find(data_dir, "go"), genes, report)
    return report


def _find(data_dir: Path, stem: str) -> Path:
    """Accept either the gzipped or plain form of an export."""
    for candidate in (data_dir / f"{stem}.tsv.gz", data_dir / f"{stem}.tsv"):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"neither {stem}.tsv.gz nor {stem}.tsv found in {data_dir}")
