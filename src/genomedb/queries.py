"""The analysis query set, and how to run it."""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .db import Connection


@dataclass(frozen=True)
class Query:
    """One named analysis query, with its default parameters."""

    name: str
    title: str
    defaults: tuple[Any, ...] = ()

    @property
    def filename(self) -> str:
        return f"{self.name}.sql"

    def sql(self, sql_dir: Path) -> str:
        return (sql_dir / "queries" / self.filename).read_text(encoding="utf-8")


# Parameterised where the original hardcoded a literal, so the same query can be
# pointed at a different gene or GO term without editing the file.
QUERIES: tuple[Query, ...] = (
    Query("Q1", "Transcripts of a named gene", ("DEFB125",)),
    Query("Q2", "Most heavily transcribed genes"),
    Query("Q3", "Transcripts with the most exons"),
    Query("Q4", "Genes carrying at least 20 GO terms"),
    Query("Q5", "Annotation depth per GO namespace"),
    Query("Q6", "Structure of genes carrying a GO term", ("GO:0007186",)),
    Query("Q7", "Exon reuse across transcripts, per gene"),
)

BY_NAME = {query.name: query for query in QUERIES}


def run(
    conn: Connection,
    query: Query,
    sql_dir: Path,
    params: Sequence[Any] | None = None,
) -> tuple[list[str], list[tuple]]:
    """Execute one query and return its columns and rows."""
    statement = query.sql(sql_dir).strip().rstrip(";")
    values = tuple(params) if params is not None else query.defaults
    return conn.query(statement, values)


def write_tsv(columns: Sequence[str], rows: Sequence[tuple], path: Path) -> Path:
    """Write a result set, always with a header row.

    The original results were inconsistent: five files carried a header and two
    did not, which makes them awkward to read back programmatically.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(["" if value is None else value for value in row] for row in rows)
    return path


def run_all(conn: Connection, sql_dir: Path, results_dir: Path) -> dict[str, dict[str, Any]]:
    """Run every query and write each result to results/<name>.tsv."""
    summary: dict[str, dict[str, Any]] = {}
    for query in QUERIES:
        columns, rows = run(conn, query, sql_dir)
        write_tsv(columns, rows, results_dir / f"{query.name}.tsv")
        summary[query.name] = {
            "title": query.title,
            "columns": columns,
            "rows": len(rows),
            "parameters": list(query.defaults),
        }
    return summary


# --- ad hoc lookups, previously the interactive menu -------------------------

GENE_LOOKUP = """
SELECT g.gene_id, g.gene_name, t.transcript_id, COUNT(te.exon_id) AS n_exons
FROM gene AS g
JOIN transcript AS t ON t.gene_id = g.gene_id
LEFT JOIN transcript_exon AS te ON te.transcript_id = t.transcript_id
WHERE g.gene_name = ?
GROUP BY g.gene_id, g.gene_name, t.transcript_id
ORDER BY t.transcript_id
"""

GO_LOOKUP = """
SELECT g.gene_id, g.gene_name, COUNT(DISTINCT t.transcript_id) AS n_transcripts
FROM gene_go AS gg
JOIN gene AS g ON g.gene_id = gg.gene_id
JOIN transcript AS t ON t.gene_id = gg.gene_id
WHERE gg.go_id = ?
GROUP BY g.gene_id, g.gene_name
ORDER BY n_transcripts DESC, g.gene_id
"""


def gene(conn: Connection, gene_name: str) -> tuple[list[str], list[tuple]]:
    return conn.query(GENE_LOOKUP, (gene_name,))


def go(conn: Connection, go_id: str) -> tuple[list[str], list[tuple]]:
    return conn.query(GO_LOOKUP, (go_id,))
