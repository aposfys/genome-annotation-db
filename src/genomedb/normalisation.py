"""Formal normalisation, and what it costs.

Calling a schema "normalised" is an assertion. This module states the functional
dependencies explicitly, checks Boyce-Codd Normal Form against them mechanically,
and then measures the price of the property rather than assuming it is free.

BCNF holds when, for every non-trivial functional dependency X -> Y, X is a
superkey. That is a stricter condition than third normal form, and it is the one
worth aiming at: a schema in 3NF but not BCNF can still hold a redundancy that
lets two rows disagree with each other.

The dependencies are declared rather than inferred. Inferring them from data
would only report what happens to hold in this snapshot, which is a different
claim from what the schema guarantees.
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from .db import Connection


@dataclass(frozen=True)
class Dependency:
    """A functional dependency X -> Y within one table."""

    table: str
    determinant: tuple[str, ...]
    dependent: tuple[str, ...]

    def __str__(self) -> str:
        return f"{self.table}: {', '.join(self.determinant)} -> {', '.join(self.dependent)}"


@dataclass(frozen=True)
class Relation:
    """A table, its candidate keys, and the dependencies that hold on it."""

    name: str
    attributes: tuple[str, ...]
    candidate_keys: tuple[tuple[str, ...], ...]
    dependencies: tuple[Dependency, ...]

    def is_superkey(self, attributes: Sequence[str]) -> bool:
        """True if the given attributes contain some candidate key."""
        given = set(attributes)
        return any(set(key) <= given for key in self.candidate_keys)

    def violations(self) -> list[Dependency]:
        """Dependencies whose determinant is not a superkey.

        A trivial dependency -- one whose dependent side is contained in its
        determinant -- cannot violate BCNF and is skipped.
        """
        return [
            dependency
            for dependency in self.dependencies
            if not set(dependency.dependent) <= set(dependency.determinant)
            and not self.is_superkey(dependency.determinant)
        ]

    @property
    def in_bcnf(self) -> bool:
        return not self.violations()


# Every dependency the schema is intended to enforce. Each is keyed on the
# table's primary key, which is what makes the BCNF argument short.
SCHEMA: tuple[Relation, ...] = (
    Relation(
        name="gene",
        attributes=(
            "gene_id",
            "gene_name",
            "chrom",
            "gene_start",
            "gene_end",
            "strand",
            "biotype",
        ),
        candidate_keys=(("gene_id",),),
        dependencies=(
            Dependency(
                "gene",
                ("gene_id",),
                ("gene_name", "chrom", "gene_start", "gene_end", "strand", "biotype"),
            ),
        ),
    ),
    Relation(
        name="transcript",
        attributes=("transcript_id", "gene_id", "tx_start", "tx_end", "tx_length"),
        candidate_keys=(("transcript_id",),),
        dependencies=(
            Dependency(
                "transcript",
                ("transcript_id",),
                ("gene_id", "tx_start", "tx_end", "tx_length"),
            ),
        ),
    ),
    Relation(
        name="exon",
        attributes=("exon_id", "exon_start", "exon_end"),
        candidate_keys=(("exon_id",),),
        dependencies=(Dependency("exon", ("exon_id",), ("exon_start", "exon_end")),),
    ),
    Relation(
        name="go_term",
        attributes=("go_id", "go_name", "go_namespace"),
        candidate_keys=(("go_id",),),
        dependencies=(Dependency("go_term", ("go_id",), ("go_name", "go_namespace")),),
    ),
    Relation(
        # The junction table that makes alternative splicing representable. Its
        # only non-key attribute depends on the whole composite key: an exon's
        # rank is a property of the pairing, not of the exon alone -- the same
        # exon can be third in one transcript and first in another.
        name="transcript_exon",
        attributes=("transcript_id", "exon_id", "exon_rank"),
        candidate_keys=(("transcript_id", "exon_id"),),
        dependencies=(
            Dependency("transcript_exon", ("transcript_id", "exon_id"), ("exon_rank",)),
        ),
    ),
    Relation(
        name="gene_go",
        attributes=("gene_id", "go_id"),
        candidate_keys=(("gene_id", "go_id"),),
        dependencies=(),  # all-key relation: nothing to depend on
    ),
)


def check() -> dict[str, Any]:
    """Verify BCNF across the schema."""
    per_relation = {
        relation.name: {
            "candidate_keys": [list(key) for key in relation.candidate_keys],
            "dependencies": [str(d) for d in relation.dependencies],
            "violations": [str(v) for v in relation.violations()],
            "in_bcnf": relation.in_bcnf,
        }
        for relation in SCHEMA
    }
    return {
        "relations": per_relation,
        "all_in_bcnf": all(r.in_bcnf for r in SCHEMA),
        "violating_relations": [r.name for r in SCHEMA if not r.in_bcnf],
    }


def verify_against_data(conn: Connection) -> dict[str, Any]:
    """Check the declared dependencies actually hold in the loaded data.

    A dependency the schema does not enforce can still be violated by a load
    bug. This counts determinant values mapping to more than one dependent
    value, which is what a violation looks like in practice.
    """
    findings: dict[str, int | str] = {}
    for relation in SCHEMA:
        for dependency in relation.dependencies:
            determinant = ", ".join(dependency.determinant)
            dependent = ", ".join(dependency.dependent)
            sql = (
                f"SELECT COUNT(*) FROM (SELECT {determinant} FROM {relation.name}"
                f" GROUP BY {determinant}"
                f" HAVING COUNT(DISTINCT {dependent}) > 1) AS violations"
            )
            try:
                count = int(conn.scalar(sql) or 0)
            except Exception as error:
                findings[str(dependency)] = f"uncheckable: {error}"
                continue
            findings[str(dependency)] = count
    return {
        "checked": len(findings),
        "violations": {k: v for k, v in findings.items() if v not in (0, "0")},
        "all_hold": all(v == 0 for v in findings.values() if isinstance(v, int)),
    }


# --- what normalisation costs ------------------------------------------------

NORMALISED_EXON_COUNTS = """
SELECT t.transcript_id, COUNT(te.exon_id) AS n_exons
FROM transcript AS t
LEFT JOIN transcript_exon AS te ON te.transcript_id = t.transcript_id
GROUP BY t.transcript_id
"""

DENORMALISED_EXON_COUNTS = "SELECT transcript_id, n_exons FROM transcript_denorm"


def measure_denormalisation(conn: Connection, repeats: int = 20) -> dict[str, Any]:
    """Time an aggregate against a materialised column holding the same answer.

    Counting a transcript's exons means joining the junction table every time.
    Storing the count on ``transcript`` removes the join, and introduces the
    possibility of the two disagreeing. This measures the gain so the trade can
    be made on evidence.
    """
    conn.execute("DROP TABLE IF EXISTS transcript_denorm")
    conn.execute(
        "CREATE TABLE transcript_denorm ("
        "  transcript_id VARCHAR(30) PRIMARY KEY, n_exons INTEGER NOT NULL)"
    )
    _, counts = conn.query(NORMALISED_EXON_COUNTS)
    conn.executemany(
        "INSERT INTO transcript_denorm (transcript_id, n_exons) VALUES (?, ?)", counts
    )
    conn.commit()

    def timed(sql: str) -> float:
        durations = []
        for _ in range(repeats):
            start = time.perf_counter()
            conn.query(sql)
            durations.append(time.perf_counter() - start)
        durations.sort()
        return durations[len(durations) // 2]

    normalised = timed(NORMALISED_EXON_COUNTS)
    denormalised = timed(DENORMALISED_EXON_COUNTS)

    # The cost of the shortcut: the stored value can drift from the truth.
    drift = conn.scalar(
        "SELECT COUNT(*) FROM transcript_denorm d"
        " JOIN (" + NORMALISED_EXON_COUNTS + ") AS live"
        "   ON live.transcript_id = d.transcript_id"
        " WHERE live.n_exons <> d.n_exons"
    )

    return {
        "rows": len(counts),
        "normalised_ms": round(normalised * 1000, 3),
        "denormalised_ms": round(denormalised * 1000, 3),
        "speedup": round(normalised / denormalised, 1) if denormalised else None,
        "rows_disagreeing_now": int(drift or 0),
        "trade_off": (
            "The materialised column answers without the join, but nothing in the"
            " schema keeps it true: any write to transcript_exon that does not"
            " also update it leaves the two disagreeing, and no constraint can"
            " express that."
        ),
    }


def render(bcnf: dict[str, Any], denorm: dict[str, Any]) -> str:
    lines = [
        "| Relation | Candidate key | Dependencies | BCNF |",
        "| --- | --- | ---: | --- |",
    ]
    for name, detail in bcnf["relations"].items():
        key = ", ".join(detail["candidate_keys"][0])
        lines.append(
            f"| `{name}` | ({key}) | {len(detail['dependencies'])} "
            f"| {'yes' if detail['in_bcnf'] else 'NO'} |"
        )
    lines += [
        "",
        f"Normalised aggregate: {denorm['normalised_ms']} ms · "
        f"materialised column: {denorm['denormalised_ms']} ms · "
        f"**{denorm['speedup']}× faster**",
    ]
    return "\n".join(lines) + "\n"
