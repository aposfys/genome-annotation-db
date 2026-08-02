"""Integrity and consistency checks on the loaded database.

Foreign keys guarantee referential integrity, so orphan rows cannot exist if the
constraints were enforced. These checks confirm that they were -- SQLite accepts
`PRAGMA foreign_keys = OFF` silently, and a load run that way looks identical
until something reads it -- and then go on to the things constraints cannot
express, such as whether exon ranks are contiguous.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .db import Connection


@dataclass(frozen=True)
class Check:
    """One integrity question and the answer it must give."""

    name: str
    sql: str
    expected: int = 0
    description: str = ""

    def run(self, conn: Connection) -> CheckResult:
        observed = conn.scalar(self.sql) or 0
        return CheckResult(check=self, observed=int(observed))


@dataclass(frozen=True)
class CheckResult:
    check: Check
    observed: int

    @property
    def passed(self) -> bool:
        return self.observed == self.check.expected

    def as_dict(self) -> dict[str, Any]:
        return {
            "check": self.check.name,
            "description": self.check.description,
            "expected": self.check.expected,
            "observed": self.observed,
            "passed": self.passed,
        }


INTEGRITY_CHECKS: tuple[Check, ...] = (
    Check(
        "orphan_transcripts",
        "SELECT COUNT(*) FROM transcript t"
        " LEFT JOIN gene g ON g.gene_id = t.gene_id WHERE g.gene_id IS NULL",
        description="transcripts whose gene is missing",
    ),
    Check(
        "orphan_transcript_exon",
        "SELECT COUNT(*) FROM transcript_exon te"
        " LEFT JOIN transcript t ON t.transcript_id = te.transcript_id"
        " LEFT JOIN exon e ON e.exon_id = te.exon_id"
        " WHERE t.transcript_id IS NULL OR e.exon_id IS NULL",
        description="junction rows pointing at a missing transcript or exon",
    ),
    Check(
        "orphan_gene_go",
        "SELECT COUNT(*) FROM gene_go gg"
        " LEFT JOIN gene g ON g.gene_id = gg.gene_id"
        " LEFT JOIN go_term gt ON gt.go_id = gg.go_id"
        " WHERE g.gene_id IS NULL OR gt.go_id IS NULL",
        description="GO assignments pointing at a missing gene or term",
    ),
    Check(
        "invalid_go_namespace",
        "SELECT COUNT(*) FROM go_term WHERE go_namespace NOT IN ('BP','MF','CC')",
        description="GO terms outside the three namespaces",
    ),
    Check(
        "invalid_strand",
        "SELECT COUNT(*) FROM gene WHERE strand NOT IN (-1, 1)",
        description="genes with a strand other than +1 or -1",
    ),
    Check(
        "inverted_gene_spans",
        "SELECT COUNT(*) FROM gene WHERE gene_end < gene_start",
        description="genes whose end precedes their start",
    ),
    Check(
        "transcripts_outside_gene_span",
        "SELECT COUNT(*) FROM transcript t JOIN gene g ON g.gene_id = t.gene_id"
        " WHERE t.tx_start < g.gene_start OR t.tx_end > g.gene_end",
        description="transcripts extending beyond their gene's annotated span",
    ),
    Check(
        "transcripts_without_exons",
        "SELECT COUNT(*) FROM transcript t"
        " LEFT JOIN transcript_exon te ON te.transcript_id = t.transcript_id"
        " WHERE te.transcript_id IS NULL",
        description="transcripts with no exons recorded",
    ),
    Check(
        "non_contiguous_exon_ranks",
        "SELECT COUNT(*) FROM ("
        "  SELECT transcript_id FROM transcript_exon"
        "  GROUP BY transcript_id"
        "  HAVING MAX(exon_rank) <> COUNT(*)"
        ") AS t",
        description="transcripts whose exon ranks are not 1..n",
    ),
)

COUNTS = {
    "gene": "SELECT COUNT(*) FROM gene",
    "transcript": "SELECT COUNT(*) FROM transcript",
    "exon": "SELECT COUNT(*) FROM exon",
    "go_term": "SELECT COUNT(*) FROM go_term",
    "transcript_exon": "SELECT COUNT(*) FROM transcript_exon",
    "gene_go": "SELECT COUNT(*) FROM gene_go",
}


def table_counts(conn: Connection) -> dict[str, int]:
    return {name: int(conn.scalar(sql) or 0) for name, sql in COUNTS.items()}


def run_checks(conn: Connection) -> list[CheckResult]:
    return [check.run(conn) for check in INTEGRITY_CHECKS]


def report(conn: Connection) -> dict[str, Any]:
    results = run_checks(conn)
    return {
        "counts": table_counts(conn),
        "checks": [result.as_dict() for result in results],
        "all_passed": all(result.passed for result in results),
        "failed": [r.check.name for r in results if not r.passed],
    }
