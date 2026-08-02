"""Measure what the secondary indexes are worth.

Indexes are usually added on the assumption that they help. This measures it:
the same seven queries run against the same data with the indexes dropped and
again with them present, and the query planner is asked what it intends to do in
each case.

The comparison is fair because only the indexes change. Both passes read an
identical database, and each query is repeated so a single slow first run does
not decide the result.
"""

from __future__ import annotations

import statistics
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .db import Connection
from .queries import QUERIES, Query, run

# Repeats per query. Enough to see past scheduler noise without making the whole
# benchmark slow; the median is reported rather than the mean so one outlier
# cannot move the number.
REPEATS = 5


@dataclass(frozen=True)
class Timing:
    """How long one query took, and what it returned."""

    query: str
    title: str
    seconds: float
    rows: int
    plan: str

    @property
    def milliseconds(self) -> float:
        return self.seconds * 1000


@dataclass(frozen=True)
class Comparison:
    """One query, timed with and without indexes."""

    without: Timing
    with_indexes: Timing

    @property
    def speedup(self) -> float:
        if self.with_indexes.seconds <= 0:
            return float("inf")
        return self.without.seconds / self.with_indexes.seconds

    @property
    def uses_an_index(self) -> bool:
        """Whether the plan adopts one of *our* secondary indexes.

        Matching on "USING INDEX" alone would be satisfied by SQLite's automatic
        primary-key indexes, which are present either way, so this looks for the
        ``idx_`` prefix that indexes.sql uses.
        """
        return "idx_" in self.with_indexes.plan

    @property
    def plan_changed(self) -> bool:
        """Whether the planner did anything differently once indexes existed."""
        return self.without.plan != self.with_indexes.plan

    def as_dict(self) -> dict[str, Any]:
        return {
            "query": self.without.query,
            "title": self.without.title,
            "rows": self.with_indexes.rows,
            "without_indexes_ms": round(self.without.milliseconds, 2),
            "with_indexes_ms": round(self.with_indexes.milliseconds, 2),
            "speedup": round(self.speedup, 2),
            "plan_uses_index": self.uses_an_index,
            "plan_changed": self.plan_changed,
            "plan_with_indexes": self.with_indexes.plan,
            "plan_without_indexes": self.without.plan,
        }


def explain(conn: Connection, sql: str, params: Sequence[Any]) -> str:
    """Ask the planner how it intends to answer the query."""
    keyword = "EXPLAIN QUERY PLAN" if conn.backend.is_sqlite else "EXPLAIN"
    try:
        _, rows = conn.query(f"{keyword} {sql}", params)
    except Exception as error:
        return f"(plan unavailable: {error})"

    lines = []
    for row in rows:
        text = " ".join(str(field) for field in row if field not in (None, ""))
        lines.append(text)
    return " | ".join(lines)


def time_query(
    conn: Connection, query: Query, sql_dir: Path, repeats: int = REPEATS
) -> Timing:
    """Run a query several times and report the median elapsed time."""
    statement = query.sql(sql_dir).strip().rstrip(";")
    durations: list[float] = []
    rows: list[tuple] = []

    for _ in range(repeats):
        start = time.perf_counter()
        _, rows = run(conn, query, sql_dir)
        durations.append(time.perf_counter() - start)

    return Timing(
        query=query.name,
        title=query.title,
        seconds=statistics.median(durations),
        rows=len(rows),
        plan=explain(conn, statement, query.defaults),
    )


def drop_indexes(conn: Connection, sql_dir: Path) -> list[str]:
    """Drop every index defined in indexes.sql, returning their names."""
    text = (sql_dir / "indexes.sql").read_text(encoding="utf-8")
    names = [
        line.split("CREATE INDEX")[1].split("ON")[0].strip()
        for line in text.splitlines()
        if line.strip().upper().startswith("CREATE INDEX")
    ]
    for name in names:
        conn.execute(f"DROP INDEX IF EXISTS {name}")
    conn.commit()
    return names


def create_indexes(conn: Connection, sql_dir: Path) -> None:
    conn.executescript((sql_dir / "indexes.sql").read_text(encoding="utf-8"))
    conn.commit()


def compare(conn: Connection, sql_dir: Path, repeats: int = REPEATS) -> list[Comparison]:
    """Time every query without indexes, then with them.

    The unindexed pass runs first so the indexed pass cannot be flattered by a
    warm cache from its own earlier run.
    """
    names = drop_indexes(conn, sql_dir)
    print(f"Dropped {len(names)} indexes; timing without them...")
    without = {query.name: time_query(conn, query, sql_dir, repeats) for query in QUERIES}

    create_indexes(conn, sql_dir)
    print(f"Recreated {len(names)} indexes; timing with them...")
    with_indexes = {query.name: time_query(conn, query, sql_dir, repeats) for query in QUERIES}

    return [
        Comparison(without=without[query.name], with_indexes=with_indexes[query.name])
        for query in QUERIES
    ]


def summarise(comparisons: Sequence[Comparison]) -> dict[str, Any]:
    speedups = [c.speedup for c in comparisons if c.speedup != float("inf")]
    total_without = sum(c.without.seconds for c in comparisons)
    total_with = sum(c.with_indexes.seconds for c in comparisons)
    return {
        "queries": len(comparisons),
        "repeats": REPEATS,
        "total_without_indexes_ms": round(total_without * 1000, 2),
        "total_with_indexes_ms": round(total_with * 1000, 2),
        "overall_speedup": round(total_without / total_with, 2) if total_with else None,
        "median_speedup": round(statistics.median(speedups), 2) if speedups else None,
        "max_speedup": round(max(speedups), 2) if speedups else None,
        "queries_using_an_index": sum(1 for c in comparisons if c.uses_an_index),
        "queries_with_a_changed_plan": sum(1 for c in comparisons if c.plan_changed),
    }


def render(comparisons: Sequence[Comparison]) -> str:
    """A Markdown table of the comparison, for the README."""
    header = (
        "| Query | Rows | Without indexes | With indexes | Speed-up | Index used |\n"
        "| --- | ---: | ---: | ---: | ---: | --- |\n"
    )
    body = "".join(
        f"| {c.without.query} — {c.without.title} | {c.with_indexes.rows} "
        f"| {c.without.milliseconds:,.1f} ms | {c.with_indexes.milliseconds:,.1f} ms "
        f"| {c.speedup:.1f}× | {'yes' if c.uses_an_index else 'no'} |\n"
        for c in comparisons
    )
    return header + body
