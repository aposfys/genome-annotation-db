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
    vm_steps: int = -1

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
    def work_ratio(self) -> float:
        """Speed-up expressed in engine work rather than elapsed time.

        Wall-clock timings depend on the machine; VM-instruction counts do not,
        so this is the figure that reproduces on someone else's hardware.
        """
        if self.with_indexes.vm_steps <= 0 or self.without.vm_steps <= 0:
            return float("nan")
        return self.without.vm_steps / self.with_indexes.vm_steps

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
            "vm_steps_without_indexes": self.without.vm_steps,
            "vm_steps_with_indexes": self.with_indexes.vm_steps,
            "work_ratio": (
                None if self.work_ratio != self.work_ratio else round(self.work_ratio, 2)
            ),
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

    # Counted in a separate pass: the progress handler that makes it possible
    # also makes execution far slower, so it must never run inside a timed one.
    # The first count warms the page cache and reads a little high, so a median
    # of three is used.
    steps = statistics.median(
        [conn.count_vm_steps(statement, query.defaults) for _ in range(3)]
    )

    return Timing(
        query=query.name,
        title=query.title,
        seconds=statistics.median(durations),
        rows=len(rows),
        plan=explain(conn, statement, query.defaults),
        vm_steps=int(steps),
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
        "total_vm_steps_without_indexes": sum(
            c.without.vm_steps for c in comparisons if c.without.vm_steps > 0
        ),
        "total_vm_steps_with_indexes": sum(
            c.with_indexes.vm_steps for c in comparisons if c.with_indexes.vm_steps > 0
        ),
    }


def render(comparisons: Sequence[Comparison]) -> str:
    """A Markdown table of the comparison, for the README."""
    header = (
        "| Query | Rows | Without indexes | With indexes | Speed-up | VM steps (without → with) | Work ratio |\n"
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n"
    )
    body = "".join(
        f"| {c.without.query} — {c.without.title} | {c.with_indexes.rows} "
        f"| {c.without.milliseconds:,.1f} ms | {c.with_indexes.milliseconds:,.1f} ms "
        f"| {c.speedup:.1f}× "
        f"| {c.without.vm_steps:,} → {c.with_indexes.vm_steps:,} "
        f"| {c.work_ratio:.1f}× |\n"
        for c in comparisons
    )
    return header + body


# --- interval strategies -----------------------------------------------------


@dataclass(frozen=True)
class IntervalResult:
    """How one interval strategy performed over a set of query windows."""

    strategy: str
    description: str
    windows: int
    total_seconds: float
    rows_returned: int
    plan: str

    @property
    def per_query_ms(self) -> float:
        return (self.total_seconds / self.windows) * 1000 if self.windows else 0.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "strategy": self.strategy,
            "description": self.description,
            "windows": self.windows,
            "per_query_ms": round(self.per_query_ms, 4),
            "total_ms": round(self.total_seconds * 1000, 2),
            "rows_returned": self.rows_returned,
            "plan": self.plan,
        }


def compare_intervals(
    conn: Connection,
    window_count: int = 200,
    width: int = 1_000_000,
    repeats: int = 3,
) -> tuple[list[IntervalResult], dict[str, Any]]:
    """Time each interval strategy over the same set of random windows.

    Correctness is checked before timing: a faster structure that returns
    different rows is not a faster structure.
    """
    from . import intervals

    available = [
        s for s in intervals.STRATEGIES if not s.requires_sqlite or conn.backend.is_sqlite
    ]
    query_windows = list(intervals.windows(conn, window_count, width))

    for strategy in available:
        intervals.build(conn, strategy.name)

    # Every strategy must agree with the baseline, window for window.
    baseline = {
        (chrom, start, end): intervals.query(conn, "btree", chrom, start, end)
        for chrom, start, end in query_windows
    }
    disagreements: list[str] = []
    for strategy in available[1:]:
        for (chrom, start, end), expected in baseline.items():
            got = intervals.query(conn, strategy.name, chrom, start, end)
            if [row[0] for row in got] != [row[0] for row in expected]:
                disagreements.append(f"{strategy.name} at {chrom}:{start}-{end}")

    results: list[IntervalResult] = []
    for strategy in available:
        elapsed = 0.0
        returned = 0
        for _ in range(repeats):
            start_time = time.perf_counter()
            for chrom, start, end in query_windows:
                returned += len(intervals.query(conn, strategy.name, chrom, start, end))
            elapsed += time.perf_counter() - start_time

        sample = query_windows[0]
        results.append(
            IntervalResult(
                strategy=strategy.name,
                description=strategy.description,
                windows=len(query_windows) * repeats,
                total_seconds=elapsed,
                rows_returned=returned,
                plan=intervals.plan(conn, strategy.name, *sample),
            )
        )

    fastest = min(results, key=lambda r: r.per_query_ms)
    summary = {
        "windows": len(query_windows),
        "window_width_bp": width,
        "repeats": repeats,
        "strategies_agree": not disagreements,
        "disagreements": disagreements[:5],
        "fastest": fastest.strategy,
        "region_stats": intervals.region_stats(conn),
        "speedup_over_btree": {
            r.strategy: round(
                next(x for x in results if x.strategy == "btree").per_query_ms
                / r.per_query_ms,
                2,
            )
            for r in results
            if r.per_query_ms > 0
        },
    }
    return results, summary


def render_intervals(results: Sequence[IntervalResult]) -> str:
    header = (
        "| Strategy | Per query | Total | Rows | Structure |\n"
        "| --- | ---: | ---: | ---: | --- |\n"
    )
    body = "".join(
        f"| {r.strategy} | {r.per_query_ms:.3f} ms | {r.total_seconds * 1000:,.0f} ms "
        f"| {r.rows_returned:,} | {r.description} |\n"
        for r in results
    )
    return header + body
