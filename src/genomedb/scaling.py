"""Empirical complexity: how query cost grows with the size of the table.

A single timing says an index helped. It does not say *how* the cost grows, and
growth is what decides whether a design survives a hundredfold more data. This
module measures query time across a geometric series of table sizes and fits the
observed curve against the two candidate laws:

* a B-tree seek should be **O(log n)** -- time roughly linear in log n;
* a full scan should be **O(n)** -- time roughly linear in n.

Fitting both and comparing the residuals says which model the data actually
supports, rather than assuming the textbook answer.

Synthetic intervals are used rather than real annotation because the question is
about the structures, and it needs sizes the real chromosomes cannot supply.
They are drawn to match the real length distribution, which is what determines
how much a bounding structure can prune.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from .db import Connection

# Geometric series: each step is 4x the last, so a handful of points span three
# orders of magnitude and any power law shows up as a straight line on a log
# axis.
DEFAULT_SIZES = (1_000, 4_000, 16_000, 64_000, 256_000)

# Human chromosome 1 is about 250 Mb; using a realistic span keeps interval
# density -- and therefore how much a bin or bounding box can prune -- realistic.
COORDINATE_SPAN = 250_000_000


@dataclass(frozen=True)
class Fit:
    """A least-squares fit of time against one candidate growth law."""

    model: str
    slope: float
    intercept: float
    r_squared: float

    def as_dict(self) -> dict[str, float | str]:
        return {
            "model": self.model,
            "slope": round(self.slope, 6),
            "intercept": round(self.intercept, 6),
            "r_squared": round(self.r_squared, 4),
        }


def fit(xs: Sequence[float], ys: Sequence[float], model: str) -> Fit:
    """Least-squares fit of y = a*x + b, with the coefficient of determination."""
    mean_x, mean_y = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mean_x) ** 2 for x in xs)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys, strict=True))

    slope = sxy / sxx if sxx else 0.0
    intercept = mean_y - slope * mean_x

    predicted = [slope * x + intercept for x in xs]
    ss_res = sum((y - p) ** 2 for y, p in zip(ys, predicted, strict=True))
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    r_squared = 1 - ss_res / ss_tot if ss_tot else 1.0

    return Fit(model=model, slope=slope, intercept=intercept, r_squared=r_squared)


def classify(sizes: Sequence[int], times: Sequence[float]) -> dict[str, Any]:
    """Decide which growth law the measurements support.

    The classification is by *how much the cost grew* across the size range,
    compared with how much each candidate law predicts it should have. Over
    1k to 256k rows, O(n) predicts a 256-fold rise, O(log n) about 1.8-fold,
    and O(1) none at all — three well-separated predictions.

    This is used in preference to picking the higher R². With five points, a
    least-squares fit will always name a winner, including for a response that
    is really flat and noisy: a lookup rising 1.2-fold across a 256-fold
    increase in n is not evidence of linear growth however well a line happens
    to pass through it. The fits are still reported as supporting detail.
    """
    linear = fit(list(map(float, sizes)), times, "O(n)")
    logarithmic = fit([math.log2(s) for s in sizes], times, "O(log n)")

    observed = (max(times) / min(times)) if min(times) > 0 else float("inf")
    size_ratio = max(sizes) / min(sizes)
    predicted = {
        "O(1)": 1.0,
        "O(log n)": math.log2(max(sizes)) / math.log2(min(sizes)),
        "O(n)": size_ratio,
    }

    # Compare in log space, so being 2x out is judged the same whether the
    # prediction was 1.8 or 256.
    distances = {
        model: abs(math.log(observed) - math.log(value)) for model, value in predicted.items()
    }
    verdict = min(distances, key=lambda m: distances[m])

    # A verdict is only meaningful if the nearest prediction is clearly nearer
    # than the runner-up.
    ordered = sorted(distances.values())
    decisive = len(ordered) < 2 or (ordered[1] - ordered[0]) > 0.25

    return {
        "linear": linear.as_dict(),
        "logarithmic": logarithmic.as_dict(),
        "best_fit": max((linear, logarithmic), key=lambda f: f.r_squared).model,
        "observed_growth": round(observed, 2),
        "predicted_growth": {k: round(v, 2) for k, v in predicted.items()},
        "verdict": verdict if decisive else "indeterminate",
        "decisive": decisive,
        "growth_factor": round(observed, 2),
    }


def synthesise(count: int, seed: int = 20250101) -> list[tuple[int, str, int, int]]:
    """Generate interval rows with a realistic length distribution.

    Gene lengths are heavy-tailed -- most are a few kb, a few span megabases --
    and that skew is precisely what a hierarchical bin or a bounding box has to
    cope with, so a uniform distribution would flatter them.
    """
    rng = random.Random(seed)
    rows = []
    for identifier in range(1, count + 1):
        # Log-normal: median ~25 kb, long right tail.
        length = int(min(rng.lognormvariate(math.log(25_000), 1.2), 2_000_000))
        start = rng.randrange(0, COORDINATE_SPAN - length)
        rows.append((identifier, "chrS", start, start + length))
    return rows


def _create_synthetic(conn: Connection, rows: list[tuple[int, str, int, int]]) -> None:
    conn.execute("DROP TABLE IF EXISTS synthetic")
    conn.execute(
        "CREATE TABLE synthetic ("
        "  id INTEGER PRIMARY KEY, chrom VARCHAR(10) NOT NULL,"
        "  span_start INTEGER NOT NULL, span_end INTEGER NOT NULL)"
    )
    conn.executemany(
        "INSERT INTO synthetic (id, chrom, span_start, span_end) VALUES (?, ?, ?, ?)",
        rows,
    )
    conn.commit()


POINT_LOOKUP = "SELECT id FROM synthetic WHERE span_start = ?"


def measure_point_lookups(
    conn: Connection, sizes: Sequence[int] = DEFAULT_SIZES, probes: int = 200
) -> dict[str, Any]:
    """Time an equality lookup with and without an index, across table sizes.

    A point lookup is the case a B-tree exists for, so this is where the
    O(log n) versus O(n) distinction should be visible if it is visible anywhere.
    """
    unindexed: list[float] = []
    indexed: list[float] = []
    rng = random.Random(7)

    for size in sizes:
        rows = synthesise(size)
        _create_synthetic(conn, rows)
        targets = [rows[rng.randrange(len(rows))][2] for _ in range(probes)]

        conn.execute("DROP INDEX IF EXISTS idx_synthetic_start")
        conn.commit()
        start = time.perf_counter()
        for target in targets:
            conn.query(POINT_LOOKUP, (target,))
        unindexed.append((time.perf_counter() - start) / probes)

        conn.execute("CREATE INDEX idx_synthetic_start ON synthetic(span_start)")
        conn.commit()
        start = time.perf_counter()
        for target in targets:
            conn.query(POINT_LOOKUP, (target,))
        indexed.append((time.perf_counter() - start) / probes)

        print(
            f"  n={size:>8,}  unindexed {unindexed[-1] * 1e6:>8.1f} us"
            f"   indexed {indexed[-1] * 1e6:>7.1f} us"
            f"   ratio {unindexed[-1] / indexed[-1]:>6.1f}x"
        )

    return {
        "sizes": list(sizes),
        "unindexed_us": [round(t * 1e6, 3) for t in unindexed],
        "indexed_us": [round(t * 1e6, 3) for t in indexed],
        "speedup": [round(u / i, 2) for u, i in zip(unindexed, indexed, strict=True)],
        "unindexed_growth": classify(sizes, unindexed),
        "indexed_growth": classify(sizes, indexed),
    }


def measure_interval_strategies(
    conn: Connection,
    sizes: Sequence[int] = DEFAULT_SIZES,
    probes: int = 100,
    width: int = 1_000_000,
) -> dict[str, Any]:
    """Find where each interval structure wins, as the table grows.

    At small n the naive range scan is fastest because it has no structure to
    traverse. The question is the crossover: how many intervals before a
    hierarchical bin or an R-tree earns its overhead.
    """
    from . import intervals

    # The same inclusive convention as intervals.py. Over a 250 Mb span a random
    # window practically never lands on an interval boundary, so this does not
    # move the timings -- but two modules teaching different conventions is how
    # the wrong one gets copied.
    rng = random.Random(11)
    timings: dict[str, list[float]] = {"btree": [], "binning": [], "rtree": []}

    for size in sizes:
        rows = synthesise(size)
        _create_synthetic(conn, rows)

        conn.execute("DROP INDEX IF EXISTS idx_synthetic_locus")
        conn.execute(
            "CREATE INDEX idx_synthetic_locus ON synthetic(chrom, span_start, span_end)"
        )

        conn.execute("DROP TABLE IF EXISTS synthetic_bin")
        conn.execute(
            "CREATE TABLE synthetic_bin (id INTEGER PRIMARY KEY, bin INTEGER NOT NULL,"
            " span_start INTEGER NOT NULL, span_end INTEGER NOT NULL)"
        )
        conn.executemany(
            "INSERT INTO synthetic_bin (id, bin, span_start, span_end) VALUES (?,?,?,?)",
            [
                (identifier, intervals.assign_bin(start, end), start, end)
                for identifier, _, start, end in rows
            ],
        )
        conn.execute("CREATE INDEX idx_synthetic_bin ON synthetic_bin(bin)")

        conn.execute("DROP TABLE IF EXISTS synthetic_rtree")
        conn.execute(
            "CREATE VIRTUAL TABLE synthetic_rtree USING rtree(id, min_start, max_end)"
        )
        conn.executemany(
            "INSERT INTO synthetic_rtree (id, min_start, max_end) VALUES (?, ?, ?)",
            [(identifier, start, end) for identifier, _, start, end in rows],
        )
        conn.commit()

        queries = [
            (s, s + width)
            for s in (rng.randrange(0, COORDINATE_SPAN - width) for _ in range(probes))
        ]

        plans = {
            "btree": (
                "SELECT id FROM synthetic WHERE chrom='chrS'"
                " AND span_start <= ? AND span_end >= ?",
                lambda s, e: (e, s),
            ),
            "rtree": (
                "SELECT id FROM synthetic_rtree WHERE min_start <= ? AND max_end >= ?",
                lambda s, e: (e, s),
            ),
        }

        for name, (sql, bind) in plans.items():
            start_time = time.perf_counter()
            for s, e in queries:
                conn.query(sql, bind(s, e))
            timings[name].append((time.perf_counter() - start_time) / probes)

        start_time = time.perf_counter()
        for s, e in queries:
            bins = intervals.overlapping_bins(s, e)
            placeholders = ", ".join("?" for _ in bins)
            conn.query(
                f"SELECT id FROM synthetic_bin WHERE bin IN ({placeholders})"
                " AND span_start <= ? AND span_end >= ?",
                (*bins, e, s),
            )
        timings["binning"].append((time.perf_counter() - start_time) / probes)

        fastest = min(timings, key=lambda k: timings[k][-1])
        print(
            f"  n={size:>8,}  "
            + "  ".join(
                f"{k} {timings[k][-1] * 1e6:>8.1f}us" for k in ("btree", "binning", "rtree")
            )
            + f"   fastest: {fastest}"
        )

    winners = [min(timings, key=lambda k: timings[k][index]) for index in range(len(sizes))]
    crossover = next((sizes[i] for i, w in enumerate(winners) if w != winners[0]), None)

    return {
        "sizes": list(sizes),
        "per_query_us": {k: [round(t * 1e6, 2) for t in v] for k, v in timings.items()},
        "winner_by_size": dict(zip(sizes, winners, strict=True)),
        "crossover_n": crossover,
        "growth": {k: classify(sizes, v) for k, v in timings.items()},
    }


def render(point: dict[str, Any], interval: dict[str, Any]) -> str:
    """Markdown tables for the README."""
    lines = ["### Point lookup, indexed versus not", ""]
    lines.append("| n | Unindexed | Indexed | Speed-up |")
    lines.append("| ---: | ---: | ---: | ---: |")
    for size, u, i, s in zip(
        point["sizes"],
        point["unindexed_us"],
        point["indexed_us"],
        point["speedup"],
        strict=True,
    ):
        lines.append(f"| {size:,} | {u:,.1f} µs | {i:,.1f} µs | {s:.1f}× |")

    lines += ["", "### Interval overlap, three structures", ""]
    lines.append("| n | B-tree scan | UCSC binning | R*Tree | Fastest |")
    lines.append("| ---: | ---: | ---: | ---: | --- |")
    per = interval["per_query_us"]
    for index, size in enumerate(interval["sizes"]):
        fastest = min(per, key=lambda k: per[k][index])
        lines.append(
            f"| {size:,} | {per['btree'][index]:,.1f} µs | {per['binning'][index]:,.1f} µs"
            f" | {per['rtree'][index]:,.1f} µs | {fastest} |"
        )
    return "\n".join(lines) + "\n"
