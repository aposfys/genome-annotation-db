"""Genomic interval queries, and three data structures for answering them.

The defining query in genome annotation is *what overlaps this region?* Two
intervals [a1, a2] and [b1, b2] overlap when ``a1 <= b2 AND a2 >= b1``, and that
conjunction is what makes the query hard to index.

A B-tree on ``(chrom, gene_start)`` can seek to the first candidate by start
position, but it cannot bound the *end*: a gene beginning far to the left may
still extend into the query window, so the engine has to keep reading. The
predicate is two-dimensional and a B-tree is a one-dimensional structure.

Three answers are implemented here so they can be compared on the same data:

* **btree** -- the naive baseline, a plain composite index and a range scan.
* **binning** -- the UCSC scheme (Kent et al. 2002). Each interval is assigned
  to the smallest of a hierarchy of fixed bins that fully contains it; a query
  then needs only the handful of bins its window can intersect. Pure SQL, no
  extension required, which is why the genome browsers adopted it.
* **rtree** -- SQLite's R*Tree module, a genuine multidimensional index that
  stores the interval as a one-dimensional bounding box.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from typing import Any

from .db import Connection

# --- UCSC binning ------------------------------------------------------------
#
# The scheme layers coarse-to-fine bins over the coordinate space. Level 0 is a
# single bin covering 512 Mb; each subsequent level divides by 8. An interval is
# stored in the finest level whose bin contains it whole, so a large feature
# lands in a coarse bin and a small one in a fine bin.
#
# Kent, W. J. et al. (2002) The human genome browser at UCSC.
# Genome Research 12, 996-1006.

BIN_OFFSETS = (512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0)
# The finest bin spans 2^17 = 128 kb; each coarser level multiplies by 8.
BIN_FIRST_SHIFT = 17
BIN_NEXT_SHIFT = 3


def assign_bin(start: int, end: int) -> int:
    """The finest UCSC bin that wholly contains a half-open interval.

    Args:
        start: 0-based inclusive start.
        end: exclusive end.
    """
    start_bin, end_bin = start, end - 1
    start_bin >>= BIN_FIRST_SHIFT
    end_bin >>= BIN_FIRST_SHIFT

    for offset in BIN_OFFSETS:
        if start_bin == end_bin:
            return offset + start_bin
        start_bin >>= BIN_NEXT_SHIFT
        end_bin >>= BIN_NEXT_SHIFT

    raise ValueError(f"interval {start}-{end} exceeds the binning scheme's range")


def overlapping_bins(start: int, end: int) -> list[int]:
    """Every bin an interval could intersect.

    A query must look in its own fine bins *and* in every coarser bin above
    them, because a large feature containing the window is stored higher up.
    """
    start_bin, end_bin = start, end - 1
    start_bin >>= BIN_FIRST_SHIFT
    end_bin >>= BIN_FIRST_SHIFT

    bins: list[int] = []
    for offset in BIN_OFFSETS:
        bins.extend(range(offset + start_bin, offset + end_bin + 1))
        start_bin >>= BIN_NEXT_SHIFT
        end_bin >>= BIN_NEXT_SHIFT
    return bins


# --- the three strategies ----------------------------------------------------


@dataclass(frozen=True)
class Strategy:
    """One way of answering an overlap query."""

    name: str
    description: str
    requires_sqlite: bool = False


STRATEGIES: tuple[Strategy, ...] = (
    Strategy("btree", "Composite B-tree on (chrom, gene_start), range scan"),
    Strategy("binning", "UCSC hierarchical binning scheme, pure SQL"),
    Strategy("rtree", "SQLite R*Tree multidimensional index", requires_sqlite=True),
)

BY_NAME = {strategy.name: strategy for strategy in STRATEGIES}


def build_btree(conn: Connection) -> None:
    """The baseline: a composite index and nothing else."""
    conn.execute("DROP INDEX IF EXISTS idx_gene_locus")
    conn.execute("CREATE INDEX idx_gene_locus ON gene(chrom, gene_start, gene_end)")
    conn.commit()


def build_binning(conn: Connection) -> None:
    """Materialise a bin number per gene, indexed by (chrom, bin)."""
    conn.execute("DROP TABLE IF EXISTS gene_bin")
    conn.execute(
        "CREATE TABLE gene_bin ("
        "  gene_id VARCHAR(30) PRIMARY KEY,"
        "  chrom   VARCHAR(10) NOT NULL,"
        "  bin     INTEGER NOT NULL,"
        "  gene_start INTEGER NOT NULL,"
        "  gene_end   INTEGER NOT NULL"
        ")"
    )
    _, rows = conn.query("SELECT gene_id, chrom, gene_start, gene_end FROM gene")
    conn.executemany(
        "INSERT INTO gene_bin (gene_id, chrom, bin, gene_start, gene_end)"
        " VALUES (?, ?, ?, ?, ?)",
        [
            (gene_id, chrom, assign_bin(start, end), start, end)
            for gene_id, chrom, start, end in rows
        ],
    )
    conn.execute("CREATE INDEX idx_gene_bin ON gene_bin(chrom, bin)")
    conn.commit()


def build_rtree(conn: Connection) -> None:
    """Store each gene's span in an R*Tree virtual table."""
    if not conn.backend.is_sqlite:
        raise RuntimeError("the R*Tree strategy needs SQLite")

    conn.execute("DROP TABLE IF EXISTS gene_rtree")
    conn.execute("DROP TABLE IF EXISTS gene_rtree_map")
    # R*Tree keys are integers, so a mapping table carries the gene identifier
    # and the chromosome, which the tree itself cannot store.
    conn.execute("CREATE VIRTUAL TABLE gene_rtree USING rtree(id, min_start, max_end)")
    conn.execute(
        "CREATE TABLE gene_rtree_map ("
        "  id INTEGER PRIMARY KEY, gene_id VARCHAR(30) NOT NULL, chrom VARCHAR(10) NOT NULL"
        ")"
    )

    _, rows = conn.query(
        "SELECT gene_id, chrom, gene_start, gene_end FROM gene ORDER BY gene_id"
    )
    conn.executemany(
        "INSERT INTO gene_rtree (id, min_start, max_end) VALUES (?, ?, ?)",
        [(index, start, end) for index, (_, _, start, end) in enumerate(rows, 1)],
    )
    conn.executemany(
        "INSERT INTO gene_rtree_map (id, gene_id, chrom) VALUES (?, ?, ?)",
        [(index, gene_id, chrom) for index, (gene_id, chrom, _, _) in enumerate(rows, 1)],
    )
    conn.execute("CREATE INDEX idx_rtree_map_chrom ON gene_rtree_map(chrom)")
    conn.commit()


BUILDERS = {"btree": build_btree, "binning": build_binning, "rtree": build_rtree}


def build(conn: Connection, strategy: str) -> None:
    if strategy not in BUILDERS:
        raise KeyError(f"Unknown strategy {strategy!r}; expected one of {sorted(BUILDERS)}")
    BUILDERS[strategy](conn)


# --- querying ----------------------------------------------------------------

# Half-open convention throughout: [start, end). Two intervals overlap when each
# begins before the other ends.
BTREE_SQL = """
SELECT gene_id, gene_name, gene_start, gene_end
FROM gene
WHERE chrom = ? AND gene_start < ? AND gene_end > ?
ORDER BY gene_start
"""

RTREE_SQL = """
SELECT m.gene_id, g.gene_name, r.min_start, r.max_end
FROM gene_rtree AS r
JOIN gene_rtree_map AS m ON m.id = r.id
JOIN gene AS g ON g.gene_id = m.gene_id
WHERE r.min_start < ? AND r.max_end > ? AND m.chrom = ?
ORDER BY r.min_start
"""


def binning_sql(bins: list[int]) -> str:
    placeholders = ", ".join("?" for _ in bins)
    return f"""
SELECT b.gene_id, g.gene_name, b.gene_start, b.gene_end
FROM gene_bin AS b
JOIN gene AS g ON g.gene_id = b.gene_id
WHERE b.chrom = ? AND b.bin IN ({placeholders})
  AND b.gene_start < ? AND b.gene_end > ?
ORDER BY b.gene_start
"""


def query(conn: Connection, strategy: str, chrom: str, start: int, end: int) -> list[tuple]:
    """Genes overlapping ``chrom:start-end`` under one strategy.

    Every strategy returns the same rows in the same order; only the route
    differs. That is what makes the comparison meaningful, and a test asserts it.
    """
    if strategy == "btree":
        _, rows = conn.query(BTREE_SQL, (chrom, end, start))
    elif strategy == "binning":
        bins = overlapping_bins(start, end)
        _, rows = conn.query(binning_sql(bins), (chrom, *bins, end, start))
    elif strategy == "rtree":
        _, rows = conn.query(RTREE_SQL, (end, start, chrom))
    else:
        raise KeyError(f"Unknown strategy {strategy!r}")
    return rows


def plan(conn: Connection, strategy: str, chrom: str, start: int, end: int) -> str:
    """The planner's route for one strategy."""
    keyword = "EXPLAIN QUERY PLAN" if conn.backend.is_sqlite else "EXPLAIN"
    sql: str
    params: tuple[Any, ...]
    if strategy == "btree":
        sql, params = BTREE_SQL, (chrom, end, start)
    elif strategy == "binning":
        bins = overlapping_bins(start, end)
        sql, params = binning_sql(bins), (chrom, *bins, end, start)
    else:
        sql, params = RTREE_SQL, (end, start, chrom)

    try:
        _, rows = conn.query(f"{keyword} {sql}", params)
    except Exception as error:
        return f"(plan unavailable: {error})"
    return " | ".join(" ".join(str(f) for f in row if f not in (None, "")) for row in rows)


def windows(
    conn: Connection, count: int = 200, width: int = 1_000_000, seed: int = 20250101
) -> Iterator[tuple[str, int, int]]:
    """Random query windows drawn from the coordinate range actually present.

    Windows are sampled from the span each chromosome occupies rather than from
    an arbitrary range, so the benchmark measures realistic hit rates instead of
    mostly-empty lookups.
    """
    import random

    _, extents = conn.query(
        "SELECT chrom, MIN(gene_start), MAX(gene_end) FROM gene GROUP BY chrom"
    )
    rng = random.Random(seed)
    for _ in range(count):
        chrom, low, high = extents[rng.randrange(len(extents))]
        if high - low <= width:
            yield chrom, low, high
            continue
        start = rng.randrange(low, high - width)
        yield chrom, start, start + width


def region_stats(conn: Connection) -> dict[str, Any]:
    """Facts about the interval data that explain the benchmark's shape."""
    _, rows = conn.query(
        "SELECT COUNT(*), MIN(gene_end - gene_start), MAX(gene_end - gene_start),"
        " AVG(gene_end - gene_start) FROM gene"
    )
    count, shortest, longest, mean = rows[0]
    return {
        "genes": int(count),
        "shortest_bp": int(shortest),
        "longest_bp": int(longest),
        "mean_bp": round(float(mean), 1),
    }
