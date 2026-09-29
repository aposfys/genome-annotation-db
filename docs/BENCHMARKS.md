# Benchmarks

Every number on this page is taken from the committed files in [`results/`](../results).
Wall-clock times come from one run on one machine and will differ on yours. VM step counts
and query plans are deterministic and reproduce exactly.

## Indexes are not a blanket win

`genomedb benchmark` drops every secondary index, times all seven queries, recreates the
indexes and times them again. Same data, same queries, only the indexes change.

Two figures are reported. Wall-clock time is what a user feels, but it depends on the
machine. **VM steps, the virtual-machine instructions SQLite executes, measure the work the
query actually does and are identical on any hardware**, so the work ratio is the number
that reproduces.

| Query | Rows | Wall time | Speed-up | VM steps (without → with) | Work ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Q1, transcripts of a named gene | 2 | 2.31 → 0.05 ms | 44× | 59,746 → **61** | **979×** |
| Q2, most heavily transcribed genes | 20 | 11.6 → 6.7 ms | 1.7× | 281,488 → 230,199 | 1.2× |
| Q3, transcripts with the most exons | 30 | 96.1 → 93.5 ms | 1.0× | 2,810,638 → 2,810,638 | **1.0×** |
| Q4, genes carrying ≥ 20 GO terms | 219 | 7.5 → 7.2 ms | 1.0× | 290,146 → 290,146 | **1.0×** |
| Q5, annotation depth per GO namespace | 3 | 9.0 → 8.3 ms | 1.1× | 258,238 → 176,553 | 1.5× |
| Q6, structure of genes carrying a GO term | 13 | 9.5 → 1.8 ms | 5.1× | 176,706 → 79,462 | 2.2× |
| Q7, exon reuse across transcripts | 702 | 129.4 → 114.2 ms | 1.1× | 3,394,126 → 3,394,126 | **1.0×** |

**Q1 does 979× less work with the index, and runs 44× faster by the clock.** The wall-clock
figure includes fixed per-query overhead (parsing, planning, returning rows) that no index
can remove. The work ratio shows the index eliminating 99.9% of the instructions. The
wall-clock ratio also moves between runs and machines, which is why the work ratio is the
one to quote.

**Q3, Q4 and Q7 execute the same number of instructions with and without indexes, and their
plans do not change.** These queries aggregate over every row, so there is nothing for an
index to skip. Whatever separates their wall-clock times, it is not the index.

<details>
<summary>Q1's query plan, before and after</summary>

```
without: SCAN t | SEARCH g USING INDEX sqlite_autoindex_gene_1 (gene_id=?) | USE TEMP B-TREE FOR ORDER BY
with:    SEARCH g USING INDEX idx_gene_name (gene_name=?) | SEARCH t USING INDEX idx_transcript_gene (gene_id=?)
```
</details>

## Interval queries, three structures

The defining query in genome annotation is *what overlaps this region?* Two intervals
overlap when `a_start <= b_end AND a_end >= b_start`, and that conjunction is what makes it
hard to index. A B-tree on `(chrom, start)` can seek to the first candidate but cannot bound
the *end*, because a gene starting far to the left may still reach into the window. The
predicate is two-dimensional and a B-tree is not.

Three answers are implemented and compared on identical data. A plain **B-tree** range
scan, the **UCSC binning scheme** (Kent et al. 2002, pure SQL, which is what the genome
browsers use), and SQLite's **R\*Tree** module. All three must return the same genes. A test
asserts it on small coordinates and on coordinates above 2^24, and the nightly Pipeline
checks it on random windows over the real data.

### What the real database shows

`genomedb intervals --windows 300` times 300 random 1 Mb windows, three times each, against
the 769-gene database ([`results/intervals.json`](../results/intervals.json)).

| Strategy | µs per query | Rows | Plan |
| --- | ---: | ---: | --- |
| B-tree | **21.8** | 7,386 | `SEARCH gene USING INDEX idx_gene_locus (chrom=? AND gene_start<?)` |
| UCSC binning | 30.6 | 7,386 | `SEARCH b USING INDEX idx_gene_bin (chrom=? AND bin=?)`, then the gene by key |
| R\*Tree | 363.6 | 7,386 | `SEARCH m USING INDEX idx_rtree_map_chrom (chrom=?)`, then `SCAN r VIRTUAL TABLE INDEX 1:` |

The B-tree is fastest and the R\*Tree is last, 17× slower. The plan says why. An R\*Tree
stores only an integer key and the box, so the chromosome and gene identifier live in a
mapping table that has to be joined. Given that join, SQLite chose to start from the
chromosome index on the mapping table and then fetch each gene's box from the R\*Tree by row
id. `INDEX 1` is that row-id lookup. The R\*Tree's spatial search is never used, so every
gene on the chromosome is visited for every window.

The 17× is therefore a measurement of that join order, not of the R\*Tree as a structure.
Writing the query so the R\*Tree drives the join (for example with `CROSS JOIN` or `+m.chrom`)
is the obvious next measurement, and I have not committed it. What the result does show is
that a structure which wins in a micro-benchmark can lose once it sits inside a schema, and
that the plan, not the timing, is what tells you why.

The timings above were measured when the R\*Tree still stored 32-bit float boxes (see
[validation](#validation-against-bedtools)). With the integer table the same windows return
the same 7,386 rows through the same plan.

### A portability trap worth knowing

**MySQL/InnoDB creates an index for every foreign key if none exists. SQLite creates nothing
beyond primary keys.** So a schema with no explicit indexes is already indexed on its join
columns under MySQL, and does full scans for the identical joins under SQLite. A design that
performs acceptably on one engine can be unusable on the other, and nothing in the DDL hints
at it. Every index this schema relies on is therefore declared explicitly in
[`sql/indexes.sql`](../sql/indexes.sql). The benchmark and interval subcommands themselves
use SQLite-only DDL, and CI runs everything on SQLite only.

## Validation against bedtools

The three interval strategies are checked against each other, which catches a mistake in any
one of them but **not a mistake they share**. A coordinate-convention error is exactly that
kind. [bedtools](https://bedtools.readthedocs.io/) is an entirely independent
implementation. Agreeing with it is evidence, agreeing with yourself is not.

| Windows | B-tree agreement with bedtools |
| --- | --- |
| 200 random | **200 / 200** |
| 400 placed exactly on gene boundaries | **400 / 400** |

**Setting this up found a real off-by-one.** Ensembl coordinates are **1-based inclusive**,
so `DEFB125` spans 87,250–97,094, which is 9,845 bases. The overlap predicate was written
with strict `<` and `>`, the *half-open* rule. A gene ending exactly where a window starts
shares one base with it and does overlap, and every strategy was silently excluding those.

The fix is `<=` and `>=`, plus an explicit conversion when exporting to BED, which really is
0-based half-open. A 1-based `[87250, 97094]` becomes BED `[87249, 97094)`.

Two details make the check worth trusting.

- **The boundary windows are deliberate.** Random windows essentially never land on a gene
  edge, so they cannot detect this class of bug. The 200 random windows agreed even
  *before* the fix.
- **The check can fail.** Each boundary gene gets four windows, two that touch its edge and
  two that stop one base short. The strict comparison misses every touching window, so it
  fails on 200 of the 400.

**bedtools has only checked the B-tree strategy**, and the other two were compared with it
only on random windows. That gap hid a second boundary bug. The R\*Tree stored its boxes in
SQLite's default `rtree` module, which uses 32-bit floats and so holds integers exactly only
up to 2^24 = 16.7 Mb. 636 of the 769 genes lie beyond that, and their boxes were rounded
outward. On 3,076 windows placed one base either side of every gene boundary the R\*Tree
disagreed with the B-tree on 2,141, while agreeing on every random window. It now uses
`rtree_i32` and disagrees on none, and a test places a gene above 2^24.

## Empirical complexity, on real coordinates

`genomedb scaling` measures how cost grows across a geometric series of table sizes, using
**real gene coordinates for every gene in the human genome, 78,733 of them**, subsampled to
each size ([`results/scaling.json`](../results/scaling.json)). Chromosomes 20 and 21 hold
only 769 genes between them. The whole genome supplies two orders of magnitude more, with
the clustering and heavy-tailed length distribution that decide how much a bounding
structure can actually prune.

Growth is classified by comparing observed growth against what each law *predicts* over the
range, rather than by picking whichever least-squares fit scores higher. With five points a
fit always names a winner, including in noise.

### Point lookup

| n genes | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 30.4 µs | 7.9 µs | 3.8× |
| 4,000 | 105.6 µs | 7.0 µs | 15.0× |
| 16,000 | 389.4 µs | 7.5 µs | 51.6× |
| 64,000 | 1,520.4 µs | 7.5 µs | 203.7× |
| 78,733 | 2,336.0 µs | 8.7 µs | **270.0×** |

- **Unindexed, O(n).** It grew 76.8× across a 78.7-fold increase in rows, and a linear fit
  gives R² = 0.979.
- **Indexed, indeterminate.** It grew 1.23× over the same range. O(log n) predicts 1.63× and
  O(1) predicts 1.0, and 1.23 is too close to either to call, so the classifier says so
  rather than picking one. At a few microseconds per lookup, run-to-run noise is the same
  size as the effect.

### Interval overlap, three structures

| n genes | B-tree scan | UCSC binning | R\*Tree |
| ---: | ---: | ---: | ---: |
| 1,000 | 30.3 µs | 16.4 µs | **8.8 µs** |
| 4,000 | 95.8 µs | 17.6 µs | **9.8 µs** |
| 16,000 | 364.7 µs | 24.1 µs | **11.6 µs** |
| 64,000 | 1,383.4 µs | 36.6 µs | **20.1 µs** |
| 78,733 | 1,615.1 µs | 40.6 µs | **22.0 µs** |

The B-tree grows 53-fold across the range, and binning and the R\*Tree grow about 2.5-fold.
At 78,733 genes the R\*Tree answers 73× faster than the B-tree, and the gap widens with n.

Three things limit how far this carries.

- **The B-tree is at its worst here.** The chromosomes are laid end to end on one 3.1 Gb
  axis with no chromosome predicate, so the B-tree can bound only the start and reads every
  gene that starts before the window. The real database filters by chromosome first.
- **These are bare tables.** The R\*Tree query here reads the R\*Tree alone, with no join, so
  it does not face the join-order problem above.
- **The R\*Tree here is the float32 module and its results are not checked.** On a 3.1 Gb
  axis its boxes are approximate, and the scaling run times the three strategies without
  comparing what they return. `rtree_i32` cannot hold coordinates beyond 2^31, so the fix
  used for the real database does not carry over directly.

### A ceiling the classic scheme has

Running this genome-wide exposed a real property of UCSC binning. **The classic five-level
scheme reaches only 2²⁹ = 512 Mb.** The human genome is 3.1 Gb, so laying the chromosomes end
to end onto one axis overflows it. That is why UCSC applies the scheme per chromosome. The
longest human chromosome is 249 Mb, so one always fits and a genome never does.

The repository implements both the classic scheme and UCSC's extended six-level form
(2³² = 4.3 Gb), and picks the narrowest that holds the data. The scheme is chosen **once per
dataset, not per interval**. A small feature binned under the classic offsets is invisible to
a query using extended ones, because the two numberings do not correspond. A test asserts
that failure mode directly.
