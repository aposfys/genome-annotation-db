# Benchmarks

## Indexes are not a blanket win

`genomedb benchmark` drops every secondary index, times all seven queries, recreates the
indexes and times them again. Same data, same queries, only the indexes change.

Two figures are reported. Wall-clock time is what a user feels, but it depends on the
machine. **VM steps — the virtual-machine instructions SQLite executes — measure the work
the query actually does, and are identical on any hardware**, so the work ratio is the
number that reproduces.

| Query | Rows | Wall time | Speed-up | VM steps (without → with) | Work ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Q1 — Transcripts of a named gene | 2 | 1.2 → 0.02 ms | 46× | 59,746 → **61** | **979×** |
| Q2 — Most heavily transcribed genes | 20 | 6.9 → 3.4 ms | 2.0× | 281,488 → 230,199 | 1.2× |
| Q3 — Transcripts with the most exons | 30 | 50.0 → 47.6 ms | 1.0× | 2,810,638 → 2,810,638 | **1.0×** |
| Q4 — Genes carrying ≥ 20 GO terms | 219 | 3.9 → 3.8 ms | 1.0× | 290,146 → 290,146 | **1.0×** |
| Q5 — Annotation depth per GO namespace | 3 | 4.5 → 4.4 ms | 1.0× | 258,238 → 176,553 | 1.5× |
| Q6 — Structure of genes carrying a GO term | 13 | 4.9 → 1.0 ms | 5.1× | 176,706 → 79,462 | 2.2× |
| Q7 — Exon reuse across transcripts | 702 | 66.6 → 58.5 ms | 1.1× | 3,394,126 → 3,394,126 | **1.0×** |

**Q1's real gain is 979×, not 46×.** The wall-clock figure is floored by fixed per-query
overhead — parsing, planning, returning rows — which the index cannot remove. Strip that
away and the index eliminates 99.9% of the work.

**Q3, Q4 and Q7 execute a byte-identical number of instructions with and without indexes.**
That is not "we measured no difference"; it is proof the plan did not change. These queries
aggregate over every row, so there is nothing for an index to skip, and Q4 is marginally
*slower* with them — the honest shape of that trade-off.

<details>
<summary>Q1's query plan, before and after</summary>

```
without: SCAN t | SEARCH g USING INDEX sqlite_autoindex_gene_1 (gene_id=?) | USE TEMP B-TREE FOR ORDER BY
with:    SEARCH g USING INDEX idx_gene_name (gene_name=?) | SEARCH t USING INDEX idx_transcript_gene (gene_id=?)
```
</details>

## Interval queries: three structures

The defining query in genome annotation is *what overlaps this region?* Two intervals
overlap when `a_start < b_end AND a_end > b_start`, and that conjunction is what makes it
hard to index: a B-tree on `(chrom, start)` can seek to the first candidate but cannot bound
the *end*, because a gene starting far to the left may still reach into the window. The
predicate is two-dimensional; a B-tree is not.

Three answers are implemented and compared on identical data — a plain **B-tree** range
scan, the **UCSC binning scheme** (Kent et al. 2002, pure SQL, which is what the genome
browsers use), and SQLite's **R\*Tree** module. All three return byte-identical results; a
test asserts it, because a faster structure that returns different rows is not a faster
structure.

### The result that matters

Run the comparison against the *real* database and the R\*Tree comes last, 17× slower than
the naive B-tree:

| On the real 769-gene database | µs per query |
| --- | ---: |
| R\*Tree, bare search | **12.8** |
| B-tree, single-table scan | 21.8 |
| R\*Tree **plus the two joins back to `gene`** | **382.5** |

**The R\*Tree search is the fastest of the three. The joins around it cost 30× more than the
search does.** An R\*Tree stores integer keys, so recovering the gene identifier and
filtering by chromosome needs a mapping table and a join back to `gene`.

Two things follow, and neither is visible from a single benchmark:

- **Micro-benchmarks of a data structure can invert once it is embedded in a schema.**
- **The crossover depends on the join cost, not just on n.** Binning wins in practice at
  this scale precisely because it keeps the coordinates and the key in one indexed table, so
  it needs no mapping hop. That is very likely why UCSC chose it.

### A portability trap worth knowing

**MySQL/InnoDB silently creates an index for every foreign key. SQLite creates nothing
beyond primary keys.** So a schema with no explicit indexes is already indexed on its join
columns under MySQL, and does full scans for the identical joins under SQLite. A design that
performs acceptably on one engine can be unusable on the other, and nothing in the DDL hints
at it. Every index this schema relies on is therefore declared explicitly in
[`sql/indexes.sql`](../sql/indexes.sql).

## Validation against bedtools

The three interval strategies are checked against each other, which catches a mistake in any
one of them but **not a mistake they share**. A coordinate-convention error is exactly that
kind. [bedtools](https://bedtools.readthedocs.io/) is an entirely independent
implementation; agreeing with it is evidence, agreeing with yourself is not.

| Windows | Agreement |
| --- | --- |
| 200 random | **200 / 200** |
| 400 placed exactly on gene boundaries | **400 / 400** |

**Setting this up found a real off-by-one.** Ensembl coordinates are **1-based inclusive** —
`DEFB125` spans 87,250–97,094, which is 9,845 bases. The overlap predicate was written with
strict `<` and `>`, the *half-open* rule. A gene ending exactly where a window starts shares
one base with it and does overlap, and every strategy was silently excluding those.

The fix is `<=` and `>=`, plus an explicit conversion when exporting to BED, which really is
0-based half-open: a 1-based `[87250, 97094]` becomes BED `[87249, 97094)`.

Two details make the check worth trusting:

- **The boundary windows are deliberate.** Random windows essentially never land on a gene
  edge, so they cannot detect this class of bug — the 200 random windows agreed even
  *before* the fix.
- **The check was verified to be capable of failing.** Reintroducing the strict comparison
  makes it fail on **200 of 400** boundary windows. A validation that has never been seen to
  fail is not yet a validation.

## Empirical complexity, on real coordinates

`genomedb scaling` measures how cost grows across a geometric series of table sizes, using
**real gene coordinates — every gene in the human genome, 78,733 of them**, subsampled to
each size. Chromosomes 20 and 21 hold only 769 genes between them, which is why an earlier
version generated intervals instead; the whole genome supplies two orders of magnitude more,
with the clustering and heavy-tailed length distribution that decide how much a bounding
structure can actually prune.

Growth is classified by comparing observed growth against what each law *predicts* over the
range, rather than by picking whichever least-squares fit scores higher — with five points a
fit always names a winner, including in noise.

### Point lookup

| n genes | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 17.7 µs | 3.6 µs | 4.9× |
| 4,000 | 54.0 µs | 3.7 µs | 14.4× |
| 16,000 | 195.6 µs | 3.8 µs | 51.5× |
| 64,000 | 781.7 µs | 4.2 µs | 185.2× |
| 78,733 | 1,179.4 µs | 3.9 µs | **305.4×** |

- **Unindexed: O(n).** Linear fit at R² = 1.000.
- **Indexed: indeterminate — and that is the honest answer.** It grew 1.16× across a 79-fold
  increase in rows. O(log n) predicts 1.58× over that range and O(1) predicts 1.0; 1.16 sits
  between them, too close to either to call. The classifier says so rather than picking one.

### Interval overlap, three structures

| n genes | B-tree scan | UCSC binning | R\*Tree |
| ---: | ---: | ---: | ---: |
| 1,000 | 14.1 µs | 8.2 µs | **4.4 µs** |
| 4,000 | 46.9 µs | 9.0 µs | **5.0 µs** |
| 16,000 | 190.4 µs | 10.6 µs | **5.9 µs** |
| 64,000 | 742.9 µs | 19.7 µs | **10.2 µs** |
| 78,733 | 885.4 µs | 20.7 µs | **11.1 µs** |

The B-tree grows 63-fold across the range; binning and the R\*Tree grow 2.5-fold. **At every
gene in the genome the R\*Tree answers 80× faster than the B-tree**, and the gap widens with
n.

### A ceiling the classic scheme has

Running this genome-wide exposed a real property of UCSC binning: **the classic five-level
scheme reaches only 2²⁹ = 512 Mb.** The human genome is 3.1 Gb, so laying the chromosomes
end to end onto one axis overflows it. That is exactly why UCSC applies the scheme per
chromosome — the longest human chromosome is 249 Mb, so one always fits and a genome never
does.

The repository implements both the classic scheme and UCSC's extended six-level form
(2³² = 4.3 Gb), and picks the narrowest that holds the data. The scheme is chosen **once per
dataset, not per interval**: a small feature binned under the classic offsets is invisible to
a query using extended ones, because the two numberings do not correspond. A test asserts
that failure mode directly.
