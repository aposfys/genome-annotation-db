# Genome Annotation Database: Schema Design and Index Benchmarking

[![CI](https://github.com/aposfys/genome-annotation-db/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/genome-annotation-db/actions/workflows/ci.yml)
[![Pipeline](https://github.com/aposfys/genome-annotation-db/actions/workflows/pipeline.yml/badge.svg)](https://github.com/aposfys/genome-annotation-db/actions/workflows/pipeline.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Ruff](https://img.shields.io/badge/lint-ruff-261230)](https://docs.astral.sh/ruff/)

A normalised relational database over Ensembl gene annotation — genes, their transcripts, the exons those transcripts are assembled from, and Gene Ontology terms — with a validating loader, seven analysis queries, and a benchmark that measures what the indexes are actually worth rather than assuming.

Clone it and you have a working 148,000-row database in under a second. No server, no credentials.

| | |
| --- | --- |
| **Data** | Ensembl BioMart release 113, human chromosomes 20 and 21; whole-genome gene set for the scaling study |
| **Scale** | 769 genes · 9,950 transcripts · 26,970 exons · 97,742 transcript–exon links · 12,596 GO assignments |
| **Engines** | SQLite by default; the same schema and queries run on MySQL 8 |
| **Headline** | Indexes buy **48×** on a point lookup and **nothing** on four of seven queries |
| **Also** | An R\*Tree beats a B-tree 80× on interval search — then loses 17× once the joins around it are counted |
| **Validated** | Overlap results agree with `bedtools` on 600/600 windows, including 400 placed on gene boundaries |

## The result: indexes are not a blanket win

`genomedb benchmark` drops every secondary index, times all seven queries, recreates the indexes and times them again. Same data, same queries, only the indexes change.

Two figures are reported. Wall-clock time is what a user feels, but it depends on the machine. **VM steps — the virtual-machine instructions SQLite executes — measure the work the query actually does, and are identical on any hardware**, so the work ratio is the number that reproduces.

| Query | Rows | Wall time | Speed-up | VM steps (without → with) | Work ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Q1 — Transcripts of a named gene | 2 | 1.2 → 0.02 ms | 46× | 59,746 → **61** | **979×** |
| Q2 — Most heavily transcribed genes | 20 | 6.9 → 3.4 ms | 2.0× | 281,488 → 230,199 | 1.2× |
| Q3 — Transcripts with the most exons | 30 | 50.0 → 47.6 ms | 1.0× | 2,810,638 → 2,810,638 | **1.0×** |
| Q4 — Genes carrying ≥ 20 GO terms | 219 | 3.9 → 3.8 ms | 1.0× | 290,146 → 290,146 | **1.0×** |
| Q5 — Annotation depth per GO namespace | 3 | 4.5 → 4.4 ms | 1.0× | 258,238 → 176,553 | 1.5× |
| Q6 — Structure of genes carrying a GO term | 13 | 4.9 → 1.0 ms | 5.1× | 176,706 → 79,462 | 2.2× |
| Q7 — Exon reuse across transcripts | 702 | 66.6 → 58.5 ms | 1.1× | 3,394,126 → 3,394,126 | **1.0×** |

The VM counts say more than the clock does.

**Q1's real gain is 979×, not 46×.** The wall-clock figure is floored by fixed per-query overhead — parsing, planning, returning rows — which the index cannot remove. Strip that away and the index eliminates 99.9% of the work.

**Q3, Q4 and Q7 execute a byte-identical number of instructions with and without indexes.** That is not "we measured no difference"; it is proof the plan did not change. These queries aggregate over every row, so there is nothing for an index to skip, and Q4 is marginally *slower* with them — the honest shape of that trade-off.

An index earns its keep in proportion to how much of the table it lets you skip. A report that touches everything is not a candidate.

<details>
<summary>Q1's query plan, before and after</summary>

```
without: SCAN t | SEARCH g USING INDEX sqlite_autoindex_gene_1 (gene_id=?) | USE TEMP B-TREE FOR ORDER BY
with:    SEARCH g USING INDEX idx_gene_name (gene_name=?) | SEARCH t USING INDEX idx_transcript_gene (gene_id=?)
```
</details>

## Interval queries: three structures, and a lesson about where cost actually lives

The defining query in genome annotation is *what overlaps this region?* Two intervals overlap when `a_start < b_end AND a_end > b_start`, and that conjunction is what makes it hard to index: a B-tree on `(chrom, start)` can seek to the first candidate but cannot bound the *end*, because a gene starting far to the left may still reach into the window. The predicate is two-dimensional; a B-tree is not.

Three answers are implemented and compared on identical data — a plain **B-tree** range scan, the **UCSC binning scheme** (Kent et al. 2002, pure SQL, which is what the genome browsers use), and SQLite's **R\*Tree** module. All three return byte-identical results; a test asserts it, because a faster structure that returns different rows is not a faster structure.

```bash
genomedb region 20 1000000 2000000 --strategy rtree
genomedb intervals --windows 300
```

### How they scale

Measured on real gene coordinates up to the whole genome — see [Empirical complexity](#empirical-complexity-on-real-coordinates) for the full table. The B-tree grows 63-fold across the range while binning and the R\*Tree grow 2.5-fold, leaving the R\*Tree **80× faster at 78,733 genes**.

### The result that matters

Run the same comparison against the *real* database and the R\*Tree comes last, 17× slower than the naive B-tree. That contradiction is the most instructive thing in this project, so it was measured rather than explained away:

| On the real 769-gene database | µs per query |
| --- | ---: |
| R\*Tree, bare search | **12.8** |
| B-tree, single-table scan | 21.8 |
| R\*Tree **plus the two joins back to `gene`** | **382.5** |

**The R\*Tree search is the fastest of the three. The joins around it cost 30× more than the search does.**

An R\*Tree stores integer keys, so recovering the gene identifier and filtering by chromosome needs a mapping table and a join back to `gene` — and at this scale that dominates completely. The index was never the bottleneck; the normalisation around it was.

Two things follow, and neither is visible from a single benchmark:

- **Micro-benchmarks of a data structure can invert once it is embedded in a schema.** The structure that wins in isolation lost by 17× in place.
- **The crossover depends on the join cost, not just on n.** Binning wins in practice at this scale precisely because it keeps the coordinates and the key in one indexed table, so it needs no mapping hop. That is very likely why UCSC chose it.

### A portability trap worth knowing

**MySQL/InnoDB silently creates an index for every foreign key. SQLite creates nothing beyond primary keys.**

So a schema with no explicit indexes is already indexed on its join columns under MySQL, and does full scans for the identical joins under SQLite. A design that performs acceptably on one engine can be unusable on the other, and nothing in the DDL hints at it. Every index this schema relies on is therefore declared explicitly in [`sql/indexes.sql`](sql/indexes.sql) rather than left to the engine.

## Validation against bedtools — and the bug it found

The three interval strategies are checked against each other, which catches a mistake in any one of them but **not a mistake they share**. A coordinate-convention error is exactly that kind: all three read the same columns and apply the same comparison, so all three would agree and all three would be wrong.

[bedtools](https://bedtools.readthedocs.io/) is the standard toolkit for genomic interval arithmetic and an entirely independent implementation. Agreeing with it is evidence; agreeing with yourself is not.

```bash
make validate
```

| Windows | Agreement |
| --- | --- |
| 200 random | **200 / 200** |
| 400 placed exactly on gene boundaries | **400 / 400** |

**Setting this up found a real off-by-one.** Ensembl coordinates are **1-based inclusive** — `DEFB125` spans 87,250–97,094, which is 9,845 bases. The overlap predicate was written with strict `<` and `>`, the *half-open* rule. A gene ending exactly where a window starts shares one base with it and does overlap, and every strategy was silently excluding those.

The fix is `<=` and `>=`, plus an explicit conversion when exporting to BED, which really is 0-based half-open: a 1-based `[87250, 97094]` becomes BED `[87249, 97094)` — the start moves back one, the end does not.

Two details make the check worth trusting:

- **The boundary windows are deliberate.** Random windows essentially never land on a gene edge, so they cannot detect this class of bug — the 200 random windows agreed even *before* the fix. The 400 boundary windows place each query exactly on a gene's first and last base, where an off-by-one has to show.
- **The check was verified to be capable of failing.** Reintroducing the strict comparison makes it fail on **200 of 400** boundary windows. A validation that has never been seen to fail is not yet a validation.

## Normalisation: proved, then priced

`make normalisation` states the functional dependencies, checks Boyce-Codd Normal Form against them mechanically, and confirms the dependencies actually hold in the loaded data.

| Relation | Candidate key | BCNF |
| --- | --- | --- |
| `gene` | (gene_id) | yes |
| `transcript` | (transcript_id) | yes |
| `exon` | (exon_id) | yes |
| `go_term` | (go_id) | yes |
| `transcript_exon` | (transcript_id, exon_id) | yes |
| `gene_go` | (gene_id, go_id) | yes |

BCNF holds when every non-trivial dependency has a superkey on its left-hand side. Here each table's dependencies are keyed on its primary key, which makes the argument short — and `gene_go` is an all-key relation with no non-key attribute to depend on anything.

The interesting case is `transcript_exon`. `exon_rank` depends on the **whole** composite key, not on `exon_id` alone: the same exon can be third in one transcript and first in another. Had rank been stored on `exon`, the schema would have been wrong in a way that only shows up on alternatively spliced genes.

### What that property costs

Counting a transcript's exons means joining the junction table every time. Materialising the count removes the join:

| | 9,950 transcripts |
| --- | ---: |
| Join every time (normalised) | 8.51 ms |
| Materialised column | 1.94 ms |
| | **4.4× faster** |

So normalisation costs about 6.6 ms on this query. What it buys is that the answer cannot be wrong: nothing in the schema can keep a materialised count true, and any write to `transcript_exon` that forgets to update it leaves the two disagreeing — a constraint cannot express that dependency. The 4.4× is the price of that guarantee, stated rather than assumed.

## Empirical complexity, on real coordinates

`genomedb scaling` measures how cost grows across a geometric series of table sizes, using **real gene coordinates — every gene in the human genome, 78,733 of them**, subsampled to each size. Chromosomes 20 and 21 hold only 769 genes between them, which is why an earlier version generated intervals instead; the whole genome supplies two orders of magnitude more, with the clustering and heavy-tailed length distribution that decide how much a bounding structure can actually prune.

Growth is classified by comparing observed growth against what each law *predicts* over the range, rather than by picking whichever least-squares fit scores higher — with five points a fit always names a winner, including in noise.

### Point lookup

| n genes | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 17.7 µs | 3.6 µs | 4.9× |
| 4,000 | 54.0 µs | 3.7 µs | 14.4× |
| 16,000 | 195.6 µs | 3.8 µs | 51.5× |
| 64,000 | 781.7 µs | 4.2 µs | 185.2× |
| 78,733 | 1,179.4 µs | 3.9 µs | **305.4×** |

- **Unindexed: O(n).** Linear fit at R² = 1.000.
- **Indexed: indeterminate — and that is the honest answer.** It grew 1.16× across a 79-fold increase in rows. O(log n) predicts 1.58× over that range and O(1) predicts 1.0; 1.16 sits between them, too close to either to call. The classifier says so rather than picking one.

### Interval overlap, three structures

| n genes | B-tree scan | UCSC binning | R\*Tree |
| ---: | ---: | ---: | ---: |
| 1,000 | 14.1 µs | 8.2 µs | **4.4 µs** |
| 4,000 | 46.9 µs | 9.0 µs | **5.0 µs** |
| 16,000 | 190.4 µs | 10.6 µs | **5.9 µs** |
| 64,000 | 742.9 µs | 19.7 µs | **10.2 µs** |
| 78,733 | 885.4 µs | 20.7 µs | **11.1 µs** |

The B-tree grows 63-fold across the range; binning and the R\*Tree grow 2.5-fold. **At every gene in the genome the R\*Tree answers 80× faster than the B-tree**, and the gap widens with n — which is the whole argument for a structure that can bound both ends of an interval rather than only one.

### A ceiling the classic scheme has

Running this genome-wide exposed a real property of UCSC binning: **the classic five-level scheme reaches only 2²⁹ = 512 Mb.** The human genome is 3.1 Gb, so laying the chromosomes end to end onto one axis overflows it. That is exactly why UCSC applies the scheme per chromosome — the longest human chromosome is 249 Mb, so one always fits and a genome never does.

The repository implements both the classic scheme and UCSC's extended six-level form (2³² = 4.3 Gb), and picks the narrowest that holds the data. The scheme is chosen **once per dataset, not per interval**: a small feature binned under the classic offsets is invisible to a query using extended ones, because the two numberings do not correspond. A test asserts that failure mode directly.

## The schema

```
gene ──< transcript ──< transcript_exon >── exon
  │
  └──< gene_go >── go_term
```

Six tables, third normal form, with two junction tables carrying the many-to-many relationships.

**The exon model is the design decision that matters.** An exon is shared between transcripts rather than duplicated per transcript — 26,970 distinct exons support 97,742 transcript–exon links, so the average exon is used by 3.6 transcripts. Modelling that as a junction table is what makes Q7 answerable at all: *how much of a gene's exon set is constitutive, and how much is alternative?* PLCB4 turns out to have 172 distinct exons of which 98 are shared across two or more transcripts.

Duplicating exons per transcript would have been simpler to load and would have made that question unanswerable without string-matching coordinates.

## Quick start

```bash
pip install -e ".[dev]"

make build       # create the schema and load the data (~0.6 s)
make queries     # run Q1-Q7, write results/
make check       # integrity and consistency checks
make benchmark   # time the queries with and without indexes
make intervals   # B-tree vs UCSC binning vs R*Tree on overlap queries
make scaling     # growth of query cost with table size
make validate    # cross-check overlap results against bedtools
make normalisation  # BCNF check and the cost of normalisation
make test

genomedb fetch   # regenerate the Ensembl exports at the pinned release
```

Ad-hoc lookups:

```bash
genomedb gene DEFB125              # transcripts of a gene
genomedb go GO:0007186             # genes carrying a GO term
genomedb query Q1 --params PLCB4   # any query, re-parameterised
genomedb query Q7 --limit 10 --save
```

Against MySQL instead of SQLite:

```bash
export GENOMEDB_URL='mysql://user:password@host/database'
pip install -e ".[mysql]"
make build
```

The core package has **no third-party dependencies** — SQLite is in the standard library. `mysql-connector-python` is an optional extra, imported only if a MySQL URL is used.

## Data quality

The loader validates as it goes and reports what it refused, because silently dropping malformed rows leaves a database that looks complete and is not:

```
Inserted:
  exon                 26,970      gene                    769
  gene_go              12,596      go_term               3,817
  transcript            9,950      transcript_exon      97,742
Rejected:
  go: unmapped or missing domain      3,483   e.g. GO:0042421 domain='empty'
  go: no accession                    3,352   e.g. gene ENSG00000260861
```

**6,835 of the 69,501 GO rows — just under 10% — are unusable**, carrying an accession with no domain or no accession at all. That is a property of the BioMart export, not a bug, but a loader that discards them without saying so is hiding the fact that a tenth of the annotation never arrived.

Nine integrity checks then run against the loaded database: orphan rows in each junction table, namespaces outside BP/MF/CC, invalid strands, inverted coordinate spans, transcripts extending beyond their gene, transcripts with no exons, and non-contiguous exon ranks. All nine pass.

The foreign-key check is not redundant, incidentally: **SQLite accepts `PRAGMA foreign_keys = OFF` silently**, and a database loaded that way looks fine until something reads it. A test asserts the constraints are actually enforced.

## The discrepancy, explained

An earlier result file for the GO namespace query contains a single row — biological process only, 769 genes, 13,939 annotations. The rebuilt database reports all three namespaces:

| Namespace | Genes | Annotations | Mean per gene |
| --- | ---: | ---: | ---: |
| Biological process | 639 | 5,638 | 8.82 |
| Cellular component | 728 | 3,602 | 4.95 |
| Molecular function | 701 | 3,356 | 4.79 |

Rather than leave that unexplained, the old numbers were reproduced from the shipped export. **Ignoring the GO domain entirely — collapsing every term into one namespace — gives 757 genes and 13,437 annotations**, against the original's 769 and 13,939. That is within 1.6% and 3.6%, and it reproduces the *shape* exactly: one namespace, essentially every gene, and an annotation count larger than the 12,596 distinct gene–term pairs the data actually contains.

So there were two causes, not one:

- **The namespace was collapsed.** Respecting the GO domain gives 12,596 pairs across three namespaces; ignoring it gives 13,437 in one. The original is the second shape.
- **The source export was slightly larger.** The residual — 12 genes and 502 annotations — is what remains after the logic is accounted for, and is consistent with the original having been loaded from a marginally different download.

The second half is why `genomedb fetch` exists. The exports were originally produced by hand through the BioMart web interface, so there was no way to tell whether a disagreement came from the code or the data. The queries are now in [`biomart.py`](src/genomedb/biomart.py), pinned to **Ensembl release 113** through the archive URL, and any export can be regenerated byte-for-byte:

```bash
genomedb fetch                      # all four exports, pinned release
genomedb fetch --export genes_all   # just the whole-genome gene set
```

## Design decisions

- **Credentials come from the environment, never the source.** `$GENOMEDB_URL` or `--db-url`. A host, user and password committed to a file are in version control permanently, and are exactly what an ETL script tends to accumulate.
- **Queries are parameterised.** Q1 and Q6 take their gene symbol and GO term as bound parameters rather than SQL literals, so they are reusable and cannot be broken by a quoting mistake.
- **`CHECK` constraints instead of MySQL's `ENUM`.** Equivalent expressiveness, and it runs on both engines. Constraints also encode the domain rules — strand is ±1, spans are non-inverted, exon ranks start at 1 — so bad data is rejected by the database rather than by convention.
- **No columns the loader cannot fill.** `transcript_name` and `biotype` were declared but never populated, showing up as blank columns in every result. Dropped. `tx_length` was present in the source and unused; added.
- **Every result file has a header row.** Two of the seven previously did not, which makes them awkward to read back programmatically.
- **The benchmark drops indexes before timing them, in that order.** Timing the indexed pass second means it cannot be flattered by a cache the unindexed pass warmed.

## Repository layout

```
sql/
  schema.sql          Six tables, portable across SQLite and MySQL
  indexes.sql         Secondary indexes, separate so they can be measured
  queries/Q1-Q7.sql   The analysis queries, parameterised
src/genomedb/
  db.py         Connection handling and the SQLite/MySQL dialect shim
  load.py       Validating ETL with rejection reporting
  queries.py    Query registry, runner and TSV output
  quality.py    Integrity and consistency checks
  biomart.py    Pinned Ensembl queries, so every export can be regenerated
  intervals.py  UCSC binning (classic and extended), R*Tree, B-tree
  external.py   bedtools cross-validation and BED coordinate conversion
  normalisation.py  Functional dependencies, BCNF check, denormalisation cost
  scaling.py    Growth measurement and growth-law classification
  benchmark.py  Index timing and query-plan capture
  cli.py        Subcommands: build, query, gene, go, check, benchmark
data/           Ensembl BioMart exports, gzipped (1.7 MB)
results/        Query output, load report, quality report, benchmark
tests/          pytest suite (55 tests)
```

## Data sources and licences

The MIT licence covers the code in this repository only.

| Source | Used for | Licence |
| --- | --- | --- |
| [Ensembl / BioMart](https://www.ensembl.org/info/about/legal/disclaimer.html) | Gene, transcript, exon and GO annotation for human chr20–21 | No restrictions on use; EMBL-EBI terms |
| [Gene Ontology](http://geneontology.org/docs/go-citation-policy/) | GO term names and namespaces, via BioMart | CC BY 4.0 |
| [bedtools](https://github.com/arq5x/bedtools2) | Independent validation of the overlap results (optional) | MIT |

**What this repository ships.** `data/` contains three unmodified BioMart exports, gzipped, so the database builds offline and reproducibly. Everything in `results/` is generated from them by the code here. If you reuse these results, cite Ensembl and the Gene Ontology as well.

## References

1. Harrison, P. W. *et al.* (2024). Ensembl 2024. *Nucleic Acids Research* **52**, D891–D899.
2. The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. *Genetics* **224**, iyad031.
3. Codd, E. F. (1970). A relational model of data for large shared data banks. *Communications of the ACM* **13**, 377–387.
4. Kent, W. J. *et al.* (2002). The human genome browser at UCSC. *Genome Research* **12**, 996–1006. (the binning scheme)
5. Guttman, A. (1984). R-trees: a dynamic index structure for spatial searching. *SIGMOD* **14**, 47–57.
6. Beckmann, N. *et al.* (1990). The R*-tree: an efficient and robust access method for points and rectangles. *SIGMOD* **19**, 322–331.
7. Alekseyenko, A. V. & Lee, C. J. (2007). Nested Containment List (NCList). *Bioinformatics* **23**, 1386–1393.
8. Hipp, D. R. *et al.* SQLite. https://www.sqlite.org/

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
