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
| **Data** | Ensembl BioMart, human chromosomes 20 and 21 |
| **Scale** | 769 genes · 9,950 transcripts · 26,970 exons · 97,742 transcript–exon links · 12,596 GO assignments |
| **Engines** | SQLite by default; the same schema and queries run on MySQL 8 |
| **Headline** | Indexes buy **48×** on a point lookup and **nothing** on four of seven queries |
| **Also** | An R\*Tree beats a B-tree 14× on interval search — then loses 16× once the joins around it are counted |

## The result: indexes are not a blanket win

`genomedb benchmark` drops every secondary index, times all seven queries, recreates the indexes and times them again. Same data, same queries, only the indexes change. Medians of seven runs on an Apple M4:

| Query | Rows | Without indexes | With indexes | Speed-up | Index used |
| --- | ---: | ---: | ---: | ---: | --- |
| Q1 — Transcripts of a named gene | 2 | 1.1 ms | 0.02 ms | **47.1×** | yes |
| Q2 — Most heavily transcribed genes | 20 | 5.6 ms | 3.3 ms | 1.7× | yes |
| Q3 — Transcripts with the most exons | 30 | 46.4 ms | 47.3 ms | 1.0× | no |
| Q4 — Genes carrying ≥ 20 GO terms | 219 | 3.8 ms | 4.1 ms | 0.9× | no |
| Q5 — Annotation depth per GO namespace | 3 | 4.4 ms | 4.3 ms | 1.0× | yes |
| Q6 — Structure of genes carrying a GO term | 13 | 4.8 ms | 1.0 ms | **4.9×** | yes |
| Q7 — Exon reuse across transcripts | 702 | 65.7 ms | 63.4 ms | 1.0× | no |

**The split is entirely predictable once you look at the query plans**, which the benchmark captures alongside the timings.

Q1 and Q6 are *selective*: they find a handful of rows matching one value. Without an index the planner scans a whole table; with one it seeks straight to the match. Q1's plan goes from `SCAN t` to `SEARCH g USING INDEX idx_gene_name (gene_name=?)`, and 1.1 ms becomes 0.02 ms.

Q3, Q4 and Q7 aggregate over *every* row. There is nothing to seek to — the query has to visit the whole table either way — so the index adds maintenance cost and returns nothing. Q4 is marginally slower with indexes than without, which is the honest shape of that trade-off.

The lesson generalises: an index earns its keep in proportion to how much of the table it lets you skip. A report that touches everything is not a candidate.

The absolute figures are hardware-dependent — the same benchmark on a GitHub Actions runner gives Q1 a 21× speed-up rather than 48× — but the *shape* is not. Selective queries gain; full aggregations do not. Re-run `make benchmark` and you get your own numbers.

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

Measured on synthetic intervals drawn from a realistic log-normal length distribution, since the real chromosomes cannot supply the sizes:

| n intervals | B-tree scan | UCSC binning | R\*Tree | Fastest |
| ---: | ---: | ---: | ---: | --- |
| 1,000 | 26.9 µs | 18.7 µs | **10.3 µs** | R\*Tree |
| 4,000 | 81.0 µs | 30.0 µs | **15.5 µs** | R\*Tree |
| 16,000 | 352.7 µs | 74.7 µs | **32.7 µs** | R\*Tree |
| 64,000 | 1,359.8 µs | 288.4 µs | **106.3 µs** | R\*Tree |
| 256,000 | 6,171.0 µs | 2,406.8 µs | **443.3 µs** | R\*Tree |

The B-tree degrades linearly, as predicted — it is scanning. Binning is roughly 2.5× better. The R\*Tree is **14× faster than the B-tree at 256,000 intervals**, and its lead widens with n.

### The result that matters

Run the same comparison against the *real* database and the R\*Tree comes last, 16× slower than the naive B-tree. That contradiction is the most instructive thing in this project, so it was measured rather than explained away:

| On the real 769-gene database | µs per query |
| --- | ---: |
| R\*Tree, bare search | **13.3** |
| B-tree, single-table scan | 21.9 |
| R\*Tree **plus the two joins back to `gene`** | **369.3** |

**The R\*Tree search is the fastest of the three. The joins around it cost 27× more than the search does.**

An R\*Tree stores integer keys, so recovering the gene identifier and filtering by chromosome needs a mapping table and a join back to `gene` — and at this scale that dominates completely. The index was never the bottleneck; the normalisation around it was.

Two things follow, and neither is visible from a single benchmark:

- **Micro-benchmarks of a data structure can invert once it is embedded in a schema.** The structure that wins in isolation lost by 16× in place.
- **The crossover depends on the join cost, not just on n.** Binning wins in practice at this scale precisely because it keeps the coordinates and the key in one indexed table, so it needs no mapping hop. That is very likely why UCSC chose it.

### A portability trap worth knowing

**MySQL/InnoDB silently creates an index for every foreign key. SQLite creates nothing beyond primary keys.**

So a schema with no explicit indexes is already indexed on its join columns under MySQL, and does full scans for the identical joins under SQLite. A design that performs acceptably on one engine can be unusable on the other, and nothing in the DDL hints at it. Every index this schema relies on is therefore declared explicitly in [`sql/indexes.sql`](sql/indexes.sql) rather than left to the engine.

## Empirical complexity, not a single measurement

`genomedb scaling` times a point lookup across a geometric series of table sizes and classifies the growth by comparing what was *observed* against what each candidate law *predicts* over that range — O(n) predicts a 256-fold rise from 1k to 256k rows, O(log n) about 1.8-fold, O(1) none.

| n rows | Unindexed | Indexed | Speed-up |
| ---: | ---: | ---: | ---: |
| 1,000 | 37.4 µs | 7.1 µs | 5.3× |
| 4,000 | 128.5 µs | 7.1 µs | 18.1× |
| 16,000 | 496.1 µs | 7.5 µs | 66.2× |
| 64,000 | 1,986.1 µs | 7.3 µs | 271.9× |
| 256,000 | 9,649.0 µs | 8.2 µs | **1,180.7×** |

- **Unindexed: O(n).** Grew 258-fold across a 256-fold increase in rows — the textbook scan, with a linear fit at R² = 0.998.
- **Indexed: O(1) over this range.** Grew 1.15-fold. A B-tree seek is O(log n) in theory, but from 1k to 256k rows that predicts only a 1.8-fold rise, and the observed 1.15 is nearer to flat. The honest statement is that the measurement cannot separate O(log n) from O(1) here, not that logarithmic growth was disproved.

That second point is why the classifier compares against predicted growth rather than picking whichever least-squares fit scores a higher R². With five points a fit will always name a winner, including for a response that is really flat and noisy — and it did, calling a 1.2-fold rise "linear", before the method was changed.

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
make test
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

## A discrepancy worth recording

The rebuilt database reports GO annotation across all three namespaces:

| Namespace | Genes | Annotations | Mean per gene |
| --- | ---: | ---: | ---: |
| Biological process | 639 | 5,638 | 8.82 |
| Cellular component | 728 | 3,602 | 4.95 |
| Molecular function | 701 | 3,356 | 4.79 |

These figures reconcile exactly with an independent recount of the source export — 12,596 distinct gene–GO pairs, confirmed outside the loader entirely. The earlier result file for this query contained a single row, biological process only, with 13,939 annotations: more BP annotations alone than there are distinct pairs in the data shipped alongside it.

I could not reproduce that from these inputs and have not guessed at a cause; the most likely explanation is that those results were produced against a database loaded from a different snapshot. It is recorded here because a result set that cannot be regenerated from its stated inputs is worth flagging rather than quietly replacing.

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
  intervals.py  UCSC binning, R*Tree and B-tree overlap strategies
  scaling.py    Growth measurement and growth-law classification
  benchmark.py  Index timing and query-plan capture
  cli.py        Subcommands: build, query, gene, go, check, benchmark
data/           Ensembl BioMart exports, gzipped (1.7 MB)
results/        Query output, load report, quality report, benchmark
tests/          pytest suite (35 tests)
```

## Data sources and licences

The MIT licence covers the code in this repository only.

| Source | Used for | Licence |
| --- | --- | --- |
| [Ensembl / BioMart](https://www.ensembl.org/info/about/legal/disclaimer.html) | Gene, transcript, exon and GO annotation for human chr20–21 | No restrictions on use; EMBL-EBI terms |
| [Gene Ontology](http://geneontology.org/docs/go-citation-policy/) | GO term names and namespaces, via BioMart | CC BY 4.0 |

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
