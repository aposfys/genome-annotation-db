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
| **Headline** | Indexes buy **47×** on a point lookup and **nothing at all** on four of seven queries |

## The result: indexes are not a blanket win

`genomedb benchmark` drops every secondary index, times all seven queries, recreates the indexes and times them again. Same data, same queries, only the indexes change. Medians of seven runs:

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

<details>
<summary>Q1's query plan, before and after</summary>

```
without: SCAN t | SEARCH g USING INDEX sqlite_autoindex_gene_1 (gene_id=?) | USE TEMP B-TREE FOR ORDER BY
with:    SEARCH g USING INDEX idx_gene_name (gene_name=?) | SEARCH t USING INDEX idx_transcript_gene (gene_id=?)
```
</details>

### A portability trap worth knowing

**MySQL/InnoDB silently creates an index for every foreign key. SQLite creates nothing beyond primary keys.**

So a schema with no explicit indexes is already indexed on its join columns under MySQL, and does full scans for the identical joins under SQLite. A design that performs acceptably on one engine can be unusable on the other, and nothing in the DDL hints at it. Every index this schema relies on is therefore declared explicitly in [`sql/indexes.sql`](sql/indexes.sql) rather than left to the engine.

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
  benchmark.py  Index timing and query-plan capture
  cli.py        Subcommands: build, query, gene, go, check, benchmark
data/           Ensembl BioMart exports, gzipped (1.7 MB)
results/        Query output, load report, quality report, benchmark
tests/          pytest suite (24 tests)
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
4. Hipp, D. R. *et al.* SQLite. https://www.sqlite.org/

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
