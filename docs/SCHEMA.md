# Schema, normalisation and data quality

## The schema

```
gene ──< transcript ──< transcript_exon >── exon
  │
  └──< gene_go >── go_term
```

Six tables, third normal form, with two junction tables carrying the many-to-many
relationships.

**The exon model is the design decision that matters.** An exon is shared between
transcripts rather than duplicated per transcript — 26,970 distinct exons support 97,742
transcript–exon links, so the average exon is used by 3.6 transcripts. Modelling that as a
junction table is what makes Q7 answerable at all: *how much of a gene's exon set is
constitutive, and how much is alternative?* PLCB4 turns out to have 172 distinct exons of
which 98 are shared across two or more transcripts.

Duplicating exons per transcript would have been simpler to load and would have made that
question unanswerable without string-matching coordinates.

## Normalisation: proved, then priced

`make normalisation` states the functional dependencies, checks Boyce-Codd Normal Form
against them mechanically, and confirms the dependencies actually hold in the loaded data.

| Relation | Candidate key | BCNF |
| --- | --- | --- |
| `gene` | (gene_id) | yes |
| `transcript` | (transcript_id) | yes |
| `exon` | (exon_id) | yes |
| `go_term` | (go_id) | yes |
| `transcript_exon` | (transcript_id, exon_id) | yes |
| `gene_go` | (gene_id, go_id) | yes |

BCNF holds when every non-trivial dependency has a superkey on its left-hand side. Here each
table's dependencies are keyed on its primary key, and `gene_go` is an all-key relation with
no non-key attribute to depend on anything.

The interesting case is `transcript_exon`. `exon_rank` depends on the **whole** composite
key, not on `exon_id` alone: the same exon can be third in one transcript and first in
another. Had rank been stored on `exon`, the schema would have been wrong in a way that only
shows up on alternatively spliced genes.

### What that property costs

Counting a transcript's exons means joining the junction table every time:

| | 9,950 transcripts |
| --- | ---: |
| Join every time (normalised) | 8.51 ms |
| Materialised column | 1.94 ms |
| | **4.4× faster** |

So normalisation costs about 6.6 ms on this query. What it buys is that the answer cannot be
wrong: nothing in the schema can keep a materialised count true, and any write to
`transcript_exon` that forgets to update it leaves the two disagreeing — a constraint cannot
express that dependency. The 4.4× is the price of that guarantee, stated rather than
assumed.

## Data quality

The loader validates as it goes and reports what it refused, because silently dropping
malformed rows leaves a database that looks complete and is not:

```
Inserted:
  exon                 26,970      gene                    769
  gene_go              12,596      go_term               3,817
  transcript            9,950      transcript_exon      97,742
Rejected:
  go: unmapped or missing domain      3,483   e.g. GO:0042421 domain='empty'
  go: no accession                    3,352   e.g. gene ENSG00000260861
```

**6,835 of the 69,501 GO rows — just under 10% — are unusable**, carrying an accession with
no domain or no accession at all. That is a property of the BioMart export, not a bug, but a
loader that discards them without saying so is hiding the fact that a tenth of the
annotation never arrived.

Nine integrity checks then run against the loaded database: orphan rows in each junction
table, namespaces outside BP/MF/CC, invalid strands, inverted coordinate spans, transcripts
extending beyond their gene, transcripts with no exons, and non-contiguous exon ranks. All
nine pass.

The foreign-key check is not redundant, incidentally: **SQLite accepts
`PRAGMA foreign_keys = OFF` silently**, and a database loaded that way looks fine until
something reads it. A test asserts the constraints are actually enforced.

## The discrepancy, explained

An earlier result file for the GO namespace query contains a single row — biological process
only, 769 genes, 13,939 annotations. The rebuilt database reports all three namespaces:

| Namespace | Genes | Annotations | Mean per gene |
| --- | ---: | ---: | ---: |
| Biological process | 639 | 5,638 | 8.82 |
| Cellular component | 728 | 3,602 | 4.95 |
| Molecular function | 701 | 3,356 | 4.79 |

Rather than leave that unexplained, the old numbers were reproduced from the shipped export.
**Ignoring the GO domain entirely — collapsing every term into one namespace — gives 757
genes and 13,437 annotations**, against the original's 769 and 13,939. That is within 1.6%
and 3.6%, and it reproduces the *shape* exactly.

So there were two causes, not one:

- **The namespace was collapsed.** Respecting the GO domain gives 12,596 pairs across three
  namespaces; ignoring it gives 13,437 in one. The original is the second shape.
- **The source export was slightly larger.** The residual — 12 genes and 502 annotations —
  is consistent with the original having been loaded from a marginally different download.

The second half is why `genomedb fetch` exists. The exports were originally produced by hand
through the BioMart web interface, so there was no way to tell whether a disagreement came
from the code or the data. The queries are now in
[`biomart.py`](../src/genomedb/biomart.py), pinned to **Ensembl release 113** through the
archive URL, and any export can be regenerated byte-for-byte.

## Design decisions

- **Credentials come from the environment, never the source.** `$GENOMEDB_URL` or
  `--db-url`. A host, user and password committed to a file are in version control
  permanently, and are exactly what an ETL script tends to accumulate.
- **Queries are parameterised.** Q1 and Q6 take their gene symbol and GO term as bound
  parameters rather than SQL literals, so they are reusable and cannot be broken by a
  quoting mistake.
- **`CHECK` constraints instead of MySQL's `ENUM`.** Equivalent expressiveness, and it runs
  on both engines. Constraints also encode the domain rules — strand is ±1, spans are
  non-inverted, exon ranks start at 1 — so bad data is rejected by the database rather than
  by convention.
- **No columns the loader cannot fill.** `transcript_name` and `biotype` were declared but
  never populated, showing up as blank columns in every result. Dropped. `tx_length` was
  present in the source and unused; added.
- **Every result file has a header row.** Two of the seven previously did not, which makes
  them awkward to read back programmatically.
- **The benchmark drops indexes before timing them, in that order.** Timing the indexed pass
  second means it cannot be flattered by a cache the unindexed pass warmed.

## Running against MySQL

```
export GENOMEDB_URL='mysql://user:password@host/database'
pip install -e ".[mysql]"
make build
```

`mysql-connector-python` is an optional extra, imported only if a MySQL URL is used.

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
