# Genome Annotation Database: Schema Design and Index Benchmarking
A normalised relational database over Ensembl gene annotation, with a validating loader and a benchmark measuring what the indexes are actually worth.

[![Pipeline](https://github.com/aposfys/genome-annotation-db/actions/workflows/pipeline.yml/badge.svg)](https://github.com/aposfys/genome-annotation-db/actions/workflows/pipeline.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Genes, their transcripts, the exons those transcripts are assembled from, and Gene Ontology terms — 769 genes, 9,950 transcripts, 26,970 exons, 97,742 transcript–exon links, 12,596 GO assignments, from Ensembl BioMart release 113. Clone it and you have a working 148,000-row database in under a second. No server, no credentials. SQLite by default; the same schema and queries run on MySQL 8.

### Three results

**Indexes are not a blanket win.** Across seven analysis queries, indexes buy 48× on a point lookup and *nothing* on four of them. Q3, Q4 and Q7 execute a byte-identical number of SQLite VM instructions with and without indexes — not "no measurable difference" but proof the plan did not change, because those queries aggregate over every row and there is nothing to skip. An index earns its keep in proportion to how much of the table it lets you skip.

**A structure that wins in isolation can lose in place.** On interval overlap, an R\*Tree beats a B-tree 80× at whole-genome scale (78,733 genes). Against the real database it comes *last*, 17× slower than the naive B-tree: the bare R\*Tree search is the fastest of the three at 12.8 µs, but the two joins back to `gene` needed to recover the identifier cost 382.5 µs. The index was never the bottleneck; the normalisation around it was. That is very likely why UCSC binning, which keeps coordinates and key in one table, is what the genome browsers actually use.

**Validating against bedtools found a real off-by-one.** Ensembl coordinates are 1-based inclusive; the overlap predicate had been written with the half-open `<`/`>` rule, so genes ending exactly where a window starts were silently excluded. All three interval strategies agreed with each other and all three were wrong. Agreement is now 600/600 windows, 400 of them placed deliberately on gene boundaries — random windows essentially never land on an edge and agreed even before the fix.

### Quick start

```
pip install -e ".[dev]"

make build       # create the schema and load the data (~0.6 s)
make queries     # run Q1-Q7, write results/
make benchmark   # time the queries with and without indexes
make intervals   # B-tree vs UCSC binning vs R*Tree on overlap queries
make validate    # cross-check overlap results against bedtools
make test
```

The core package has **no third-party dependencies** — SQLite is in the standard library.

```
genomedb gene DEFB125              # transcripts of a gene
genomedb go GO:0007186             # genes carrying a GO term
genomedb region 20 1000000 2000000 --strategy rtree
genomedb fetch                     # regenerate the Ensembl exports, pinned release
```

### Prior work

Genomic interval indexing is a solved and well-published problem, and none of the three
strategies compared here is original to this repository:

- **UCSC binning** — the hierarchical bin scheme behind the UCSC Genome Browser, BEDTools and
  SAMtools, which keeps the bin alongside the coordinate so an overlap query pre-filters on
  bins. Its known weakness is exactly the one seen here: when query intervals share the same
  bias, few bins are examined and each contains many intervals to enumerate.
- **Binary Interval Search (BITS)**, *Bioinformatics* 2013, and **augmented range trees**,
  *Scientific Reports* 2019 — the scalable alternatives, both benchmarked at genome scale.
- **Segment trees** (Segtor, *PLOS One* 2011) for annotating coordinates and variants.

The R\*Tree-loses-in-place result is therefore an implementation observation rather than a new
algorithmic finding, and its explanation — that keeping coordinates and key in one table
avoids the joins that dominate — is the reason the published schemes are shaped the way they
are. At 769 genes this is also far below the scale those papers benchmark at.

Read it as a schema-design and benchmarking exercise that measures a known trade-off
carefully, and validates its interval logic against bedtools, catching a real off-by-one in
the process. That validation is the part worth keeping.

### More

- [Benchmarks in full: index timings, interval structures, empirical complexity](docs/BENCHMARKS.md)
- [Schema, normalisation, data quality and design decisions](docs/SCHEMA.md)
- [Data sources, licences and references](docs/DATA.md)

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
