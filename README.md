# Genome Annotation Database: Schema Design and Index Benchmarking
A normalised relational database over Ensembl gene annotation, with a validating loader and a benchmark measuring what the indexes are actually worth.

[![Pipeline](https://github.com/aposfys/genome-annotation-db/actions/workflows/pipeline.yml/badge.svg)](https://github.com/aposfys/genome-annotation-db/actions/workflows/pipeline.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Genes, their transcripts, the exons those transcripts are assembled from, and Gene Ontology terms for human chromosomes 20 and 21. That is 769 protein-coding genes, 9,950 transcripts, 26,970 exons, 97,742 transcript–exon links and 12,596 GO assignments, exported from Ensembl BioMart. Clone it and you have a working 151,844-row SQLite database in under a second, with no server and no credentials. A MySQL 8 connection path exists for building and querying, but CI does not run it.

### Three results

**Indexes are not a blanket win.** Across seven analysis queries, the indexes cut the work of a point lookup from 59,746 SQLite VM instructions to 61 (979× less work, 44× in wall-clock time in the committed run). Q3, Q4 and Q7 execute the same number of VM instructions with and without indexes, and their plans do not change, because they aggregate over every row and there is nothing to skip. Q5 does 1.5× less work for a 1.1× wall-clock gain. An index earns its keep in proportion to how much of the table it lets a query skip.

**A structure that wins in isolation lost in place, and the query plan shows why.** On a bare interval table of every human gene (78,733, laid end to end on one axis with no chromosome filter, which is the B-tree's worst case), the R\*Tree answered 73× faster than a B-tree scan (22.0 vs 1,615.1 µs). Against the real database it came last, at 363.6 µs per query against 21.8 µs for the B-tree. The committed plan explains it. SQLite starts from the chromosome index on the R\*Tree's mapping table and then looks up each gene in the R\*Tree by its row id, so the spatial search is never used. The loss measures that join order, not the R\*Tree. I have not yet measured the same query with the R\*Tree forced to drive the join.

**Validating against bedtools found a real off-by-one.** Ensembl coordinates are 1-based inclusive, but the overlap predicate had been written with the half-open `<`/`>` rule, so genes ending exactly where a window starts were silently excluded. All three interval strategies agreed with each other and all three were wrong. The B-tree strategy now agrees with bedtools on 600 of 600 windows, 400 of them placed deliberately on gene boundaries. Random windows almost never land on an edge and agreed even before the fix. The same lesson caught a second bug later. The R\*Tree stored its boxes as 32-bit floats, which round coordinates above 16.7 Mb, so it disagreed with the B-tree on boundary windows while agreeing on random ones. It now stores integers, and a test covers it.

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

The core package has **no third-party dependencies**. SQLite is in the standard library.

```
genomedb gene DEFB125              # transcripts of a gene
genomedb go GO:0007186             # genes carrying a GO term
genomedb region 20 1000000 2000000 --strategy rtree
```

### Prior work

Genomic interval indexing is a solved and well-published problem, and none of the three
strategies compared here is original to this repository.

- **UCSC binning**, the hierarchical bin scheme behind the UCSC Genome Browser, BEDTools and
  SAMtools, which keeps the bin alongside the coordinate so an overlap query pre-filters on
  bins.
- **Binary Interval Search (BITS)**, *Bioinformatics* 2013, and **augmented range trees**,
  *Scientific Reports* 2019, the scalable alternatives, both benchmarked at genome scale.
- **Segment trees** (Segtor, *PLOS One* 2011) for annotating coordinates and variants.

The R\*Tree result above is an observation about one SQL query and the plan SQLite chose for
it, not a finding about interval structures. At 769 genes this is also far below the scale
those papers benchmark at.

Read it as a schema-design and benchmarking exercise that measures a known trade-off and
validates its interval logic against bedtools, catching a real off-by-one in the process.
That validation is the part worth keeping.

### More

- [Benchmarks in full: index timings, interval structures, empirical complexity](docs/BENCHMARKS.md)
- [Schema, normalisation, data quality and design decisions](docs/SCHEMA.md)
- [Data sources, provenance, licences and references](docs/DATA.md)

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
