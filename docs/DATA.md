# Data sources, licences and references

The MIT licence covers the code in this repository only.

| Source | Used for | Licence |
| --- | --- | --- |
| [Ensembl / BioMart](https://www.ensembl.org/info/about/legal/disclaimer.html) | Gene, transcript, exon and GO annotation for human chr20–21 | No restrictions on use; EMBL-EBI terms |
| [Gene Ontology](http://geneontology.org/docs/go-citation-policy/) | GO term names and namespaces, via BioMart | CC BY 4.0 |
| [bedtools](https://github.com/arq5x/bedtools2) | Independent validation of the overlap results (optional) | MIT |

**What this repository ships.** `data/` contains four unmodified BioMart exports, gzipped,
3.1 MB in all, so the database builds offline and reproducibly. Three of them (`genes`,
`transcripts_exons` and `go`) cover the 769 protein-coding genes on chromosomes 20 and 21 and
feed the database. The fourth, `genes_all`, lists all 78,733 genes of every biotype on the
assembled human chromosomes and feeds only the scaling study. Everything in `results/` is
generated from these files by the code here. If you reuse these results, cite Ensembl and the Gene Ontology as well.

## Provenance

The exports were downloaded by hand through the BioMart web interface, and the Ensembl
release was not recorded at the time. What the files themselves show is limited.

- `genes` and `genes_all` share all 769 chromosome 20 and 21 gene IDs, but 387 of those genes
  have different start or end coordinates in the two files. They therefore come from
  different Ensembl releases.
- `genomedb fetch` holds the queries in [`biomart.py`](../src/genomedb/biomart.py), but it does
  not reproduce these files. Its pinned archive URL does not resolve to a BioMart service, its
  chromosome 20 and 21 gene query has no protein-coding filter, and it asks for unique rows,
  while the committed GO export repeats each gene–GO pair about five times.

Every number in this repository comes from the committed files as they are, so the results
reproduce from a clone. Regenerating the exports from one recorded release would close the
provenance gap, and would change the row counts.

## References

1. Harrison, P. W. *et al.* (2024). Ensembl 2024. *Nucleic Acids Research* **52**, D891–D899.
2. The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. *Genetics*
   **224**, iyad031.
3. Codd, E. F. (1970). A relational model of data for large shared data banks.
   *Communications of the ACM* **13**, 377–387.
4. Kent, W. J. *et al.* (2002). The human genome browser at UCSC. *Genome Research* **12**,
   996–1006. (the binning scheme)
5. Guttman, A. (1984). R-trees: a dynamic index structure for spatial searching. *SIGMOD*
   **14**, 47–57.
6. Beckmann, N. *et al.* (1990). The R*-tree: an efficient and robust access method for
   points and rectangles. *SIGMOD* **19**, 322–331.
7. Alekseyenko, A. V. & Lee, C. J. (2007). Nested Containment List (NCList).
   *Bioinformatics* **23**, 1386–1393.
8. Hipp, D. R. *et al.* SQLite. https://www.sqlite.org/
