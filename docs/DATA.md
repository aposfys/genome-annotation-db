# Data sources, licences and references

The MIT licence covers the code in this repository only.

| Source | Used for | Licence |
| --- | --- | --- |
| [Ensembl / BioMart](https://www.ensembl.org/info/about/legal/disclaimer.html) | Gene, transcript, exon and GO annotation for human chr20–21 | No restrictions on use; EMBL-EBI terms |
| [Gene Ontology](http://geneontology.org/docs/go-citation-policy/) | GO term names and namespaces, via BioMart | CC BY 4.0 |
| [bedtools](https://github.com/arq5x/bedtools2) | Independent validation of the overlap results (optional) | MIT |

**What this repository ships.** `data/` contains three unmodified BioMart exports, gzipped,
so the database builds offline and reproducibly. Everything in `results/` is generated from
them by the code here. If you reuse these results, cite Ensembl and the Gene Ontology as
well.

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
