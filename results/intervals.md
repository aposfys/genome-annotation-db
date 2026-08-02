| Strategy | Per query | Total | Rows | Structure |
| --- | ---: | ---: | ---: | --- |
| btree | 0.011 ms | 10 ms | 7,386 | Composite B-tree on (chrom, gene_start), range scan |
| binning | 0.015 ms | 14 ms | 7,386 | UCSC hierarchical binning scheme, pure SQL |
| rtree | 0.201 ms | 181 ms | 7,386 | SQLite R*Tree multidimensional index |
