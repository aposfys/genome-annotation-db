| Strategy | Per query | Total | Rows | Structure |
| --- | ---: | ---: | ---: | --- |
| btree | 0.023 ms | 21 ms | 7,386 | Composite B-tree on (chrom, gene_start), range scan |
| binning | 0.029 ms | 27 ms | 7,386 | UCSC hierarchical binning scheme, pure SQL |
| rtree | 0.377 ms | 339 ms | 7,386 | SQLite R*Tree multidimensional index |
