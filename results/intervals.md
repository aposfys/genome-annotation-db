| Strategy | Per query | Total | Rows | Structure |
| --- | ---: | ---: | ---: | --- |
| btree | 0.022 ms | 20 ms | 7,386 | Composite B-tree on (chrom, gene_start), range scan |
| binning | 0.031 ms | 27 ms | 7,386 | UCSC hierarchical binning scheme, pure SQL |
| rtree | 0.364 ms | 327 ms | 7,386 | SQLite R*Tree multidimensional index |
