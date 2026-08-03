| Strategy | Per query | Total | Rows | Structure |
| --- | ---: | ---: | ---: | --- |
| btree | 0.011 ms | 10 ms | 7,386 | Composite B-tree on (chrom, gene_start), range scan |
| binning | 0.016 ms | 14 ms | 7,386 | UCSC hierarchical binning scheme, pure SQL |
| rtree | 0.192 ms | 173 ms | 7,386 | SQLite R*Tree multidimensional index |
