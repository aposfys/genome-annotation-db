| Query | Rows | Without indexes | With indexes | Speed-up | Index used |
| --- | ---: | ---: | ---: | ---: | --- |
| Q1 — Transcripts of a named gene | 2 | 1.2 ms | 0.0 ms | 46.5× | yes |
| Q2 — Most heavily transcribed genes | 20 | 6.0 ms | 3.5 ms | 1.7× | yes |
| Q3 — Transcripts with the most exons | 30 | 49.6 ms | 47.9 ms | 1.0× | no |
| Q4 — Genes carrying at least 20 GO terms | 219 | 3.8 ms | 3.8 ms | 1.0× | no |
| Q5 — Annotation depth per GO namespace | 3 | 4.5 ms | 4.3 ms | 1.0× | yes |
| Q6 — Structure of genes carrying a GO term | 13 | 4.8 ms | 1.0 ms | 5.0× | yes |
| Q7 — Exon reuse across transcripts, per gene | 702 | 67.8 ms | 59.9 ms | 1.1× | no |
