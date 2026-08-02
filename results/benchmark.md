| Query | Rows | Without indexes | With indexes | Speed-up | Index used |
| --- | ---: | ---: | ---: | ---: | --- |
| Q1 — Transcripts of a named gene | 2 | 1.1 ms | 0.0 ms | 46.1× | yes |
| Q2 — Most heavily transcribed genes | 20 | 5.8 ms | 3.3 ms | 1.7× | yes |
| Q3 — Transcripts with the most exons | 30 | 48.5 ms | 45.8 ms | 1.1× | no |
| Q4 — Genes carrying at least 20 GO terms | 219 | 3.7 ms | 3.6 ms | 1.0× | no |
| Q5 — Annotation depth per GO namespace | 3 | 4.4 ms | 4.2 ms | 1.1× | yes |
| Q6 — Structure of genes carrying a GO term | 13 | 4.7 ms | 1.0 ms | 4.7× | yes |
| Q7 — Exon reuse across transcripts, per gene | 702 | 65.3 ms | 58.0 ms | 1.1× | no |
