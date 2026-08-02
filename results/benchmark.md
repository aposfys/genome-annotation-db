| Query | Rows | Without indexes | With indexes | Speed-up | Index used |
| --- | ---: | ---: | ---: | ---: | --- |
| Q1 — Transcripts of a named gene | 2 | 1.1 ms | 0.0 ms | 47.1× | yes |
| Q2 — Most heavily transcribed genes | 20 | 5.6 ms | 3.3 ms | 1.7× | yes |
| Q3 — Transcripts with the most exons | 30 | 46.4 ms | 47.3 ms | 1.0× | no |
| Q4 — Genes carrying at least 20 GO terms | 219 | 3.8 ms | 4.1 ms | 0.9× | no |
| Q5 — Annotation depth per GO namespace | 3 | 4.4 ms | 4.3 ms | 1.0× | yes |
| Q6 — Structure of genes carrying a GO term | 13 | 4.8 ms | 1.0 ms | 4.9× | yes |
| Q7 — Exon reuse across transcripts, per gene | 702 | 65.7 ms | 63.4 ms | 1.0× | no |
