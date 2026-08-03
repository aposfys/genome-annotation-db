| Query | Rows | Without indexes | With indexes | Speed-up | VM steps (without → with) | Work ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Q1 — Transcripts of a named gene | 2 | 1.2 ms | 0.0 ms | 45.9× | 59,746 → 61 | 979.4× |
| Q2 — Most heavily transcribed genes | 20 | 6.0 ms | 3.4 ms | 1.7× | 281,488 → 230,199 | 1.2× |
| Q3 — Transcripts with the most exons | 30 | 48.9 ms | 47.5 ms | 1.0× | 2,810,638 → 2,810,638 | 1.0× |
| Q4 — Genes carrying at least 20 GO terms | 219 | 3.9 ms | 3.8 ms | 1.0× | 290,146 → 290,146 | 1.0× |
| Q5 — Annotation depth per GO namespace | 3 | 4.6 ms | 4.5 ms | 1.0× | 258,238 → 176,553 | 1.5× |
| Q6 — Structure of genes carrying a GO term | 13 | 4.9 ms | 1.0 ms | 5.0× | 176,706 → 79,462 | 2.2× |
| Q7 — Exon reuse across transcripts, per gene | 702 | 66.8 ms | 59.6 ms | 1.1× | 3,394,126 → 3,394,126 | 1.0× |
