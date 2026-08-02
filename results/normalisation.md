| Relation | Candidate key | Dependencies | BCNF |
| --- | --- | ---: | --- |
| `gene` | (gene_id) | 1 | yes |
| `transcript` | (transcript_id) | 1 | yes |
| `exon` | (exon_id) | 1 | yes |
| `go_term` | (go_id) | 1 | yes |
| `transcript_exon` | (transcript_id, exon_id) | 1 | yes |
| `gene_go` | (gene_id, go_id) | 0 | yes |

Normalised aggregate: 8.507 ms · materialised column: 1.943 ms · **4.4× faster**
