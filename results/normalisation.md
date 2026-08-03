| Relation | Candidate key | Dependencies | BCNF |
| --- | --- | ---: | --- |
| `gene` | (gene_id) | 1 | yes |
| `transcript` | (transcript_id) | 1 | yes |
| `exon` | (exon_id) | 1 | yes |
| `go_term` | (go_id) | 1 | yes |
| `transcript_exon` | (transcript_id, exon_id) | 1 | yes |
| `gene_go` | (gene_id, go_id) | 0 | yes |

Normalised aggregate: 17.025 ms · materialised column: 3.746 ms · **4.5× faster**
