-- Q2. The most heavily transcribed genes on a chromosome.
-- An inner join would be equivalent given the HAVING clause, but LEFT JOIN
-- states the intent: count transcripts per gene, then keep those that have any.
SELECT
    g.chrom,
    g.gene_id,
    g.gene_name,
    COUNT(t.transcript_id) AS n_transcripts
FROM gene AS g
LEFT JOIN transcript AS t ON t.gene_id = g.gene_id
GROUP BY g.chrom, g.gene_id, g.gene_name
HAVING COUNT(t.transcript_id) > 0
ORDER BY n_transcripts DESC, g.gene_id
LIMIT 20;
