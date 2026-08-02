-- Q3. Transcripts with the most exons.
-- n_exons and max_exon_rank should agree; where they do not, the exon ranks in
-- the source are not contiguous, so reporting both is a data-quality check as
-- much as a result.
SELECT
    t.transcript_id,
    t.gene_id,
    g.gene_name,
    COUNT(te.exon_id)  AS n_exons,
    MAX(te.exon_rank)  AS max_exon_rank
FROM transcript AS t
JOIN gene AS g ON g.gene_id = t.gene_id
LEFT JOIN transcript_exon AS te ON te.transcript_id = t.transcript_id
GROUP BY t.transcript_id, t.gene_id, g.gene_name
ORDER BY n_exons DESC, t.transcript_id
LIMIT 30;
