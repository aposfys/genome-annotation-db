-- Q6. Transcript and exon structure of the genes carrying a given GO term.
-- Defaults to GO:0007186, G protein-coupled receptor signalling pathway.
WITH gene_transcript_exons AS (
    SELECT
        gg.gene_id,
        t.transcript_id,
        COUNT(te.exon_id) AS n_exons
    FROM gene_go AS gg
    JOIN transcript AS t ON t.gene_id = gg.gene_id
    LEFT JOIN transcript_exon AS te ON te.transcript_id = t.transcript_id
    WHERE gg.go_id = ?
    GROUP BY gg.gene_id, t.transcript_id
)
SELECT
    g.gene_id,
    g.gene_name,
    COUNT(DISTINCT gte.transcript_id)  AS n_transcripts,
    ROUND(AVG(CAST(gte.n_exons AS REAL)), 2) AS avg_exons_per_transcript
FROM gene AS g
JOIN gene_transcript_exons AS gte ON gte.gene_id = g.gene_id
GROUP BY g.gene_id, g.gene_name
ORDER BY n_transcripts DESC, g.gene_id;
