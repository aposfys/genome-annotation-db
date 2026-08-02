-- Q7. How much each gene reuses its exons across transcripts.
-- An exon used by two or more transcripts of the same gene is a constitutive
-- exon; one used by a single transcript is where alternative splicing differs.
-- This is the query the many-to-many exon model exists to make possible.
WITH gene_exon_usage AS (
    SELECT
        t.gene_id,
        te.exon_id,
        COUNT(DISTINCT te.transcript_id) AS n_transcripts_using_exon
    FROM transcript AS t
    JOIN transcript_exon AS te ON te.transcript_id = t.transcript_id
    GROUP BY t.gene_id, te.exon_id
)
SELECT
    g.gene_id,
    g.gene_name,
    COUNT(DISTINCT geu.exon_id) AS total_unique_exons,
    SUM(CASE WHEN geu.n_transcripts_using_exon >= 2 THEN 1 ELSE 0 END) AS shared_exons
FROM gene AS g
JOIN gene_exon_usage AS geu ON geu.gene_id = g.gene_id
GROUP BY g.gene_id, g.gene_name
HAVING COUNT(DISTINCT geu.exon_id) > 1
ORDER BY shared_exons DESC, total_unique_exons DESC;
