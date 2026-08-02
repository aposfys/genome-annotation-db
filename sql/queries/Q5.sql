-- Q5. Annotation depth per GO namespace.
-- Biological process dominates in every GO release; the comparison is a check
-- that the load preserved that expected shape.
SELECT
    gt.go_namespace,
    COUNT(DISTINCT gg.gene_id)                                       AS n_genes,
    COUNT(gg.go_id)                                                  AS total_annotations,
    ROUND(CAST(COUNT(gg.go_id) AS REAL) / COUNT(DISTINCT gg.gene_id), 2) AS avg_go_per_gene
FROM go_term AS gt
JOIN gene_go AS gg ON gg.go_id = gt.go_id
GROUP BY gt.go_namespace
ORDER BY n_genes DESC;
