-- Q4. The most densely annotated genes: those carrying at least 20 GO terms.
SELECT
    g.gene_id,
    g.gene_name,
    COUNT(DISTINCT gg.go_id) AS n_go_terms
FROM gene AS g
JOIN gene_go AS gg ON gg.gene_id = g.gene_id
GROUP BY g.gene_id, g.gene_name
HAVING COUNT(DISTINCT gg.go_id) >= 20
ORDER BY n_go_terms DESC, g.gene_id;
