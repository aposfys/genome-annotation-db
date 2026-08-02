-- Q1. Every transcript of a named gene, in positional order.
-- The gene symbol is a parameter rather than a literal, so the query is
-- reusable and cannot be broken by a quoting mistake.
SELECT
    g.gene_id,
    g.gene_name,
    t.transcript_id,
    t.tx_start,
    t.tx_end,
    t.tx_length
FROM gene AS g
JOIN transcript AS t ON t.gene_id = g.gene_id
WHERE g.gene_name = ?
ORDER BY t.tx_start;
