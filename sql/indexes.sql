-- Secondary indexes, kept separate from the table definitions so that
-- `genomedb benchmark` can drop and recreate them and measure what they buy.
--
-- The two engines do not start from the same place. MySQL/InnoDB creates an
-- index for every foreign key automatically, so a schema with no explicit
-- indexes is still indexed on its join columns there. SQLite creates nothing
-- beyond primary keys, so the identical schema does full table scans for the
-- same joins. A schema that performs acceptably on one engine can therefore be
-- unusable on the other, which is the case this project started from.

-- Point lookups by gene symbol: Q1 and the `gene` CLI command filter on it.
CREATE INDEX idx_gene_name ON gene(gene_name);

-- Q2 filters genes by chromosome, then groups.
CREATE INDEX idx_gene_chrom ON gene(chrom);

-- The gene -> transcript join, used by nearly every query.
CREATE INDEX idx_transcript_gene ON transcript(gene_id);

-- transcript_exon's primary key already covers lookups by transcript_id, but
-- not the reverse. Q7 walks exons back to the transcripts that use them.
CREATE INDEX idx_transcript_exon_exon ON transcript_exon(exon_id);

-- gene_go's primary key covers gene -> GO. Q6 and the `go` CLI command go the
-- other way, from a GO term to the genes carrying it.
CREATE INDEX idx_gene_go_term ON gene_go(go_id);

-- Q5 aggregates annotations per namespace.
CREATE INDEX idx_go_term_namespace ON go_term(go_namespace);
