-- Genome annotation schema: genes, their transcripts, the exons those
-- transcripts are built from, and Gene Ontology annotation.
--
-- Written to run unchanged on SQLite and MySQL 8. Two choices make that work:
-- CHECK constraints instead of MySQL's ENUM, and VARCHAR rather than any
-- engine-specific string type.
--
-- Indexes live in indexes.sql, deliberately separate, so their effect can be
-- measured rather than assumed. See `genomedb benchmark`.

CREATE TABLE gene (
    gene_id     VARCHAR(30)  PRIMARY KEY,
    gene_name   VARCHAR(255),                      -- nullable: not every Ensembl gene is named
    chrom       VARCHAR(10)  NOT NULL,
    gene_start  INTEGER      NOT NULL,
    gene_end    INTEGER      NOT NULL,
    strand      INTEGER      NOT NULL,
    biotype     VARCHAR(50)  NOT NULL,

    CONSTRAINT ck_gene_strand   CHECK (strand IN (-1, 1)),
    CONSTRAINT ck_gene_span     CHECK (gene_end >= gene_start)
);

CREATE TABLE transcript (
    transcript_id VARCHAR(30) PRIMARY KEY,
    gene_id       VARCHAR(30) NOT NULL,
    tx_start      INTEGER     NOT NULL,
    tx_end        INTEGER     NOT NULL,
    tx_length     INTEGER,                         -- including UTRs and CDS

    CONSTRAINT fk_transcript_gene
        FOREIGN KEY (gene_id) REFERENCES gene(gene_id) ON DELETE CASCADE,
    CONSTRAINT ck_transcript_span CHECK (tx_end >= tx_start)
);

CREATE TABLE exon (
    exon_id    VARCHAR(30) PRIMARY KEY,
    exon_start INTEGER     NOT NULL,
    exon_end   INTEGER     NOT NULL,

    CONSTRAINT ck_exon_span CHECK (exon_end >= exon_start)
);

CREATE TABLE go_term (
    go_id        VARCHAR(20)  PRIMARY KEY,
    go_name      VARCHAR(255) NOT NULL,
    go_namespace VARCHAR(2)   NOT NULL,

    -- MySQL would express this as ENUM('BP','MF','CC'); a CHECK is equivalent
    -- here and runs on both engines.
    CONSTRAINT ck_go_namespace CHECK (go_namespace IN ('BP', 'MF', 'CC'))
);

-- An exon can belong to several transcripts of the same gene, which is what
-- alternative splicing is. Modelling that as a junction table rather than
-- duplicating exon rows per transcript is what makes Q7 -- how much exon
-- sharing a gene shows -- answerable at all.
CREATE TABLE transcript_exon (
    transcript_id VARCHAR(30) NOT NULL,
    exon_id       VARCHAR(30) NOT NULL,
    exon_rank     INTEGER     NOT NULL,

    PRIMARY KEY (transcript_id, exon_id),
    CONSTRAINT fk_tx_exon_transcript
        FOREIGN KEY (transcript_id) REFERENCES transcript(transcript_id) ON DELETE CASCADE,
    CONSTRAINT fk_tx_exon_exon
        FOREIGN KEY (exon_id) REFERENCES exon(exon_id) ON DELETE CASCADE,
    CONSTRAINT ck_exon_rank CHECK (exon_rank >= 1)
);

CREATE TABLE gene_go (
    gene_id VARCHAR(30) NOT NULL,
    go_id   VARCHAR(20) NOT NULL,

    PRIMARY KEY (gene_id, go_id),
    CONSTRAINT fk_gene_go_gene
        FOREIGN KEY (gene_id) REFERENCES gene(gene_id) ON DELETE CASCADE,
    CONSTRAINT fk_gene_go_term
        FOREIGN KEY (go_id) REFERENCES go_term(go_id) ON DELETE CASCADE
);
