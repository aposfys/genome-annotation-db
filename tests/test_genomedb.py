"""Tests for the schema, loader, queries and integrity checks."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from genomedb import benchmark, load, quality, queries
from genomedb.db import MYSQL, SQLITE, adapt, connect, split_statements

ROOT = Path(__file__).resolve().parents[1]
SQL_DIR = ROOT / "sql"
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"


# --- dialect handling --------------------------------------------------------


def test_placeholders_are_rewritten_for_mysql():
    sql = "SELECT * FROM gene WHERE gene_name = ? AND chrom = ?"
    assert adapt(sql, SQLITE) == sql
    assert adapt(sql, MYSQL) == ("SELECT * FROM gene WHERE gene_name = %s AND chrom = %s")


def test_statement_splitting_ignores_semicolons_in_strings():
    script = "SELECT ';' AS a; SELECT 2;"
    assert split_statements(script) == ["SELECT ';' AS a", "SELECT 2"]


def test_statement_splitting_ignores_semicolons_in_comments():
    script = "-- a comment; with a semicolon\nSELECT 1;\nSELECT 2;"
    statements = split_statements(script)
    assert len(statements) == 2
    assert statements[1] == "SELECT 2"


def test_statement_splitting_keeps_a_trailing_statement_without_semicolon():
    assert split_statements("SELECT 1; SELECT 2") == ["SELECT 1", "SELECT 2"]


# --- a small database built from fixtures ------------------------------------


def _write_tsv(
    path: Path, header: list[str], rows: list[list[str]], gzipped: bool = False
) -> Path:
    text = "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n"
    if gzipped:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def tiny_db(tmp_path):
    """A three-gene database, small enough to reason about exactly."""
    data = tmp_path / "data"
    data.mkdir()

    _write_tsv(
        data / "genes.tsv",
        [
            "Gene stable ID",
            "Gene name",
            "Chromosome/scaffold name",
            "Gene start (bp)",
            "Gene end (bp)",
            "Strand",
            "Gene type",
        ],
        [
            ["G1", "AAA", "20", "100", "900", "1", "protein_coding"],
            ["G2", "BBB", "21", "200", "800", "-1", "protein_coding"],
            # Strand 0 is invalid and must be rejected, not silently coerced.
            ["G3", "CCC", "21", "300", "700", "0", "lncRNA"],
        ],
    )
    _write_tsv(
        data / "transcripts_exons.tsv",
        [
            "Gene stable ID",
            "Transcript stable ID",
            "Transcript start (bp)",
            "Transcript end (bp)",
            "Transcript length (including UTRs and CDS)",
            "Exon stable ID",
            "Exon rank in transcript",
            "Exon region start (bp)",
            "Exon region end (bp)",
        ],
        [
            # T1 and T2 both use E1: the shared exon Q7 exists to find.
            ["G1", "T1", "100", "500", "400", "E1", "1", "100", "200"],
            ["G1", "T1", "100", "500", "400", "E2", "2", "400", "500"],
            ["G1", "T2", "100", "900", "800", "E1", "1", "100", "200"],
            ["G1", "T2", "100", "900", "800", "E3", "2", "800", "900"],
            ["G2", "T3", "200", "800", "600", "E4", "1", "200", "800"],
            # Belongs to the rejected gene G3, so must be rejected too.
            ["G3", "T4", "300", "700", "400", "E5", "1", "300", "700"],
        ],
    )
    _write_tsv(
        data / "go.tsv",
        [
            "Gene stable ID",
            "GO term accession",
            "GO term name",
            "GO term definition",
            "GO domain",
        ],
        [
            ["G1", "GO:0001", "alpha", '"some definition" [GOC:x]', "biological_process"],
            ["G1", "GO:0002", "beta", "another", "molecular_function"],
            ["G2", "GO:0001", "alpha", "def", "biological_process"],
            ["G2", "GO:0003", "gamma", "def", "cellular_component"],
            # Missing domain: the class of row the original silently dropped.
            ["G2", "GO:0004", "delta", "def", ""],
            # Missing accession entirely.
            ["G1", "", "", "", "biological_process"],
        ],
    )

    conn = connect(sqlite_path=tmp_path / "test.sqlite")
    report = load.build(conn, data, SQL_DIR)
    yield conn, report
    conn.close()


def test_invalid_strand_is_rejected(tiny_db):
    _, report = tiny_db
    assert report.inserted["gene"] == 2
    assert report.rejected["gene: strand not 1 or -1"] == 1


def test_rows_depending_on_a_rejected_gene_are_also_rejected(tiny_db):
    """A transcript of a rejected gene would otherwise violate the foreign key."""
    _, report = tiny_db
    assert report.rejected["transcript: gene_id not in gene table"] == 1
    assert report.inserted["transcript"] == 3


def test_incomplete_go_rows_are_counted_not_silently_dropped(tiny_db):
    _, report = tiny_db
    assert report.rejected["go: unmapped or missing domain"] == 1
    assert report.rejected["go: no accession"] == 1
    assert report.total_rejected == 4


def test_every_rejection_carries_an_example(tiny_db):
    _, report = tiny_db
    assert set(report.rejected) == set(report.examples)


def test_exons_are_shared_between_transcripts(tiny_db):
    """The many-to-many model is the point; E1 belongs to both T1 and T2."""
    conn, _ = tiny_db
    assert conn.scalar("SELECT COUNT(*) FROM exon") == 4
    assert conn.scalar("SELECT COUNT(*) FROM transcript_exon") == 5
    assert (
        conn.scalar(
            "SELECT COUNT(DISTINCT transcript_id) FROM transcript_exon WHERE exon_id = 'E1'"
        )
        == 2
    )


def test_transcript_length_is_populated(tiny_db):
    """The original schema declared columns the loader never filled."""
    conn, _ = tiny_db
    assert conn.scalar("SELECT COUNT(*) FROM transcript WHERE tx_length IS NULL") == 0


def test_integrity_checks_pass_on_a_clean_load(tiny_db):
    conn, _ = tiny_db
    report = quality.report(conn)
    assert report["all_passed"], report["failed"]


def test_foreign_keys_are_enforced(tiny_db):
    """SQLite disables foreign keys by default; the connection must turn them on."""
    import sqlite3

    conn, _ = tiny_db
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute(
            "INSERT INTO transcript (transcript_id, gene_id, tx_start, tx_end)"
            " VALUES ('T9', 'NO_SUCH_GENE', 1, 2)"
        )


def test_check_constraint_rejects_an_invalid_namespace(tiny_db):
    import sqlite3

    conn, _ = tiny_db
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute(
            "INSERT INTO go_term (go_id, go_name, go_namespace) VALUES ('GO:9', 'x', 'ZZ')"
        )


def test_check_constraint_rejects_an_inverted_span(tiny_db):
    import sqlite3

    conn, _ = tiny_db
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute(
            "INSERT INTO gene (gene_id, gene_name, chrom, gene_start, gene_end,"
            " strand, biotype) VALUES ('G9', 'x', '1', 500, 100, 1, 'protein_coding')"
        )


# --- queries -----------------------------------------------------------------


def test_every_query_has_a_file():
    for query in queries.QUERIES:
        assert (SQL_DIR / "queries" / query.filename).exists(), query.name


def test_parameterised_queries_accept_their_parameter(tiny_db):
    conn, _ = tiny_db
    _, rows = queries.run(conn, queries.BY_NAME["Q1"], SQL_DIR, ("AAA",))
    assert {row[2] for row in rows} == {"T1", "T2"}

    _, none = queries.run(conn, queries.BY_NAME["Q1"], SQL_DIR, ("NOT_A_GENE",))
    assert none == []


def test_q7_finds_the_shared_exon(tiny_db):
    conn, _ = tiny_db
    columns, rows = queries.run(conn, queries.BY_NAME["Q7"], SQL_DIR)
    by_gene = {row[0]: dict(zip(columns, row, strict=True)) for row in rows}
    assert by_gene["G1"]["total_unique_exons"] == 3
    assert by_gene["G1"]["shared_exons"] == 1


def test_all_queries_run_and_write_headers(tiny_db, tmp_path):
    conn, _ = tiny_db
    summary = queries.run_all(conn, SQL_DIR, tmp_path / "out")
    assert set(summary) == {q.name for q in queries.QUERIES}

    for query in queries.QUERIES:
        path = tmp_path / "out" / f"{query.name}.tsv"
        header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert header and all(header), f"{query.name} wrote an empty header"


def test_written_results_round_trip(tiny_db, tmp_path):
    conn, _ = tiny_db
    columns, rows = queries.run(conn, queries.BY_NAME["Q5"], SQL_DIR)
    path = queries.write_tsv(columns, rows, tmp_path / "Q5.tsv")
    with path.open(encoding="utf-8") as handle:
        read_back = list(csv.DictReader(handle, delimiter="\t"))
    assert len(read_back) == len(rows)
    assert set(read_back[0]) == set(columns)


# --- benchmark ---------------------------------------------------------------


def test_indexes_can_be_dropped_and_recreated(tiny_db):
    conn, _ = tiny_db
    names = benchmark.drop_indexes(conn, SQL_DIR)
    assert names, "indexes.sql declared no indexes"
    assert all(name.startswith("idx_") for name in names)
    benchmark.create_indexes(conn, SQL_DIR)


def test_benchmark_reports_a_speedup_for_every_query(tiny_db):
    conn, _ = tiny_db
    comparisons = benchmark.compare(conn, SQL_DIR, repeats=1)
    assert len(comparisons) == len(queries.QUERIES)
    assert all(c.speedup > 0 for c in comparisons)

    summary = benchmark.summarise(comparisons)
    assert summary["queries"] == len(queries.QUERIES)


def test_automatic_primary_key_indexes_do_not_count_as_index_use():
    """SQLite names its implicit indexes sqlite_autoindex_*; only idx_* is ours."""
    timing = benchmark.Timing(
        "Q", "t", 0.01, 1, "SEARCH g USING INDEX sqlite_autoindex_gene_1"
    )
    indexed = benchmark.Timing("Q", "t", 0.01, 1, "SEARCH g USING INDEX idx_gene_name")
    assert not benchmark.Comparison(timing, timing).uses_an_index
    assert benchmark.Comparison(timing, indexed).uses_an_index


# --- the shipped database ----------------------------------------------------


@pytest.mark.skipif(
    not (RESULTS_DIR / "quality_report.json").exists(), reason="run `make build` first"
)
def test_shipped_quality_report_passed():
    import json

    report = json.loads((RESULTS_DIR / "quality_report.json").read_text())
    assert report["all_passed"], report["failed"]
    assert report["counts"]["gene"] > 700


@pytest.mark.skipif(not (RESULTS_DIR / "Q5.tsv").exists(), reason="run `make queries` first")
def test_all_three_go_namespaces_are_present():
    """The original result file contained only BP; all three should be loaded."""
    with (RESULTS_DIR / "Q5.tsv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["go_namespace"] for row in rows} == {"BP", "MF", "CC"}


# --- CLI argument handling ---------------------------------------------------


def test_shared_options_are_accepted_after_the_subcommand():
    """`genomedb build --results-dir out` must parse.

    Declaring the shared options only on the top-level parser makes that an
    error and forces `genomedb --results-dir out build`, which is not the order
    anyone reaches for. A CI run caught exactly that.
    """
    from genomedb.cli import build_parser

    parser = build_parser()
    for command in ("build", "query", "check", "benchmark", "gene X", "go GO:1"):
        argv = [*command.split(), "--results-dir", "out", "--sql-dir", "sql"]
        args = parser.parse_args(argv)
        assert str(args.results_dir) == "out"
        assert str(args.sql_dir) == "sql"


def test_every_subcommand_binds_a_handler():
    from genomedb.cli import build_parser

    parser = build_parser()
    for command in ("build", "query", "check", "benchmark", "gene X", "go GO:1"):
        args = parser.parse_args(command.split())
        assert callable(args.func), command
