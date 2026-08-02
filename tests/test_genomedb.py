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


# --- interval structures -----------------------------------------------------


def test_bin_assignment_contains_the_interval():
    """An interval's own bin must be among the bins a query for it searches."""
    from genomedb import intervals

    for start, end in [
        (0, 1000),
        (100_000, 200_000),
        (1_000_000, 1_000_100),
        (0, 200_000_000),
    ]:
        assigned = intervals.assign_bin(start, end)
        assert assigned in intervals.overlapping_bins(start, end)


def test_larger_intervals_land_in_coarser_bins():
    """The scheme's whole point: a feature too big for a fine bin moves up."""
    from genomedb import intervals

    small = intervals.assign_bin(1_000_000, 1_000_100)
    huge = intervals.assign_bin(0, 200_000_000)
    assert huge < small  # coarser levels have lower bin numbers


def test_a_query_touches_few_bins():
    from genomedb import intervals

    assert len(intervals.overlapping_bins(1_000_000, 2_000_000)) < 20


@pytest.fixture
def interval_db(tiny_db):
    from genomedb import intervals

    conn, _ = tiny_db
    for strategy in intervals.STRATEGIES:
        intervals.build(conn, strategy.name)
    return conn


def test_all_strategies_return_identical_results(interval_db):
    """A faster structure that returns different rows is not a faster structure."""
    from genomedb import intervals

    windows = [("20", 0, 1000), ("20", 150, 450), ("21", 0, 1000), ("20", 5000, 6000)]
    for chrom, start, end in windows:
        results = {}
        for strategy in intervals.STRATEGIES:
            rows = intervals.query(interval_db, strategy.name, chrom, start, end)
            results[strategy.name] = tuple(row[0] for row in rows)
        assert len(set(results.values())) == 1, (chrom, start, end, results)


# --- scaling -----------------------------------------------------------------


def test_a_flat_response_is_not_classified_as_a_growth_law():
    """Fitting a law to noise is the failure mode this guards against."""
    from genomedb import scaling

    sizes = [1_000, 4_000, 16_000, 64_000, 256_000]
    flat = [7.0, 7.3, 7.2, 7.5, 8.0]
    assert scaling.classify(sizes, flat)["verdict"] == "O(1)"


def test_linear_growth_is_recognised():
    from genomedb import scaling

    sizes = [1_000, 4_000, 16_000, 64_000, 256_000]
    linear = [float(s) for s in sizes]
    result = scaling.classify(sizes, linear)
    assert result["verdict"] == "O(n)"
    assert result["linear"]["r_squared"] > 0.99


def test_logarithmic_growth_is_recognised():
    import math

    from genomedb import scaling

    sizes = [1_000, 4_000, 16_000, 64_000, 256_000]
    logarithmic = [math.log2(s) for s in sizes]
    result = scaling.classify(sizes, logarithmic)
    assert result["verdict"] == "O(log n)"


def test_synthetic_intervals_are_valid_and_reproducible():
    from genomedb import scaling

    first = scaling.synthesise(500)
    assert first == scaling.synthesise(500)
    assert len(first) == 500
    for _identifier, _chrom, start, end in first:
        assert end > start
        assert 0 <= start < scaling.COORDINATE_SPAN


# --- coordinate conventions and external validation --------------------------


def test_bed_conversion_is_lossless():
    """1-based inclusive to 0-based half-open and back."""
    from genomedb import external

    for chrom, start, end in [("20", 87250, 97094), ("21", 1, 1), ("X", 500, 1500)]:
        assert external.from_bed(*external.to_bed(chrom, start, end)) == (chrom, start, end)


def test_bed_conversion_moves_only_the_start():
    """A 1-based inclusive [87250, 97094] is BED [87249, 97094)."""
    from genomedb import external

    assert external.to_bed("20", 87250, 97094) == ("20", 87249, 97094)


def test_overlap_is_inclusive_not_half_open(interval_db):
    """Ensembl coordinates are 1-based inclusive.

    A gene ending exactly where a window begins shares one base with it and
    does overlap. Treating the coordinates as half-open silently drops those,
    which is the classic genomics off-by-one.
    """
    from genomedb import intervals

    # G1 spans 100-900 inclusive on chromosome 20.
    assert intervals.query(interval_db, "btree", "20", 900, 1000) != []
    assert intervals.query(interval_db, "btree", "20", 901, 1000) == []
    assert intervals.query(interval_db, "btree", "20", 50, 100) != []
    assert intervals.query(interval_db, "btree", "20", 50, 99) == []


def test_all_strategies_share_the_inclusive_convention(interval_db):
    """An off-by-one in the shared predicate would pass the agreement test."""
    from genomedb import intervals

    for strategy in intervals.STRATEGIES:
        touching = intervals.query(interval_db, strategy.name, "20", 900, 1000)
        assert [r[0] for r in touching] == ["G1"], strategy.name


# --- normalisation -----------------------------------------------------------


def test_every_relation_is_in_bcnf():
    from genomedb import normalisation

    result = normalisation.check()
    assert result["all_in_bcnf"], result["violating_relations"]


def test_a_non_superkey_determinant_is_caught():
    """The check must be capable of failing, or it proves nothing."""
    from genomedb.normalisation import Dependency, Relation

    # chrom -> biotype is not implied by any key: a violation by construction.
    broken = Relation(
        name="broken",
        attributes=("gene_id", "chrom", "biotype"),
        candidate_keys=(("gene_id",),),
        dependencies=(Dependency("broken", ("chrom",), ("biotype",)),),
    )
    assert not broken.in_bcnf
    assert len(broken.violations()) == 1


def test_trivial_dependencies_do_not_violate_bcnf():
    from genomedb.normalisation import Dependency, Relation

    trivial = Relation(
        name="t",
        attributes=("a", "b"),
        candidate_keys=(("a",),),
        dependencies=(Dependency("t", ("a", "b"), ("b",)),),
    )
    assert trivial.in_bcnf


def test_exon_rank_depends_on_the_whole_composite_key():
    """The same exon can be ranked differently in different transcripts, so
    rank is a property of the pairing, not of the exon."""
    from genomedb import normalisation

    junction = next(r for r in normalisation.SCHEMA if r.name == "transcript_exon")
    assert junction.candidate_keys == (("transcript_id", "exon_id"),)
    assert junction.in_bcnf


def test_declared_dependencies_hold_in_the_data(tiny_db):
    from genomedb import normalisation

    conn, _ = tiny_db
    result = normalisation.verify_against_data(conn)
    assert result["all_hold"], result["violations"]


def test_denormalisation_measures_a_real_trade_off(tiny_db):
    from genomedb import normalisation

    conn, _ = tiny_db
    result = normalisation.measure_denormalisation(conn, repeats=3)
    assert result["rows"] > 0
    assert result["rows_disagreeing_now"] == 0
    assert result["normalised_ms"] > 0 and result["denormalised_ms"] > 0


@pytest.mark.skipif(
    __import__("shutil").which("bedtools") is None, reason="bedtools not installed"
)
def test_agrees_with_bedtools_on_boundary_windows(interval_db):
    """The check that caught the off-by-one.

    Boundary windows sit exactly on a gene's first and last base, which is
    where a coordinate-convention error has to show. Random windows agreed even
    while the bug was present.
    """
    from genomedb import external

    windows = external.boundary_windows(interval_db, limit=3)
    comparison = external.validate(interval_db, windows)
    assert comparison.agrees, (
        comparison.only_ours,
        comparison.only_bedtools,
    )


@pytest.mark.skipif(
    __import__("shutil").which("bedtools") is None, reason="bedtools not installed"
)
def test_the_bedtools_check_can_actually_fail(interval_db, monkeypatch):
    """A validation never seen to fail is not yet a validation.

    Reverting the inclusive comparison to the half-open one must make bedtools
    disagree; if it does not, the check has no power.
    """
    from genomedb import external, intervals

    monkeypatch.setattr(
        intervals,
        "BTREE_SQL",
        intervals.BTREE_SQL.replace("<= ?", "< ?").replace(">= ?", "> ?"),
    )
    windows = external.boundary_windows(interval_db, limit=3)
    assert not external.validate(interval_db, windows).agrees
