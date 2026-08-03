"""Command-line interface to the genome annotation database.

Subcommands rather than an interactive menu, so every operation is scriptable
and can be driven from a Makefile or CI. Connection details come from
``$GENOMEDB_URL`` or ``--db-url``; nothing is hardcoded.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from . import (
    benchmark,
    biomart,
    external,
    intervals,
    load,
    normalisation,
    quality,
    queries,
    scaling,
)
from .db import DB_URL_VAR, connection

DEFAULT_SQL_DIR = Path("sql")
DEFAULT_DATA_DIR = Path("data")
DEFAULT_RESULTS_DIR = Path("results")


def _print_table(
    columns: Sequence[str], rows: Sequence[tuple], limit: int | None = None
) -> None:
    shown = rows[:limit] if limit else rows
    if columns:
        print("\t".join(columns))
    for row in shown:
        print("\t".join("" if value is None else str(value) for value in row))
    if limit and len(rows) > limit:
        print(f"... {len(rows) - limit:,} more rows")


def _write_json(payload: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def cmd_build(args: argparse.Namespace) -> int:
    """Create the schema and load every table."""
    with connection(args.db_url, args.sqlite_path) as conn:
        print(f"Building on {conn.backend.name}...")
        report = load.build(
            conn, args.data_dir, args.sql_dir, with_indexes=not args.no_indexes
        )
        print(report.render())
        _write_json(report.as_dict(), args.results_dir / "load_report.json")

        checks = quality.report(conn)
        print("\nIntegrity checks:", "all passed" if checks["all_passed"] else "FAILED")
        if not checks["all_passed"]:
            for name in checks["failed"]:
                print(f"  failed: {name}")
        _write_json(checks, args.results_dir / "quality_report.json")
        return 0 if checks["all_passed"] else 1


def cmd_query(args: argparse.Namespace) -> int:
    """Run one query, or all of them."""
    with connection(args.db_url, args.sqlite_path) as conn:
        if args.name:
            query = queries.BY_NAME.get(args.name.upper())
            if query is None:
                print(
                    f"Unknown query {args.name!r}; expected one of "
                    f"{', '.join(queries.BY_NAME)}",
                    file=sys.stderr,
                )
                return 2
            columns, rows = queries.run(conn, query, args.sql_dir, args.params or None)
            print(f"# {query.name}: {query.title}  ({len(rows):,} rows)")
            _print_table(columns, rows, args.limit)
            if args.save:
                queries.write_tsv(columns, rows, args.results_dir / f"{query.name}.tsv")
            return 0

        summary = queries.run_all(conn, args.sql_dir, args.results_dir)
        for name, detail in summary.items():
            print(f"{name}: {detail['rows']:>6,} rows   {detail['title']}")
        _write_json(summary, args.results_dir / "query_summary.json")
        print(f"\nResults written to {args.results_dir}/")
        return 0


def cmd_gene(args: argparse.Namespace) -> int:
    """Look up the transcripts of a gene by symbol."""
    with connection(args.db_url, args.sqlite_path) as conn:
        columns, rows = queries.gene(conn, args.gene_name)
        if not rows:
            print(f"No transcripts found for gene {args.gene_name!r}")
            return 1
        _print_table(columns, rows, args.limit)
        return 0


def cmd_go(args: argparse.Namespace) -> int:
    """Look up the genes carrying a GO term."""
    with connection(args.db_url, args.sqlite_path) as conn:
        columns, rows = queries.go(conn, args.go_id)
        if not rows:
            print(f"No genes found for GO term {args.go_id!r}")
            return 1
        _print_table(columns, rows, args.limit)
        return 0


def cmd_check(args: argparse.Namespace) -> int:
    """Run the integrity checks."""
    with connection(args.db_url, args.sqlite_path) as conn:
        report = quality.report(conn)
        print("Table counts:")
        for table, count in report["counts"].items():
            print(f"  {table:18s} {count:>9,}")
        print("\nChecks:")
        for check in report["checks"]:
            status = "ok  " if check["passed"] else "FAIL"
            print(
                f"  {status} {check['check']:32s} {check['observed']:>6}"
                f"   {check['description']}"
            )
        _write_json(report, args.results_dir / "quality_report.json")
        return 0 if report["all_passed"] else 1


def cmd_benchmark(args: argparse.Namespace) -> int:
    """Time every query with and without the secondary indexes."""
    with connection(args.db_url, args.sqlite_path) as conn:
        comparisons = benchmark.compare(conn, args.sql_dir, repeats=args.repeats)
        summary = benchmark.summarise(comparisons)

        print()
        print(benchmark.render(comparisons))
        print(
            f"Overall: {summary['total_without_indexes_ms']:,.0f} ms without indexes, "
            f"{summary['total_with_indexes_ms']:,.0f} ms with "
            f"({summary['overall_speedup']}x)"
        )

        _write_json(
            {"summary": summary, "queries": [c.as_dict() for c in comparisons]},
            args.results_dir / "benchmark.json",
        )
        (args.results_dir / "benchmark.md").write_text(
            benchmark.render(comparisons), encoding="utf-8"
        )
        return 0


def _common_options() -> argparse.ArgumentParser:
    """Options every subcommand accepts.

    Attached to each subparser rather than to the top-level parser, so they can
    be written after the subcommand -- `genomedb build --results-dir out` --
    which is the order people reach for. Declaring them only at the top level
    makes that form an error.
    """
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--db-url",
        default=None,
        help=(
            "sqlite:///path or mysql://user:password@host/database. "
            f"Defaults to ${DB_URL_VAR}, then to a local SQLite file."
        ),
    )
    common.add_argument("--sqlite-path", type=Path, default=None)
    common.add_argument("--sql-dir", type=Path, default=DEFAULT_SQL_DIR)
    common.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    common.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    return common


def cmd_region(args: argparse.Namespace) -> int:
    """Genes overlapping a genomic window."""
    with connection(args.db_url, args.sqlite_path) as conn:
        if args.strategy != "btree":
            intervals.build(conn, args.strategy)
        rows = intervals.query(conn, args.strategy, args.chrom, args.start, args.end)
        print(
            f"# {args.chrom}:{args.start:,}-{args.end:,}  "
            f"({len(rows)} genes, via {args.strategy})"
        )
        _print_table(["gene_id", "gene_name", "gene_start", "gene_end"], rows, args.limit)
        return 0


def cmd_intervals(args: argparse.Namespace) -> int:
    """Benchmark the three interval-indexing strategies against each other."""
    with connection(args.db_url, args.sqlite_path) as conn:
        results, summary = benchmark.compare_intervals(
            conn, window_count=args.windows, width=args.width, repeats=args.repeats
        )
        print()
        print(benchmark.render_intervals(results))
        if not summary["strategies_agree"]:
            print("WARNING: strategies disagree:", summary["disagreements"])
        print(
            f"Fastest: {summary['fastest']}  "
            f"(speed-up over B-tree: {summary['speedup_over_btree']})"
        )

        _write_json(
            {"summary": summary, "strategies": [r.as_dict() for r in results]},
            args.results_dir / "intervals.json",
        )
        (args.results_dir / "intervals.md").write_text(
            benchmark.render_intervals(results), encoding="utf-8"
        )
        return 0 if summary["strategies_agree"] else 1


def cmd_scaling(args: argparse.Namespace) -> int:
    """Measure how query cost grows with table size, and fit the growth law."""
    requested = tuple(args.sizes) if args.sizes else scaling.DEFAULT_SIZES
    source = None if args.synthetic else args.data_dir
    sizes, available = scaling.usable_sizes(requested, source)
    if available is not None and sizes != list(requested):
        print(
            f"The genome supplies {available:,} genes; measuring at {sizes}"
            f" rather than {list(requested)}"
        )
    with connection(args.db_url, args.sqlite_path) as conn:
        origin = "synthetic intervals" if args.synthetic else "real gene coordinates"
        print(f"Point lookup, with and without an index ({origin}):")
        point = scaling.measure_point_lookups(conn, sizes, probes=args.probes, data_dir=source)
        ug, ig = point["unindexed_growth"], point["indexed_growth"]
        print(
            f"  unindexed: {ug['verdict']}"
            f"  (linear R2 {ug['linear']['r_squared']:.3f},"
            f" grew {ug['growth_factor']}x)"
        )
        print(
            f"  indexed:   {ig['verdict']}"
            f"  (grew only {ig['growth_factor']}x over the same range)"
        )

        print("\nInterval overlap, three structures:")
        interval = scaling.measure_interval_strategies(
            conn, sizes, probes=args.probes, data_dir=source
        )
        print(f"  crossover at n = {interval['crossover_n']}")

        _write_json(
            {"point_lookup": point, "intervals": interval},
            args.results_dir / "scaling.json",
        )
        (args.results_dir / "scaling.md").write_text(
            scaling.render(point, interval), encoding="utf-8"
        )
        print(f"\nWritten to {args.results_dir}/scaling.json")
        return 0


def cmd_validate(args: argparse.Namespace) -> int:
    """Check the overlap results against bedtools, an independent implementation."""
    with connection(args.db_url, args.sqlite_path) as conn:
        intervals.build(conn, args.strategy)

        random_windows = list(intervals.windows(conn, args.windows))
        boundary = external.boundary_windows(conn, limit=args.boundary_genes)

        report = {}
        for label, windows in (("random", random_windows), ("boundary", boundary)):
            comparison = external.validate(conn, windows, strategy=args.strategy)
            report[label] = comparison.as_dict()
            status = "agree" if comparison.agrees else "DISAGREE"
            print(f"{label:9s} windows: {comparison.agreed}/{comparison.windows} {status}")
            if not comparison.agrees:
                for window, genes in list(comparison.only_bedtools.items())[:3]:
                    print(f"    bedtools found {genes} at {window}, we did not")
                for window, genes in list(comparison.only_ours.items())[:3]:
                    print(f"    we found {genes} at {window}, bedtools did not")

        _write_json(report, args.results_dir / "bedtools_validation.json")
        agreed = all(section["agrees"] for section in report.values())
        print("\nbedtools agreement:", "complete" if agreed else "INCOMPLETE")
        return 0 if agreed else 1


def cmd_normalisation(args: argparse.Namespace) -> int:
    """Check BCNF, verify the dependencies in the data, and price the trade-off."""
    bcnf = normalisation.check()
    print("Boyce-Codd Normal Form:")
    for name, detail in bcnf["relations"].items():
        keys = " | ".join("(" + ", ".join(k) + ")" for k in detail["candidate_keys"])
        status = "BCNF" if detail["in_bcnf"] else "VIOLATION"
        print(f"  {status:9s} {name:16s} key {keys}")
    if not bcnf["all_in_bcnf"]:
        print("  violating:", bcnf["violating_relations"])

    with connection(args.db_url, args.sqlite_path) as conn:
        holds = normalisation.verify_against_data(conn)
        print(
            f"\nDependencies verified against the loaded data: "
            f"{holds['checked']} checked, "
            f"{'all hold' if holds['all_hold'] else 'VIOLATIONS: ' + str(holds['violations'])}"
        )

        denorm = normalisation.measure_denormalisation(conn, repeats=args.repeats)
        print(f"\nWhat normalisation costs, on {denorm['rows']:,} transcripts:")
        print(f"  join every time      {denorm['normalised_ms']:>8.3f} ms")
        print(
            f"  materialised column  {denorm['denormalised_ms']:>8.3f} ms"
            f"   ({denorm['speedup']}x faster)"
        )
        print(f"  {denorm['trade_off']}")

        _write_json(
            {"bcnf": bcnf, "dependencies_in_data": holds, "denormalisation": denorm},
            args.results_dir / "normalisation.json",
        )
        (args.results_dir / "normalisation.md").write_text(
            normalisation.render(bcnf, denorm), encoding="utf-8"
        )
        return 0 if bcnf["all_in_bcnf"] and holds["all_hold"] else 1


def cmd_fetch(args: argparse.Namespace) -> int:
    """Regenerate the Ensembl exports from BioMart."""
    wanted = args.export or list(biomart.BY_NAME)
    print(
        f"Ensembl release {biomart.ENSEMBL_RELEASE}"
        f" ({'archive, pinned' if not args.current else 'current, NOT pinned'})"
    )
    for name in wanted:
        export = biomart.BY_NAME.get(name)
        if export is None:
            print(
                f"Unknown export {name!r}; expected one of {list(biomart.BY_NAME)}",
                file=sys.stderr,
            )
            return 2
        print(f"  fetching {name}: {export.description}...")
        path = biomart.fetch(export, args.data_dir, use_archive=not args.current)
        print(f"    {path.name}  {biomart.row_count(path):,} rows")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genomedb",
        description=(
            "A relational genome annotation database over Ensembl BioMart "
            "exports: build it, query it, check it, and measure it."
        ),
    )
    common = _common_options()
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser(
        "build", parents=[common], help="create the schema and load the data"
    )
    build.add_argument(
        "--no-indexes",
        action="store_true",
        help="load without secondary indexes, to measure their effect",
    )
    build.set_defaults(func=cmd_build)

    query = sub.add_parser("query", parents=[common], help="run one query, or all of them")
    query.add_argument("name", nargs="?", help="Q1..Q7; omit to run every query")
    query.add_argument("--params", nargs="*", help="parameters for the query")
    query.add_argument("--limit", type=int, default=20, help="rows to print")
    query.add_argument("--save", action="store_true", help="also write results/<name>.tsv")
    query.set_defaults(func=cmd_query)

    gene = sub.add_parser("gene", parents=[common], help="transcripts of a gene, by symbol")
    gene.add_argument("gene_name")
    gene.add_argument("--limit", type=int, default=50)
    gene.set_defaults(func=cmd_gene)

    go = sub.add_parser("go", parents=[common], help="genes carrying a GO term")
    go.add_argument("go_id")
    go.add_argument("--limit", type=int, default=50)
    go.set_defaults(func=cmd_go)

    check = sub.add_parser("check", parents=[common], help="run the integrity checks")
    check.set_defaults(func=cmd_check)

    bench = sub.add_parser(
        "benchmark", parents=[common], help="time the queries with and without indexes"
    )
    bench.add_argument("--repeats", type=int, default=benchmark.REPEATS)
    bench.set_defaults(func=cmd_benchmark)

    region = sub.add_parser(
        "region", parents=[common], help="genes overlapping a genomic window"
    )
    region.add_argument("chrom")
    region.add_argument("start", type=int)
    region.add_argument("end", type=int)
    region.add_argument(
        "--strategy",
        choices=sorted(intervals.BY_NAME),
        default="btree",
        help="which interval index to answer through. Default: btree.",
    )
    region.add_argument("--limit", type=int, default=50)
    region.set_defaults(func=cmd_region)

    iv = sub.add_parser(
        "intervals",
        parents=[common],
        help="benchmark B-tree vs UCSC binning vs R*Tree on overlap queries",
    )
    iv.add_argument("--windows", type=int, default=200)
    iv.add_argument("--width", type=int, default=1_000_000, help="window size in bp")
    iv.add_argument("--repeats", type=int, default=3)
    iv.set_defaults(func=cmd_intervals)

    sc = sub.add_parser(
        "scaling",
        parents=[common],
        help="measure how query cost grows with table size",
    )
    sc.add_argument("--sizes", type=int, nargs="*", help="table sizes to measure")
    sc.add_argument("--probes", type=int, default=200, help="queries per size")
    sc.add_argument(
        "--synthetic",
        action="store_true",
        help="generate intervals instead of using real gene coordinates",
    )
    sc.set_defaults(func=cmd_scaling)

    va = sub.add_parser(
        "validate",
        parents=[common],
        help="check overlap results against bedtools",
    )
    va.add_argument("--windows", type=int, default=200)
    va.add_argument("--boundary-genes", type=int, default=100)
    va.add_argument("--strategy", choices=sorted(intervals.BY_NAME), default="btree")
    va.set_defaults(func=cmd_validate)

    nf = sub.add_parser(
        "normalisation",
        parents=[common],
        help="check BCNF and measure what normalisation costs",
    )
    nf.add_argument("--repeats", type=int, default=20)
    nf.set_defaults(func=cmd_normalisation)

    fe = sub.add_parser(
        "fetch",
        parents=[common],
        help="regenerate the Ensembl exports from BioMart",
    )
    fe.add_argument("--export", nargs="*", choices=sorted(biomart.BY_NAME))
    fe.add_argument(
        "--current",
        action="store_true",
        help="query the current Ensembl release rather than the pinned archive",
    )
    fe.set_defaults(func=cmd_fetch)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
