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

from . import benchmark, load, quality, queries
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genomedb",
        description=(
            "A relational genome annotation database over Ensembl BioMart "
            "exports: build it, query it, check it, and measure it."
        ),
    )
    parser.add_argument(
        "--db-url",
        default=None,
        help=(
            "sqlite:///path or mysql://user:password@host/database. "
            f"Defaults to ${DB_URL_VAR}, then to a local SQLite file."
        ),
    )
    parser.add_argument("--sqlite-path", type=Path, default=None)
    parser.add_argument("--sql-dir", type=Path, default=DEFAULT_SQL_DIR)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)

    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build", help="create the schema and load the data")
    build.add_argument(
        "--no-indexes",
        action="store_true",
        help="load without secondary indexes, to measure their effect",
    )
    build.set_defaults(func=cmd_build)

    query = sub.add_parser("query", help="run one query, or all of them")
    query.add_argument("name", nargs="?", help="Q1..Q7; omit to run every query")
    query.add_argument("--params", nargs="*", help="parameters for the query")
    query.add_argument("--limit", type=int, default=20, help="rows to print")
    query.add_argument("--save", action="store_true", help="also write results/<name>.tsv")
    query.set_defaults(func=cmd_query)

    gene = sub.add_parser("gene", help="transcripts of a gene, by symbol")
    gene.add_argument("gene_name")
    gene.add_argument("--limit", type=int, default=50)
    gene.set_defaults(func=cmd_gene)

    go = sub.add_parser("go", help="genes carrying a GO term")
    go.add_argument("go_id")
    go.add_argument("--limit", type=int, default=50)
    go.set_defaults(func=cmd_go)

    check = sub.add_parser("check", help="run the integrity checks")
    check.set_defaults(func=cmd_check)

    bench = sub.add_parser("benchmark", help="time the queries with and without indexes")
    bench.add_argument("--repeats", type=int, default=benchmark.REPEATS)
    bench.set_defaults(func=cmd_benchmark)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
