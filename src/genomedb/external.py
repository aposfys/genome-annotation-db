"""Validation against bedtools, an independent implementation.

The three interval strategies are checked against each other, which catches a
mistake in any one of them but not a mistake they share. A coordinate-convention
error is exactly the kind they would share: every strategy reads the same
columns and applies the same comparison, so all three would agree and all three
would be wrong.

bedtools is the standard toolkit for genomic interval arithmetic and an entirely
separate implementation. Agreeing with it is meaningful evidence in a way that
agreeing with oneself is not.

The conversion is where the risk sits. This database stores Ensembl's **1-based
inclusive** coordinates; BED is **0-based half-open**. A gene at 1-based
[87250, 97094] is BED [87249, 97094) — the start moves back one, the end does
not move. Getting that wrong shifts every comparison by a base and is invisible
until an external tool disagrees.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .db import Connection


def to_bed(chrom: str, start: int, end: int) -> tuple[str, int, int]:
    """Convert a 1-based inclusive interval to 0-based half-open BED."""
    return chrom, start - 1, end


def from_bed(chrom: str, start: int, end: int) -> tuple[str, int, int]:
    """Convert a 0-based half-open BED interval back to 1-based inclusive."""
    return chrom, start + 1, end


def _require() -> str:
    path = shutil.which("bedtools")
    if path is None:
        raise RuntimeError(
            "bedtools not found on PATH. Install it with: conda install -c bioconda bedtools"
        )
    return path


def export_genes(conn: Connection, path: Path) -> Path:
    """Write every gene as a BED3+ record, sorted as bedtools expects."""
    _, rows = conn.query(
        "SELECT chrom, gene_start, gene_end, gene_id FROM gene"
        " ORDER BY chrom, gene_start, gene_end"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for chrom, start, end, gene_id in rows:
            bed_chrom, bed_start, bed_end = to_bed(chrom, start, end)
            handle.write(f"{bed_chrom}\t{bed_start}\t{bed_end}\t{gene_id}\n")
    return path


def write_windows(windows: Sequence[tuple[str, int, int]], path: Path) -> Path:
    """Write query windows as BED, applying the same conversion."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for chrom, start, end in sorted(windows):
            bed_chrom, bed_start, bed_end = to_bed(chrom, start, end)
            handle.write(f"{bed_chrom}\t{bed_start}\t{bed_end}\n")
    return path


def intersect(genes_bed: Path, windows_bed: Path) -> dict[tuple[str, int, int], set[str]]:
    """Run ``bedtools intersect`` and group the gene hits by query window.

    Returns:
        A mapping from each window, in 1-based inclusive coordinates, to the
        set of gene identifiers bedtools reports as overlapping it.
    """
    result = subprocess.run(
        [
            _require(),
            "intersect",
            "-a",
            str(windows_bed),
            "-b",
            str(genes_bed),
            "-wa",
            "-wb",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    hits: dict[tuple[str, int, int], set[str]] = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        window = from_bed(fields[0], int(fields[1]), int(fields[2]))
        hits.setdefault(window, set()).add(fields[6])
    return hits


@dataclass(frozen=True)
class Comparison:
    """Agreement between this database and bedtools over a set of windows."""

    windows: int
    agreed: int
    only_ours: dict[str, list[str]]
    only_bedtools: dict[str, list[str]]

    @property
    def agrees(self) -> bool:
        return self.agreed == self.windows

    def as_dict(self) -> dict[str, Any]:
        return {
            "windows": self.windows,
            "agreed": self.agreed,
            "agreement_rate": round(self.agreed / self.windows, 4) if self.windows else 0.0,
            "windows_where_only_we_found_a_gene": self.only_ours,
            "windows_where_only_bedtools_found_a_gene": self.only_bedtools,
            "agrees": self.agrees,
        }


def validate(
    conn: Connection,
    windows: Sequence[tuple[str, int, int]],
    strategy: str = "btree",
    workspace: Path | None = None,
) -> Comparison:
    """Compare this database's overlap results against bedtools, window by window."""
    from . import intervals

    with tempfile.TemporaryDirectory() as tmp:
        directory = workspace or Path(tmp)
        genes_bed = export_genes(conn, directory / "genes.bed")
        windows_bed = write_windows(windows, directory / "windows.bed")
        theirs = intersect(genes_bed, windows_bed)

        agreed = 0
        only_ours: dict[str, list[str]] = {}
        only_theirs: dict[str, list[str]] = {}

        for window in windows:
            ours = {row[0] for row in intervals.query(conn, strategy, *window)}
            reference = theirs.get(window, set())
            if ours == reference:
                agreed += 1
                continue

            label = f"{window[0]}:{window[1]}-{window[2]}"
            if ours - reference:
                only_ours[label] = sorted(ours - reference)[:5]
            if reference - ours:
                only_theirs[label] = sorted(reference - ours)[:5]

    return Comparison(
        windows=len(windows),
        agreed=agreed,
        only_ours=only_ours,
        only_bedtools=only_theirs,
    )


def boundary_windows(conn: Connection, limit: int = 200) -> list[tuple[str, int, int]]:
    """Windows placed exactly on gene boundaries.

    Random windows almost never land on an edge, so they cannot detect an
    off-by-one. These deliberately do: for each gene, a window starting exactly
    at its last base and one ending exactly at its first.
    """
    _, rows = conn.query(
        "SELECT chrom, gene_start, gene_end FROM gene ORDER BY gene_id LIMIT ?",
        (limit,),
    )
    windows: list[tuple[str, int, int]] = []
    for chrom, start, end in rows:
        windows.append((chrom, end, end + 500))  # touches the final base
        windows.append((chrom, end + 1, end + 500))  # one past it
        windows.append((chrom, max(1, start - 500), start))  # touches the first base
        windows.append((chrom, max(1, start - 500), start - 1))
    return windows
