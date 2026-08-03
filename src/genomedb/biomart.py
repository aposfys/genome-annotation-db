"""Regenerating the Ensembl exports, so the data has verifiable provenance.

The exports in ``data/`` were originally produced by hand through the BioMart
web interface, which makes them a snapshot nobody can check: if the numbers a
query produces disagree with a previous run, there is no way to tell whether the
code changed or the download did.

This module holds the exact queries instead. Each export can be regenerated, and
because the Ensembl release is pinned in the URL, regenerating it tomorrow gives
the same rows it gave today.

That matters here for a concrete reason. An earlier result file for the GO
namespace query reports 769 genes and 13,939 annotations; the export shipped
alongside it yields 757 and 13,437 when the namespace is ignored, and 12,596
when it is respected. Ignoring the namespace reproduces the shape of the old
result, and the residual gap is consistent with the original having been loaded
from a slightly larger download. Neither half of that could be established
without knowing which query produced which file.
"""

from __future__ import annotations

import gzip
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

# Pinning the release is the point. www.ensembl.org always serves the current
# release, so a query run against it is not reproducible across releases.
ENSEMBL_RELEASE = 113
BIOMART_URL = f"https://may{ENSEMBL_RELEASE % 100}.archive.ensembl.org/biomart/martservice"
CURRENT_BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

# The chromosomes the relational database covers.
STUDY_CHROMOSOMES = "20,21"
ALL_CHROMOSOMES = ",".join([*(str(n) for n in range(1, 23)), "X", "Y", "MT"])


def _query(attributes: list[str], chromosomes: str) -> str:
    attribute_xml = "".join(f'<Attribute name="{a}"/>' for a in attributes)
    return (
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="1"'
        ' uniqueRows="1" datasetConfigVersion="0.6">'
        '<Dataset name="hsapiens_gene_ensembl" interface="default">'
        f'<Filter name="chromosome_name" value="{chromosomes}"/>'
        f"{attribute_xml}"
        "</Dataset></Query>"
    )


@dataclass(frozen=True)
class Export:
    """One BioMart download, and where it lands."""

    name: str
    attributes: tuple[str, ...]
    chromosomes: str
    description: str

    @property
    def query(self) -> str:
        return _query(list(self.attributes), self.chromosomes)

    def destination(self, data_dir: Path) -> Path:
        return data_dir / f"{self.name}.tsv.gz"


GENE_ATTRIBUTES = (
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype",
)

EXPORTS: tuple[Export, ...] = (
    Export(
        "genes",
        GENE_ATTRIBUTES,
        STUDY_CHROMOSOMES,
        "Genes on chromosomes 20 and 21 -- the relational database",
    ),
    Export(
        "transcripts_exons",
        (
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "transcript_start",
            "transcript_end",
            "transcript_length",
            "ensembl_exon_id",
            "rank",
            "exon_chrom_start",
            "exon_chrom_end",
        ),
        STUDY_CHROMOSOMES,
        "Transcript and exon structure for the same genes",
    ),
    Export(
        "go",
        ("ensembl_gene_id", "go_id", "name_1006", "definition_1006", "namespace_1003"),
        STUDY_CHROMOSOMES,
        "Gene Ontology annotation for the same genes",
    ),
    Export(
        "genes_all",
        GENE_ATTRIBUTES,
        ALL_CHROMOSOMES,
        "Every gene in the genome -- real coordinates for the scaling study",
    ),
)

BY_NAME = {export.name: export for export in EXPORTS}


def fetch(
    export: Export, data_dir: Path, use_archive: bool = True, timeout: int = 600
) -> Path:
    """Download one export and store it gzipped.

    Args:
        use_archive: Query the pinned release archive rather than the current
            release. Turn it off only to deliberately refresh against the
            latest Ensembl, which will change the row counts.
    """
    destination = export.destination(data_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)

    url = BIOMART_URL if use_archive else CURRENT_BIOMART_URL
    payload = urllib.parse.urlencode({"query": export.query}).encode()

    with urllib.request.urlopen(url, data=payload, timeout=timeout) as response:
        body = response.read()

    if not body.startswith(b"Gene stable ID"):
        head = body[:200].decode("utf-8", "replace")
        raise RuntimeError(f"BioMart returned something unexpected: {head!r}")

    with gzip.open(destination, "wb") as handle:
        handle.write(body)
    return destination


def row_count(path: Path) -> int:
    """Data rows in an export, excluding the header."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return max(sum(1 for _ in handle) - 1, 0)
    with path.open(encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)
