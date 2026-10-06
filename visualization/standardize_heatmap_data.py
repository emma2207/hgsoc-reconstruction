#!/usr/bin/env python3
"""Export parsed heatmap matrices as consistently labeled CSV files.

Example real-data export:
    python standardize_heatmap_data.py \
        --data-path ../benchmarking/data/output_data/vcf_update \
        --output-dir ../benchmarking/data/standardized_heatmaps \
        --dataset hgsoc \
        --mod1 single-cell \
        --mod2 bulk_chunk_ribo \
        --read-depth single-cell_0_bulk_chunk_ribo_0 \
        --tools HYSYS CrosscheckFingerprints Conpair

Example pseudobulk export:
    python standardize_heatmap_data.py \
        --data-path ../benchmarking/data/output_data/vcf_update \
        --output-dir ../benchmarking/data/standardized_heatmaps \
        --dataset high_grade_glioma \
        --pseudobulk \
        --ncells 50000 \
        --read-depth 0 \
        --mod1 _1 \
        --mod2 _2 \
        --tools HYSYS Vireo
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Sequence

import pandas as pd

from analysis_shared_functions import load_heatmap_data


TOOLS = (
    "Conpair",
    "CrosscheckFingerprints",
    "HYSYS",
    "NGSCheckmate",
    "ntsm",
    "OmicsPrint",
    "Peddy",
    "Somalier",
    "TimeAttackGenComp",
    "Vireo",
)

# These result sources do not have a read-depth directory in the output layout.
READ_DEPTH_INDEPENDENT_TOOLS = {"BAMixChecker", "Conpair", "Somalier", "ntsm"}


def _read_depth_component(read_depth: str) -> str:
    """Normalize a read-depth selector for use as one directory component."""
    component = str(read_depth)
    if component.startswith("read_depth_"):
        component = component.removeprefix("read_depth_")
    if not component or component in {".", ".."} or "/" in component or "\\" in component:
        raise ValueError(f"Invalid read-depth directory component: {read_depth!r}")
    return f"read_depth_{component}"


def _matrix_output_directory(
    output_dir: Path,
    dataset: str,
    pseudobulk: bool,
    ncells: int | None,
    read_depth: str,
    mod1: str,
    mod2: str,
    tool: str,
) -> Path:
    if pseudobulk:
        if ncells is None:
            raise ValueError("ncells is required for pseudobulk exports")
        relative_path = Path(dataset) / "pseudobulk" / f"ncells_{ncells}"
    else:
        relative_path = Path(dataset) / f"{mod1}_vs_{mod2}" / "ncells_null"

    if tool not in READ_DEPTH_INDEPENDENT_TOOLS:
        relative_path /= _read_depth_component(read_depth)

    return output_dir / relative_path


def _standardize_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    """Normalize labels and axis ordering without changing matrix values."""
    standardized = matrix.copy()
    standardized.index = standardized.index.map(str)
    standardized.columns = standardized.columns.map(str)

    if standardized.index.has_duplicates:
        raise ValueError("Matrix has duplicate row sample labels")
    if standardized.columns.has_duplicates:
        raise ValueError("Matrix has duplicate column sample labels")

    return standardized.sort_index(axis=0).sort_index(axis=1)


def _deduplicate_modality_label(label: object, modalities: Sequence[str]) -> str:
    """Collapse adjacent repeated copies of either requested modality label."""
    standardized = str(label)
    for modality in modalities:
        variants = list(dict.fromkeys((modality, modality.replace("-", "_"), modality.replace("_", "-"))))
        for first in variants:
            for second in variants:
                repeated = f"_{first}_{second}"
                single = f"_{first}"
                while repeated in standardized:
                    standardized = standardized.replace(repeated, single)
    return standardized


def export_heatmap_matrices(
    data_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    dataset: str,
    tools: Sequence[str],
    read_depth: str,
    mod1: str,
    mod2: str,
    *,
    pseudobulk: bool = False,
    ncells: int | None = None,
) -> list[Path]:
    """Load, standardize, and export selected tool matrices as labeled CSV files.

    Each CSV has sample IDs as both row and column labels. Missing values are
    retained, and row/column labels are sorted independently for stable exports.
    Returns the paths written; parsers that report missing source files are skipped.
    """
    data_path = Path(data_path)
    output_dir = Path(output_dir)
    written_paths = []

    for tool in tools:
        matrix = load_heatmap_data(
            str(data_path),
            pseudobulk,
            tool,
            dataset,
            ncells,
            read_depth,
            mod1,
            mod2,
        )
        if matrix is None:
            print(f"Skipping {tool}: no source matrix available")
            continue

        matrix.index = matrix.index.map(
            lambda label: _deduplicate_modality_label(label, (mod1, mod2))
        )
        matrix.columns = matrix.columns.map(
            lambda label: _deduplicate_modality_label(label, (mod1, mod2))
        )

        # Treat hyphens and underscores equivalently in modality labels.
        mod1_suffix = mod1.replace("-", "_")
        mod2_suffix = mod2.replace("-", "_")
        row_mask = matrix.index.to_series().map(
            lambda label: str(label).replace("-", "_").endswith(mod1_suffix)
        )
        column_mask = matrix.columns.to_series().map(
            lambda label: str(label).replace("-", "_").endswith(mod2_suffix)
        )
        matrix = matrix.loc[row_mask.to_numpy(), :]
        matrix = matrix.loc[:, column_mask.to_numpy()]

        matrix = _standardize_matrix(matrix)
        tool_output_dir = _matrix_output_directory(
            output_dir,
            dataset,
            pseudobulk,
            ncells,
            read_depth,
            mod1,
            mod2,
            tool,
        )
        tool_output_dir.mkdir(parents=True, exist_ok=True)
        output_path = tool_output_dir / f"{tool}.csv"
        matrix.to_csv(output_path, index_label="sample_id")
        written_paths.append(output_path)
        print(f"Wrote {output_path} ({matrix.shape[0]} x {matrix.shape[1]})")

    return written_paths


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export parsed heatmap matrices as consistently labeled CSV files."
    )
    parser.add_argument("--data-path", required=True, help="Root directory containing tool outputs")
    parser.add_argument("--output-dir", required=True, help="Directory for standardized CSV matrices")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--tools", nargs="+", required=True, choices=TOOLS)
    parser.add_argument("--read-depth", required=True, help="Read-depth tag, e.g. bulk_0_single-cell_0")
    parser.add_argument("--mod1", required=True)
    parser.add_argument("--mod2", required=True)
    parser.add_argument("--pseudobulk", action="store_true")
    parser.add_argument("--ncells", type=int, help="Required when --pseudobulk is set")
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    if args.pseudobulk and args.ncells is None:
        raise SystemExit("--ncells is required when --pseudobulk is set")

    written_paths = export_heatmap_matrices(
        args.data_path,
        args.output_dir,
        args.dataset,
        args.tools,
        args.read_depth,
        args.mod1,
        args.mod2,
        pseudobulk=args.pseudobulk,
        ncells=args.ncells,
    )
    return 0 if written_paths else 1


if __name__ == "__main__":
    raise SystemExit(main())
