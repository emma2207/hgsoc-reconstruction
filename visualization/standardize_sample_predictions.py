#!/usr/bin/env python3
"""Export standardized binary sample-match prediction matrices.

Example real-data export:
    python standardize_sample_predictions.py \
        --data-path ../benchmarking/data/output_data/vcf_update \
        --output-dir ../benchmarking/data/standardized_predictions \
        --dataset hgsoc \
        --mod1 bulk_dissociated_ribo \
        --mod2 single-cell \
        --read-depth 0 \
        --tools Conpair CrosscheckFingerprints HYSYS NGSCheckmate ntsm

Example with modality-specific read depths:
    python standardize_sample_predictions.py \
        --data-path ../benchmarking/data/output_data/vcf_update \
        --output-dir ../benchmarking/data/standardized_predictions \
        --dataset hgsoc \
        --mod1 bulk_dissociated_ribo \
        --mod2 single-cell \
        --read-depth 0 --read-depth-mod2 10 \
        --tools Vireo

Example with a custom threshold:
    python standardize_sample_predictions.py \
        --data-path ../benchmarking/data/output_data/vcf_update \
        --output-dir ../benchmarking/data/standardized_predictions \
        --dataset high_grade_glioma \
        --pseudobulk --ncells 10000 \
        --mod1 _1 --mod2 _2 --read-depth 0 \
        --tools Conpair --threshold Conpair=0.7
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd

from analysis_shared_functions import load_sample_matching_results
from standardize_heatmap_data import (
    _deduplicate_modality_label,
    _matrix_output_directory,
    _standardize_matrix,
)


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

PREDICTION_THRESHOLDS = {
    "Conpair": ("ge", 0.8),
    "OmicsPrint": ("ge", 1.7),
    "Peddy": ("le", 0.0),
    "Somalier": ("ge", 1.0),
    "TimeAttackGenComp": ("le", 0.05),
}


def _modality_mask(labels: pd.Index, modality: str) -> list[bool]:
    normalized_modality = modality.replace("-", "_")
    return [
        str(label).replace("-", "_").endswith(normalized_modality)
        for label in labels
    ]


def _parse_threshold_overrides(
    specifications: Sequence[str], selected_tools: Sequence[str]
) -> dict[str, float]:
    overrides = {}
    for specification in specifications:
        try:
            tool, raw_value = specification.split("=", maxsplit=1)
            value = float(raw_value)
        except ValueError as exc:
            raise ValueError(
                f"Invalid threshold {specification!r}; expected TOOL=VALUE"
            ) from exc
        if tool not in PREDICTION_THRESHOLDS:
            raise ValueError(
                f"{tool!r} does not use a configurable threshold; threshold tools are "
                f"{', '.join(PREDICTION_THRESHOLDS)}"
            )
        if tool not in selected_tools:
            raise ValueError(f"Threshold supplied for {tool}, but it is not selected in --tools")
        if not math.isfinite(value):
            raise ValueError(f"Threshold for {tool} must be a finite number")
        if tool in overrides:
            raise ValueError(f"Threshold for {tool} was specified more than once")
        overrides[tool] = value
    return overrides


def _format_threshold(value: float) -> str:
    return format(value, "g")


def export_sample_predictions(
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
    read_depth_mod2: str | None = None,
    threshold_overrides: Mapping[str, float] | None = None,
    standardized_heatmaps_path: str | os.PathLike[str] | None = None,
) -> list[Path]:
    """Load, standardize, and export match predictions as labeled CSV files.

    Values are 1 for predicted matches, 0 for predicted non-matches, and NaN
    where the source reports an unknown/inconclusive result. NTSM's existing
    loader infers non-matches as zero when a pair is absent from its match list.
    """
    data_path = Path(data_path)
    output_dir = Path(output_dir)
    written_paths: list[Path] = []
    threshold_overrides = dict(threshold_overrides or {})

    unknown_overrides = set(threshold_overrides) - set(PREDICTION_THRESHOLDS)
    if unknown_overrides:
        raise ValueError(
            "Threshold overrides are not supported for: "
            + ", ".join(sorted(unknown_overrides))
        )
    unselected_overrides = set(threshold_overrides) - set(tools)
    if unselected_overrides:
        raise ValueError(
            "Threshold supplied for unselected tool(s): "
            + ", ".join(sorted(unselected_overrides))
        )

    for tool in tools:
        threshold = threshold_overrides.get(tool)
        matrix = load_sample_matching_results(
            str(data_path),
            pseudobulk,
            tool,
            dataset,
            ncells,
            read_depth,
            read_depth_mod2,
            mod1,
            mod2,
            threshold=threshold,
            standardized_heatmaps_path=standardized_heatmaps_path,
        )
        if matrix is None:
            print(f"Skipping {tool}: no source predictions available")
            continue

        matrix.index = matrix.index.map(
            lambda label: _deduplicate_modality_label(label, (mod1, mod2))
        )
        matrix.columns = matrix.columns.map(
            lambda label: _deduplicate_modality_label(label, (mod1, mod2))
        )
        matrix = matrix.loc[
            _modality_mask(matrix.index, mod1),
            _modality_mask(matrix.columns, mod2),
        ]
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
        threshold_suffix = ""
        if tool in PREDICTION_THRESHOLDS:
            operator, default_threshold = PREDICTION_THRESHOLDS[tool]
            effective_threshold = threshold if threshold is not None else default_threshold
            threshold_suffix = (
                f"_threshold_{operator}_{_format_threshold(effective_threshold)}"
            )
        output_path = tool_output_dir / f"{tool}{threshold_suffix}.csv"
        matrix.to_csv(output_path, index_label="sample_id")
        written_paths.append(output_path)
        print(f"Wrote {output_path} ({matrix.shape[0]} x {matrix.shape[1]})")

    return written_paths


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export standardized binary sample-match prediction matrices."
    )
    parser.add_argument("--data-path", required=True, help="Root directory containing tool outputs")
    parser.add_argument("--output-dir", required=True, help="Directory for standardized prediction CSVs")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--tools", nargs="+", required=True, choices=TOOLS)
    parser.add_argument(
        "--threshold",
        action="append",
        default=[],
        metavar="TOOL=VALUE",
        help="Override a tool threshold, e.g. --threshold Conpair=0.7 (repeatable)",
    )
    parser.add_argument("--read-depth", required=True, help="Read-depth selector for mod1")
    parser.add_argument("--read-depth-mod2", help="Optional read-depth selector for mod2")
    parser.add_argument(
        "--heatmap-dir",
        default=str(
            Path(__file__).resolve().parent.parent
            / "benchmarking"
            / "data"
            / "standardized_heatmaps"
        ),
        help="Root containing standardized heatmap CSVs (default: benchmarking/data/standardized_heatmaps)",
    )
    parser.add_argument("--mod1", required=True)
    parser.add_argument("--mod2", required=True)
    parser.add_argument("--pseudobulk", action="store_true")
    parser.add_argument("--ncells", type=int, help="Required when --pseudobulk is set")
    return parser


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()
    if args.pseudobulk and args.ncells is None:
        raise SystemExit("--ncells is required when --pseudobulk is set")
    try:
        threshold_overrides = _parse_threshold_overrides(args.threshold, args.tools)
    except ValueError as exc:
        parser.error(str(exc))

    written_paths = export_sample_predictions(
        args.data_path,
        args.output_dir,
        args.dataset,
        args.tools,
        args.read_depth,
        args.mod1,
        args.mod2,
        pseudobulk=args.pseudobulk,
        ncells=args.ncells,
        read_depth_mod2=args.read_depth_mod2,
        threshold_overrides=threshold_overrides,
        standardized_heatmaps_path=args.heatmap_dir,
    )
    return 0 if written_paths else 1


if __name__ == "__main__":
    raise SystemExit(main())