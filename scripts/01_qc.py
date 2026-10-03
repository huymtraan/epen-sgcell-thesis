#!/usr/bin/env python3
"""
Stage 1: basic QC + MAD outlier masking + Scrublet doublet removal.

Output: filtered AnnData with var_names set to gene symbols and raw counts stored in .layers["counts"].
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input .h5ad file path")
    parser.add_argument("--output", required=True, help="Output .h5ad file path")
    parser.add_argument(
        "--metadata",
        help="Optional CSV with sample-level metadata to merge on the batch key",
    )
    parser.add_argument(
        "--batch-key",
        default="sample",
        help="obs column used for batch/Harmony integration (default: sample)",
    )
    parser.add_argument(
        "--trim-sample-prefix",
        action="store_true",
        help="Trim sample IDs before the first underscore (matches notebook helper)",
    )
    parser.add_argument(
        "--nmads",
        type=float,
        default=3.0,
        help="MAD threshold used for per-sample outlier detection",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=100,
        help="Minimum genes per cell filter",
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=3,
        help="Minimum cells per gene filter",
    )
    parser.add_argument(
        "--report",
        help="Path to write QC cell retention report (default: <output>.qc_report.tsv)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Limit threads for BLAS/OMP/Scanpy ops (e.g., 4)",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    import pandas as pd
    import scanpy as sc

    from qc_utils import (
        apply_basic_filters,
        configure_threads,
        keep_accession_prefix,
        merge_metadata,
        remove_doublets,
        stash_raw_counts_layer,
    )
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(message)s",
        level=getattr(logging, args.log_level.upper()),
    )
    configure_threads(args.threads)

    input_path = Path(args.input)
    output_path = Path(args.output)
    report_path = Path(args.report) if args.report else output_path.with_suffix(".qc_report.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Loading AnnData from %s", input_path)
    adata = sc.read_h5ad(input_path)

    report_rows = []

    def _n_samples(obj):
        return obj.obs[args.batch_key].nunique() if args.batch_key in obj.obs else None

    def log_step(step: str, before: int, after: int, obj):
        report_rows.append(
            {
                "step": step,
                "cells_retained": after,
                "removed_this_step": before - after,
                "removed_cumulative": start_cells - after,
                "samples_retained": _n_samples(obj),
            }
        )
        logging.info(
            "[%s] retained=%d removed_this_step=%d removed_cumulative=%d samples=%s",
            step,
            after,
            before - after,
            start_cells - after,
            _n_samples(obj),
        )

    start_cells = adata.n_obs
    log_step("load", start_cells, start_cells, adata)

    stash_raw_counts_layer(adata, layer_name="counts")

    if args.trim_sample_prefix:
        keep_accession_prefix(adata, col=args.batch_key)

    merge_metadata(adata, Path(args.metadata) if args.metadata else None, args.batch_key)

    adata = apply_basic_filters(
        adata,
        sample_key=args.batch_key,
        nmads=args.nmads,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
    )
    log_step("basic_qc_mad", start_cells, adata.n_obs, adata)
    before_doublet = adata.n_obs
    adata = remove_doublets(adata, batch_key=args.batch_key)
    log_step("remove_doublets", before_doublet, adata.n_obs, adata)

    pd.DataFrame(report_rows).to_csv(report_path, sep="\t", index=False)
    logging.info("Wrote QC report to %s", report_path)

    logging.info("Writing output to %s", output_path)
    adata.write_h5ad(output_path, compression="gzip")
    logging.info("Done")


if __name__ == "__main__":
    main()
