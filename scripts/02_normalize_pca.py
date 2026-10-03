#!/usr/bin/env python3
"""
Stage 2: normalization, deviance-based HVG selection (scry via rpy2), cell-cycle regression, PCA.

Input: filtered AnnData from stage 1 (raw counts in X, var_names = gene symbols).
Output: AnnData with log1p layer, HVG mask, cell-cycle regressed layer, and PCA embedding.
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
        "--batch-key",
        default="sample",
        help="obs column used for batch/Harmony integration (default: sample)",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=2000,
        help="Number of HVGs to keep (deviance)",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of principal components",
    )
    parser.add_argument(
        "--s-genes",
        help="Path to newline-delimited S-phase genes (default: cell_cycle_s_genes.txt beside qc_utils.py)",
    )
    parser.add_argument(
        "--g2m-genes",
        help="Path to newline-delimited G2M-phase genes (default: cell_cycle_g2m_genes.txt beside qc_utils.py)",
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
    parser.add_argument(
        "--report",
        help="Path to write normalization/PCA report (default: <output>.norm_pca_report.tsv)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    import pandas as pd
    import scanpy as sc

    from qc_utils import (
        configure_threads,
        normalize_and_hvg_deviance,
        run_pca,
        score_and_regress_cell_cycle,
    )
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(message)s",
        level=getattr(logging, args.log_level.upper()),
    )
    configure_threads(args.threads)

    input_path = Path(args.input)
    output_path = Path(args.output)
    report_path = Path(args.report) if args.report else output_path.with_suffix(".norm_pca_report.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Loading AnnData from %s", input_path)
    adata = sc.read_h5ad(input_path)

    adata = normalize_and_hvg_deviance(
        adata,
        n_top_genes=args.n_top_genes,
        layer_name="log1p_norm",
        batch_key=args.batch_key,
    )
    adata = score_and_regress_cell_cycle(
        adata,
        layer_name="log1p_norm",
        s_genes_path=Path(args.s_genes) if args.s_genes else None,
        g2m_genes_path=Path(args.g2m_genes) if args.g2m_genes else None,
    )
    adata = run_pca(adata, layer_name="log1p_norm", n_pcs=args.n_pcs)

    report = pd.DataFrame(
        [
            {
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "batch_key": args.batch_key,
                "batch_aware_deviance": bool(args.batch_key),
                "n_top_genes": args.n_top_genes,
                "n_hvg_selected": int(adata.var["highly_variable"].sum()),
                "n_pcs": args.n_pcs,
                "n_samples": adata.obs[args.batch_key].nunique()
                if args.batch_key in adata.obs
                else None,
            }
        ]
    )
    report.to_csv(report_path, sep="\t", index=False)
    logging.info("Wrote norm/PCA report to %s", report_path)

    logging.info("Writing output to %s", output_path)
    adata.write_h5ad(output_path, compression="gzip")
    logging.info("Done")


if __name__ == "__main__":
    main()
