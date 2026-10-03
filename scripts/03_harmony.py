#!/usr/bin/env python3
"""
Stage 3: Harmony integration followed by neighbor graph construction.

Input: AnnData with PCA results (e.g., output of stage 2).
Output: AnnData with Harmony-corrected PCs and neighbors.
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
        help="obs column used for Harmony integration (default: sample)",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
        help="Number of neighbors for graph construction after Harmony",
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
        help="Path to write Harmony/neighbors report (default: <output>.harmony_report.tsv)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    import pandas as pd
    import scanpy as sc

    from qc_utils import configure_threads
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(message)s",
        level=getattr(logging, args.log_level.upper()),
    )
    configure_threads(args.threads)

    input_path = Path(args.input)
    output_path = Path(args.output)
    report_path = Path(args.report) if args.report else output_path.with_suffix(".harmony_report.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Loading AnnData from %s", input_path)
    adata = sc.read_h5ad(input_path)

    logging.info("Running Harmony with batch_key=%s", args.batch_key)
    sc.external.pp.harmony_integrate(adata, key=args.batch_key, basis="X_pca")

    logging.info("Computing neighbors on Harmony embedding (n_neighbors=%d)", args.n_neighbors)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=args.n_neighbors)
    sc.tl.umap(adata)

    # report critical settings
    n_pcs_harmony = adata.obsm["X_pca_harmony"].shape[1] if "X_pca_harmony" in adata.obsm else None
    report = pd.DataFrame(
        [
            {
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "batch_key": args.batch_key,
                "n_pcs_harmony": n_pcs_harmony,
                "n_neighbors": args.n_neighbors,
                "n_samples": adata.obs[args.batch_key].nunique()
                if args.batch_key in adata.obs
                else None,
            }
        ]
    )
    report.to_csv(report_path, sep="\t", index=False)
    logging.info("Wrote Harmony/neighbors report to %s", report_path)

    logging.info("Writing output to %s", output_path)
    adata.write_h5ad(output_path, compression="gzip")
    logging.info("Done")


if __name__ == "__main__":
    main()
