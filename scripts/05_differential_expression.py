#!/usr/bin/env python3
"""Compute thesis cluster differential expression at Leiden resolution 0.4."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Clustered H5AD")
    parser.add_argument("--output", required=True, help="H5AD with DEG results")
    parser.add_argument("--cluster-key", default="leiden_0_4")
    parser.add_argument("--layer", default="log1p_norm")
    return parser.parse_args()


def main():
    args = parse_args()
    import scanpy as sc

    adata = sc.read_h5ad(args.input)
    result_key = f"rank_genes_groups_{args.cluster_key}"
    sc.tl.rank_genes_groups(
        adata,
        groupby=args.cluster_key,
        reference="rest",
        n_genes=None,
        method="wilcoxon",
        layer=args.layer,
        pts=True,
        key_added=result_key,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")


if __name__ == "__main__":
    main()
