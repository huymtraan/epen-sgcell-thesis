#!/usr/bin/env python3
"""Run Leiden resolutions on an existing Harmony neighbor graph."""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_RESOLUTIONS = [
    0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Harmony-integrated H5AD")
    parser.add_argument("--output", required=True, help="Clustered H5AD")
    parser.add_argument("--resolutions", nargs="+", type=float, default=DEFAULT_RESOLUTIONS)
    parser.add_argument("--threads", type=int, default=36)
    return parser.parse_args()


def main():
    args = parse_args()
    import scanpy as sc

    sc.settings.n_jobs = args.threads
    adata = sc.read_h5ad(args.input)
    for resolution in args.resolutions:
        key = f"leiden_{str(resolution).replace('.', '_')}"
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=2,
            key_added=key,
        )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")


if __name__ == "__main__":
    main()
