#!/usr/bin/env python3
"""Run infercnvpy's CopyKAT wrapper independently for each sample."""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Clustered H5AD")
    parser.add_argument("--output-dir", required=True, help="Per-sample CopyKAT directory")
    parser.add_argument("--sample-key", default="sample")
    parser.add_argument("--project-key", default="project")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--min-cells", type=int, default=200)
    parser.add_argument("--threads", type=int, default=8)
    return parser.parse_args()


def main():
    args = parse_args()
    import infercnvpy as cnv
    import pandas as pd
    import scanpy as sc

    from qc_utils import configure_threads

    configure_threads(args.threads)
    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.input)
    if args.sample_key not in adata.obs:
        raise KeyError(f"missing sample key: {args.sample_key}")

    failed_rows, skipped_rows, done_rows = [], [], []
    previous_cwd = Path.cwd()
    for sample_id in adata.obs[args.sample_key].astype(str).unique():
        sample_dir = out_root / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        subset = adata[adata.obs[args.sample_key].astype(str) == sample_id].copy()
        n_cells = subset.n_obs
        project = (
            subset.obs[args.project_key].astype(str).iloc[0]
            if args.project_key in subset.obs else "NA"
        )
        if n_cells < args.min_cells:
            skipped_rows.append({"project": project, "sample": sample_id, "n_cells": n_cells, "reason": "too_few_cells"})
            del subset
            gc.collect()
            continue
        if args.counts_layer not in subset.layers:
            skipped_rows.append({"project": project, "sample": sample_id, "n_cells": n_cells, "reason": "missing_counts_layer"})
            del subset
            gc.collect()
            continue
        try:
            os.chdir(sample_dir)
            cnv.tl.copykat(
                subset,
                key_added="cnv_cpkat",
                layer=args.counts_layer,
                n_jobs=args.threads,
                inplace=False,
                s_name=f"{sample_id}_copykat",
            )
            done_rows.append({"project": project, "sample": sample_id, "n_cells": n_cells})
        except Exception as exc:  # preserve per-sample progress when CopyKAT fails
            failed_rows.append({"project": project, "sample": sample_id, "n_cells": n_cells, "error": str(exc).replace("\n", " ")})
        finally:
            os.chdir(previous_cwd)
            del subset
            gc.collect()

    pd.DataFrame(done_rows).to_csv(out_root / "copykat_done_samples.csv", index=False)
    pd.DataFrame(skipped_rows).to_csv(out_root / "copykat_skipped_samples.csv", index=False)
    pd.DataFrame(failed_rows).to_csv(out_root / "copykat_failed_samples.csv", index=False)


if __name__ == "__main__":
    main()
