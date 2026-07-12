"""Per-sample CopyKAT with explicit result capture and coverage reporting."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)


def run_copykat_by_sample(adata, config: dict, output_dir: str | Path):
    """Run CopyKAT per sample and merge successful results into the master object.

    CopyKAT is optional and expensive. This function is never invoked during
    static reforge checks. It captures infercnvpy's returned matrix/prediction
    directly, avoiding the original script's ambiguous postprocessing failures.
    """
    import infercnvpy as cnv

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_key = config["sample_key"]
    project_key = config["project_key"]
    layer = config["counts_layer"]
    key_added = config["key_added"]
    prediction_key = config["prediction_key"]
    n_jobs = int(config["n_jobs"])

    merged_matrix = None
    merged_predictions = pd.Series(index=adata.obs_names, dtype="object")
    chrom_template = None
    summary: list[dict] = []

    for sample, positions in adata.obs.groupby(sample_key, observed=True).indices.items():
        subset = adata[positions].copy()
        project = str(subset.obs[project_key].iloc[0])
        if subset.n_obs < int(config["min_cells"]):
            summary.append(
                {"project": project, "sample": str(sample), "n_cells": subset.n_obs,
                 "status": "skipped", "reason": "below min_cells"}
            )
            continue
        if layer not in subset.layers:
            raise KeyError(f"Missing raw-count layer {layer!r}")

        # CopyKAT requires gene symbols as its input index. Only the sample copy
        # is changed; the master object's stable var_names remain untouched.
        subset.var_names = subset.var["gene_symbols"].astype(str)
        subset.var_names_make_unique()
        try:
            cnv.tl.copykat(
                subset,
                key_added=key_added,
                layer=layer,
                n_jobs=n_jobs,
                inplace=True,
                s_name=str(sample),
            )
        except Exception as exc:
            LOG.exception("CopyKAT failed for %s", sample)
            summary.append(
                {"project": project, "sample": str(sample), "n_cells": subset.n_obs,
                 "status": "failed", "reason": str(exc)}
            )
            continue

        matrix = np.asarray(subset.obsm[f"X_{key_added}"], dtype=np.float32)
        prediction = subset.obs[key_added].copy()
        sample_chrom = subset.uns[key_added]
        if chrom_template is None:
            chrom_template = sample_chrom
        elif sample_chrom != chrom_template:
            summary.append(
                {"project": project, "sample": str(sample), "n_cells": subset.n_obs,
                 "status": "failed", "reason": "CopyKAT chromosome bins differ from template"}
            )
            continue
        if merged_matrix is None:
            merged_matrix = np.full((adata.n_obs, matrix.shape[1]), np.nan, dtype=np.float32)
        if matrix.shape[1] != merged_matrix.shape[1]:
            summary.append(
                {"project": project, "sample": str(sample), "n_cells": subset.n_obs,
                 "status": "failed", "reason": "CopyKAT bin count differs from template"}
            )
            continue
        merged_matrix[positions] = matrix
        merged_predictions.loc[subset.obs_names] = prediction.reindex(subset.obs_names).to_numpy()
        sample_dir = output_dir / str(sample)
        sample_dir.mkdir(exist_ok=True)
        pd.DataFrame({prediction_key: prediction}).to_csv(
            sample_dir / "copykat_predictions.tsv", sep="\t"
        )
        summary.append(
            {"project": project, "sample": str(sample), "n_cells": subset.n_obs,
             "status": "done", "reason": ""}
        )

    if merged_matrix is None:
        raise RuntimeError("CopyKAT produced no successful sample result")
    adata.obsm[f"X_{key_added}"] = merged_matrix
    adata.uns[key_added] = chrom_template
    adata.obs[prediction_key] = merged_predictions.reindex(adata.obs_names)
    pd.DataFrame(summary).to_csv(output_dir / "copykat_summary.tsv", sep="\t", index=False)
    return adata, pd.DataFrame(summary)
