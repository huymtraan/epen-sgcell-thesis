"""Thesis-aligned Scanpy quality control and Scrublet handling."""

from __future__ import annotations

import logging

import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation

LOG = logging.getLogger(__name__)


def add_qc_metrics(adata, gene_symbols_key: str = "gene_symbols"):
    symbols = adata.var[gene_symbols_key].astype(str)
    adata.var["mt"] = symbols.str.startswith("MT-")
    adata.var["ribo"] = symbols.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = symbols.str.contains(r"^HB(?!P)", regex=True)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20],
        log1p=True,
        inplace=True,
    )
    return adata


def mad_outlier(values, nmads: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    median = np.nanmedian(values)
    mad = median_abs_deviation(values, nan_policy="omit", scale=1)
    if not np.isfinite(mad) or mad == 0:
        return np.zeros(values.shape, dtype=bool)
    return (values < median - nmads * mad) | (values > median + nmads * mad)


def flag_batch_outliers(
    adata,
    sample_key: str,
    metrics: list[str],
    nmads: float,
):
    """Flag a cell if any thesis-stated metric is outside sample median ± MAD."""
    missing = [metric for metric in metrics if metric not in adata.obs]
    if missing:
        raise KeyError(f"Missing QC metrics: {missing}")
    flags = np.zeros(adata.n_obs, dtype=bool)
    for _, positions in adata.obs.groupby(sample_key, observed=True).indices.items():
        sample_flags = np.zeros(len(positions), dtype=bool)
        for metric in metrics:
            sample_flags |= mad_outlier(adata.obs.iloc[positions][metric], nmads)
        flags[positions] = sample_flags
    adata.obs["outlier"] = flags
    return adata


def basic_filter(adata, config: dict, sample_key: str, gene_symbols_key: str):
    """Apply min genes/cells, mt ceiling, and per-sample MAD filters."""
    add_qc_metrics(adata, gene_symbols_key)
    sc.pp.filter_cells(adata, min_genes=int(config["min_genes"]))
    sc.pp.filter_genes(adata, min_cells=int(config["min_cells"]))
    adata = adata[adata.obs["pct_counts_mt"] <= float(config["max_pct_mito"])].copy()
    flag_batch_outliers(
        adata,
        sample_key=sample_key,
        metrics=list(config["mad_metrics"]),
        nmads=float(config["nmads"]),
    )
    return adata[~adata.obs["outlier"]].copy()


def remove_doublets_by_sample(
    adata,
    sample_key: str,
    requested_n_pcs: int = 30,
    random_state: int = 0,
):
    """Run Scrublet per sample with a safe PCA rank for small post-QC batches.

    The original batch-aware call failed when a sample retained fewer cells than
    the default 30 PCs. This implementation preserves per-sample operation while
    making the rank explicit and reporting samples that cannot be evaluated.
    """
    scores = np.full(adata.n_obs, np.nan, dtype=float)
    predicted = np.zeros(adata.n_obs, dtype=bool)
    evaluated = np.zeros(adata.n_obs, dtype=bool)
    failures: list[dict] = []

    for sample, positions in adata.obs.groupby(sample_key, observed=True).indices.items():
        subset = adata[positions].copy()
        n_pcs = min(requested_n_pcs, subset.n_obs - 1, subset.n_vars - 1)
        if n_pcs < 2:
            failures.append({"sample": str(sample), "reason": "fewer than 3 cells/features"})
            continue
        try:
            sc.pp.scrublet(
                subset,
                n_prin_comps=n_pcs,
                random_state=random_state,
                copy=False,
            )
        except Exception as exc:  # preserve sample and continue; report is mandatory
            failures.append({"sample": str(sample), "reason": str(exc)})
            continue
        scores[positions] = subset.obs["doublet_score"].to_numpy()
        predicted[positions] = subset.obs["predicted_doublet"].to_numpy(dtype=bool)
        evaluated[positions] = True

    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = predicted
    adata.obs["scrublet_evaluated"] = evaluated
    if failures:
        LOG.warning("Scrublet was not evaluated for %d samples", len(failures))
    return adata[~adata.obs["predicted_doublet"]].copy(), failures
