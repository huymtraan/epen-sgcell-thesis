#!/usr/bin/env python3
"""
Shared helpers for the qc_to_harmony pipeline stages.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

LOG = logging.getLogger(__name__)

DEFAULT_S_GENES_PATH = Path(__file__).with_name("cell_cycle_s_genes.txt")
DEFAULT_G2M_GENES_PATH = Path(__file__).with_name("cell_cycle_g2m_genes.txt")


def configure_threads(threads: Optional[int]):
    """Limit parallelism for BLAS/numexpr/OMP-heavy ops and Scanpy."""
    if threads is None:
        return
    if threads < 1:
        raise ValueError("threads must be >=1")
    for var in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        os.environ[var] = str(threads)
    os.environ["MKL_DYNAMIC"] = "FALSE"
    os.environ["OMP_WAIT_POLICY"] = "PASSIVE"
    sc.settings.n_jobs = threads
    LOG.info("Set thread limits to %d", threads)


def keep_accession_prefix(adata, col: str = "sample", out_col: Optional[str] = None):
    """
    Keep only the token before the first underscore in `adata.obs[col]`.
    """
    if col not in adata.obs:
        raise KeyError(f"Column '{col}' not found in adata.obs")

    series = adata.obs[col].astype(str)
    accession = series.str.extract(r"^([^_]+)")[0]
    target = col if out_col is None else out_col
    adata.obs[target] = pd.Categorical(accession.astype(str))
    return adata


def merge_metadata(adata, metadata_csv: Optional[Path], sample_key: str):
    """Optionally merge sample-level metadata CSV on the sample_key."""
    if metadata_csv is None:
        return adata
    meta = pd.read_csv(metadata_csv)
    if sample_key not in meta.columns:
        raise KeyError(f"'{sample_key}' column not found in {metadata_csv}")
    adata.obs = adata.obs.merge(
        meta.set_index(sample_key), left_on=sample_key, right_index=True, how="left"
    )
    return adata


# def set_var_names_from_gene_symbols(adata):
#     """Replace var_names with gene_symbols and make them unique."""
#     if "gene_symbols" not in adata.var:
#         raise KeyError("gene_symbols column required to set var_names")
#     adata.var_names = adata.var.gene_symbols.astype(str).to_numpy()
#     adata.var_names_make_unique()
#     return adata


def stash_raw_counts_layer(adata, layer_name: str = "counts"):
    """Stash raw counts from X into a layer if not already present."""
    if layer_name not in adata.layers:
        LOG.info("Stashing raw counts in .layers['%s']", layer_name)
        adata.layers[layer_name] = adata.X.copy()
    return adata


def add_qc_annotations(adata):
    """Add mt/ribo/hb flags and standard QC metrics."""
    adata.var["mt"] = adata.var.gene_symbols.str.startswith("MT-")
    adata.var["ribo"] = adata.var.gene_symbols.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var.gene_symbols.str.contains(r"^HB(?!P)")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    return adata


def _is_outlier(arr: Iterable[float], nmads: float) -> np.ndarray:
    """Return boolean mask of points outside median ± nmads * MAD."""
    arr = np.asarray(arr)
    med = np.median(arr)
    mad = median_abs_deviation(arr, scale=1)
    lower = med - nmads * mad
    upper = med + nmads * mad
    return (arr < lower) | (arr > upper)


def flag_outliers_by_batch(adata, sample_key: str = "sample", nmads: float = 3.0):
    """
    Recreate the per-sample outlier detection from preprocessing.ipynb.
    Adds a boolean `outlier` column to adata.obs.
    """
    if sample_key not in adata.obs:
        raise KeyError(f"Expected '{sample_key}' in adata.obs for batch QC")

    adata.obs["outlier"] = np.zeros(adata.n_obs, dtype=bool)
    # metrics = [
    #     "log1p_total_counts",
    #     "log1p_n_genes_by_counts",
    #     "pct_counts_in_top_20_genes",
    #     "pct_counts_mt",
    # ]
    for _, idx in adata.obs.groupby(sample_key).groups.items():
        subset = adata[idx]
        sample_outliers = (
            _is_outlier(subset.obs["log1p_total_counts"], nmads)
            | _is_outlier(subset.obs["log1p_n_genes_by_counts"], nmads)
            | _is_outlier(subset.obs["pct_counts_in_top_20_genes"], nmads)
        )
        adata.obs.loc[idx, "outlier"] = sample_outliers
    return adata


def apply_basic_filters(
    adata,
    sample_key: str = "sample",
    nmads: float = 3.0,
    min_genes: int = 100,
    min_cells: int = 3,
):
    """QC metrics + min gene/cell filters, then MAD-based outlier masking."""
    add_qc_annotations(adata)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if "pct_counts_mt" in adata.obs:
        before_mt = adata.n_obs
        adata = adata[adata.obs["pct_counts_mt"] <= 40].copy()
        LOG.info(
            "Removed %d cells with pct_counts_mt > 40 before MAD outlier calc",
            before_mt - adata.n_obs,
        )
    flag_outliers_by_batch(adata, sample_key=sample_key, nmads=nmads)
    adata = adata[~adata.obs["outlier"]].copy()
    return adata


def remove_doublets(adata, batch_key: str):
    """Scrublet batch-aware doublet detection and filtering."""
    LOG.info("Running Scrublet with batch_key=%s", batch_key)
    sc.pp.scrublet(adata, batch_key=batch_key)
    if "predicted_doublet" in adata.obs:
        before = adata.n_obs
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        LOG.info("Removed %d predicted doublets", before - adata.n_obs)
    else:
        LOG.warning("Scrublet did not add 'predicted_doublet'; skipping removal")
    return adata


def load_gene_list(path: Path) -> list[str]:
    """Load a simple newline-delimited gene list (ignoring blank/comment lines)."""
    genes = []
    with open(path) as handle:
        for line in handle:
            gene = line.strip()
            if not gene or gene.startswith("#"):
                continue
            genes.append(gene.upper())
    return genes


def deviance_feature_selection(
    adata,
    n_top_genes: int = 2000,
    batch_key: Optional[str] = None,
):
    """
    Run scry::devianceFeatureSelection via rpy2 to mark highly deviant genes.

    Parameters
    ----------
    adata : AnnData
        AnnData object with raw counts in .X.
    n_top_genes : int, default=2000
        Number of top deviant genes to keep globally.
    batch_key : str or None, default=None
        Column in adata.obs to use as batch factor.
        If None, deviance is computed ignoring batch (original behavior).
    """
    try:
        import anndata2ri
        import rpy2.rinterface_lib.callbacks as rcb
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
    except ImportError as exc:
        raise RuntimeError(
            "rpy2 and anndata2ri are required for deviance feature selection; "
            "install them alongside the R package 'scry'."
        ) from exc

    if batch_key is not None and batch_key not in adata.obs.columns:
        raise ValueError(
            f"batch_key='{batch_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    rcb.logger.setLevel(logging.ERROR)

    with ro.conversion.localconverter(
        ro.default_converter + pandas2ri.converter + anndata2ri.converter
    ):
        ro.globalenv["adata"] = adata
        ro.r(
            "suppressPackageStartupMessages({"
            "  library(SingleCellExperiment);"
            "  library(scry)"
            "})"
        )

        if batch_key is None:
            ro.r('sce <- devianceFeatureSelection(adata, assay = "X")')
        else:
            ro.globalenv["batch_key"] = batch_key
            ro.r(
                'colData(adata)$batch_for_dev <- as.factor(colData(adata)[[batch_key]])'
            )
            ro.r(
                'sce <- devianceFeatureSelection('
                '  adata, assay = "X", '
                '  batch = colData(adata)$batch_for_dev'
                ')'
            )

        dev = np.asarray(ro.r("rowData(sce)$binomial_deviance")).ravel()

    if dev.shape[0] != adata.n_vars:
        raise ValueError(
            f"binomial_deviance length ({dev.shape[0]}) does not match n_vars ({adata.n_vars})"
        )
    if n_top_genes > adata.n_vars:
        raise ValueError(
            f"n_top_genes={n_top_genes} is larger than number of genes ({adata.n_vars})"
        )

    idx = np.argsort(dev)[-n_top_genes:]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[idx] = True
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = dev
    return adata


def normalize_and_hvg_deviance(
    adata,
    n_top_genes: int = 2000,
    layer_name: str = "log1p_norm",
    batch_key: Optional[str] = None,
):
    """Normalize counts, log-transform, and select HVGs with R/scry deviance."""
    LOG.info("Normalizing counts and selecting deviant HVGs via scry (rpy2)")
    scaled = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers[layer_name] = sc.pp.log1p(scaled["X"], copy=True)

    deviance_feature_selection(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    sc.pp.highly_variable_genes(adata, layer=layer_name)
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    return adata


def score_and_regress_cell_cycle(
    adata,
    layer_name: str = "log1p_norm",
    s_genes_path: Optional[Path] = None,
    g2m_genes_path: Optional[Path] = None,
):
    """
    Score S/G2M cell cycle on the given layer.

    Regression is intentionally deferred to run_pca to avoid densifying the
    full matrix; here we only add S_score and G2M_score to adata.obs.
    """
    s_path = Path(s_genes_path) if s_genes_path else DEFAULT_S_GENES_PATH
    g2m_path = Path(g2m_genes_path) if g2m_genes_path else DEFAULT_G2M_GENES_PATH
    try:
        s_genes_raw = load_gene_list(s_path)
        g2m_genes_raw = load_gene_list(g2m_path)
    except FileNotFoundError:
        LOG.warning(
            "Cell cycle gene lists not found (%s, %s); using Scanpy CYCLING_GENES",
            s_path,
            g2m_path,
        )
        from scanpy.tools._score_genes_cell_cycle import CYCLING_GENES

        s_genes_raw = [g.upper() for g in CYCLING_GENES["S"]]
        g2m_genes_raw = [g.upper() for g in CYCLING_GENES["G2M"]]

    var_upper = pd.Index(adata.var_names).str.upper()
    upper_to_var = dict(zip(var_upper, adata.var_names))
    s_genes = [upper_to_var[g] for g in s_genes_raw if g in upper_to_var]
    g2m_genes = [upper_to_var[g] for g in g2m_genes_raw if g in upper_to_var]

    if len(s_genes) == 0 or len(g2m_genes) == 0:
        LOG.warning("No cell cycle genes found in var_names; skipping regression")
        return adata

    LOG.info(
        "Scoring cell cycle (%d S genes, %d G2M genes) on layer '%s'",
        len(s_genes),
        len(g2m_genes),
        layer_name,
    )
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes,
        g2m_genes=g2m_genes,
        layer=layer_name,
    )
    return adata


def run_pca(adata, layer_name: str = "log1p_norm", n_pcs: int = 50):
    """
    Run PCA using deviance-based HVGs while keeping the main object sparse.

    Workflow:
    - build a temporary HVG-only AnnData with X = adata.layers[layer_name][:, HVGs]
    - optionally regress out S_score/G2M_score on that smaller matrix
    - scale and run PCA
    - copy the PCA embedding/loadings back, keeping counts sparse in the main adata
    """
    if layer_name not in adata.layers:
        raise ValueError(
            f"Layer '{layer_name}' not found in adata.layers. "
            "Make sure normalize_and_hvg_deviance has been run."
        )
    if "highly_variable" not in adata.var:
        raise ValueError(
            "No 'highly_variable' column in adata.var. "
            "Run normalize_and_hvg_deviance before run_pca."
        )

    hvg_mask = adata.var["highly_variable"].values
    n_hvg = int(hvg_mask.sum())
    if n_hvg == 0:
        raise ValueError("No highly variable genes selected ('highly_variable' sum is 0).")

    LOG.info("Building temporary HVG-only AnnData for PCA (n_hvg=%d)", n_hvg)
    adata_hvg = adata[:, hvg_mask].copy()
    adata_hvg.X = adata.layers[layer_name][:, hvg_mask].copy()
    adata_hvg.var["highly_variable"] = True

    if "S_score" in adata_hvg.obs and "G2M_score" in adata_hvg.obs:
        LOG.info("Regressing out S_score and G2M_score on HVGs only")
        sc.pp.regress_out(adata_hvg, ["S_score", "G2M_score","pct_counts_mt"])
    else:
        LOG.warning(
            "S_score / G2M_score not found in adata.obs; skipping regression in PCA step."
        )

    LOG.info("Scaling HVGs for PCA")
    sc.pp.scale(adata_hvg)

    LOG.info("Running PCA (n_pcs=%d) on HVG-only object", n_pcs)
    sc.pp.pca(adata_hvg, svd_solver="arpack", n_comps=n_pcs)

    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
    adata.uns["pca"] = adata_hvg.uns["pca"].copy()

    PCs_hvg = adata_hvg.varm["PCs"]
    PCs_full = np.zeros((adata.n_vars, PCs_hvg.shape[1]), dtype=PCs_hvg.dtype)
    PCs_full[hvg_mask, :] = PCs_hvg
    adata.varm["PCs"] = PCs_full

    LOG.info(
        "PCA complete; results stored in adata.obsm['X_pca'], adata.varm['PCs'], adata.uns['pca']"
    )
    return adata
