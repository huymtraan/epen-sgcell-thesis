"""Shifted-log normalization, deviance selection, regression, and PCA."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

LOG = logging.getLogger(__name__)


def shifted_log_normalize(adata, layer: str, target_sum: float | None = None):
    scaled = sc.pp.normalize_total(adata, target_sum=target_sum, inplace=False)
    adata.layers[layer] = sc.pp.log1p(scaled["X"], copy=True)
    return adata


def select_deviant_genes(adata, n_top_genes: int, batch_key: str | None = None):
    """Select genes with binomial deviance using R/scry via rpy2."""
    try:
        import anndata2ri
        import rpy2.rinterface_lib.callbacks as rcb
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
    except ImportError as exc:
        raise RuntimeError("Deviance selection requires rpy2, anndata2ri, and R/scry") from exc

    rcb.logger.setLevel(logging.ERROR)
    with ro.conversion.localconverter(
        ro.default_converter + pandas2ri.converter + anndata2ri.converter
    ):
        ro.globalenv["adata"] = adata
        ro.r("suppressPackageStartupMessages({library(SingleCellExperiment); library(scry)})")
        if batch_key:
            if batch_key not in adata.obs:
                raise KeyError(f"Missing batch key: {batch_key}")
            ro.globalenv["batch_key"] = batch_key
            ro.r('colData(adata)$batch_for_dev <- as.factor(colData(adata)[[batch_key]])')
            ro.r(
                'sce <- devianceFeatureSelection(adata, assay="X", '
                'batch=colData(adata)$batch_for_dev)'
            )
        else:
            ro.r('sce <- devianceFeatureSelection(adata, assay="X")')
        deviance = np.asarray(ro.r("rowData(sce)$binomial_deviance")).ravel()

    if not 0 < n_top_genes <= adata.n_vars:
        raise ValueError(f"n_top_genes must be within 1..{adata.n_vars}")
    selected = np.zeros(adata.n_vars, dtype=bool)
    selected[np.argsort(deviance)[-n_top_genes:]] = True
    adata.var["binomial_deviance"] = deviance
    adata.var["highly_deviant"] = selected
    adata.var["highly_variable"] = selected
    return adata


def _load_genes(path: str | Path) -> list[str]:
    return [
        line.strip().upper()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def score_cell_cycle(adata, layer: str, s_genes_path: str | Path, g2m_genes_path: str | Path):
    symbol_index = pd.Index(adata.var_names.astype(str).str.upper())
    lookup = dict(zip(symbol_index, adata.var_names))
    s_genes = [lookup[g] for g in _load_genes(s_genes_path) if g in lookup]
    g2m_genes = [lookup[g] for g in _load_genes(g2m_genes_path) if g in lookup]
    if not s_genes or not g2m_genes:
        raise ValueError("Cell-cycle gene lists do not overlap adata.var_names")
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, layer=layer)
    return adata


def pca_on_selected_genes(
    adata,
    layer: str,
    n_pcs: int,
    solver: str,
    regress_covariates: list[str],
):
    selected = adata.var["highly_variable"].to_numpy(dtype=bool)
    work = ad.AnnData(
        X=adata.layers[layer][:, selected].copy(),
        obs=adata.obs.copy(),
        var=adata.var.loc[selected].copy(),
    )
    if regress_covariates:
        missing = [key for key in regress_covariates if key not in work.obs]
        if missing:
            raise KeyError(f"Missing regression covariates: {missing}")
        sc.pp.regress_out(work, regress_covariates)
    sc.pp.scale(work)
    sc.pp.pca(work, n_comps=n_pcs, svd_solver=solver)
    adata.obsm["X_pca"] = work.obsm["X_pca"].copy()
    adata.uns["pca"] = work.uns["pca"].copy()
    loadings = np.zeros((adata.n_vars, n_pcs), dtype=work.varm["PCs"].dtype)
    loadings[selected] = work.varm["PCs"]
    adata.varm["PCs"] = loadings
    return adata
