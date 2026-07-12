"""Small figure functions corresponding to thesis Figures 2–5."""

from __future__ import annotations

import numpy as np
import scanpy as sc


def leiden_umap(adata, cluster_key: str = "leiden_0_4"):
    return sc.pl.umap(adata, color=cluster_key, legend_loc="on data", show=False)


def annotation_umap(adata, label_key: str):
    return sc.pl.umap(
        adata,
        color=label_key,
        title="Cell Type",
        legend_loc="right margin",
        size=5,
        alpha=0.2,
        frameon=False,
        show=False,
    )


def copykat_umap(adata, prediction_key: str = "copykat_pred"):
    return sc.pl.umap(adata, color=prediction_key, title="CopyKAT prediction", show=False)


def deterministic_copykat_subset(
    adata,
    sample_key: str = "sample",
    prediction_key: str = "copykat_pred",
    fraction: float = 0.10,
    random_state: int = 0,
):
    """Recreate the thesis heatmap sampling while retaining both classes when possible."""
    rng = np.random.default_rng(random_state)
    keep = []
    valid = adata.obs[prediction_key].isin(["aneuploid", "diploid"])
    for _, group in adata.obs.loc[valid].groupby(sample_key, observed=True):
        target = max(1, int(np.ceil(len(group) * fraction)))
        chosen = []
        for label in ("aneuploid", "diploid"):
            members = group.index[group[prediction_key] == label]
            if len(members):
                chosen.append(rng.choice(members))
        target = max(target, len(chosen))
        remaining = group.drop(index=chosen)
        if target > len(chosen) and len(remaining):
            chosen.extend(rng.choice(remaining.index, size=min(target - len(chosen), len(remaining)), replace=False))
        keep.extend(chosen)
    return adata[keep].copy()
