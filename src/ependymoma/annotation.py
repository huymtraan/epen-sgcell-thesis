"""Deterministic final cluster annotation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


FINAL_COLORS = {
    "ependymal-like malignant": "#d6bcc0",
    "cycling cells": "#e07b91",
    "oligodendrocyte": "#bb7784",
    "neuronal-like malignant": "#4a6fe3",
    "translation-high malignant": "#8595e1",
    "stressed malignant": "#e6afb9",
    "neuron": "#f0b98d",
    "ependymal": "#9cded6",
    "astrocyte": "#d5eae7",
    "mural cells": "#bec1d4",
    "endothelial": "#b5bbe3",
    "myeloid-mimicry malignant": "#c6dec7",
    "myeloid": "#7d87b9",
    "T cell": "#8e063b",
}


def apply_annotation(adata, mapping_path: str | Path, cluster_key: str, label_key: str):
    mapping = pd.read_csv(mapping_path, sep="\t", dtype={cluster_key: str})
    if list(mapping.columns) != [cluster_key, "cell_type"]:
        raise ValueError(f"Annotation table must contain {cluster_key!r} and 'cell_type'")
    if mapping[cluster_key].duplicated().any():
        raise ValueError("Annotation table contains duplicate cluster IDs")
    lookup = mapping.set_index(cluster_key)["cell_type"]
    labels = adata.obs[cluster_key].astype(str).map(lookup)
    if labels.isna().any():
        missing = sorted(adata.obs.loc[labels.isna(), cluster_key].astype(str).unique())
        raise ValueError(f"Unmapped clusters: {missing}")
    categories = list(dict.fromkeys(mapping["cell_type"]))
    adata.obs[label_key] = pd.Categorical(labels, categories=categories)
    adata.uns[f"{label_key}_colors"] = [FINAL_COLORS[label] for label in categories]
    return adata
