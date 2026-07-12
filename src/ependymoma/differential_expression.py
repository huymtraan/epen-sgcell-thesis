"""Cluster differential expression matching the final Leiden key."""

from __future__ import annotations

import scanpy as sc


def rank_cluster_genes(adata, cluster_key: str, config: dict):
    if cluster_key not in adata.obs:
        raise KeyError(f"Missing cluster key: {cluster_key}")
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        reference=config["reference"],
        n_genes=None,
        method=config["method"],
        layer=config["layer"],
        pts=bool(config["pts"]),
        key_added=config["key"],
    )
    return adata
