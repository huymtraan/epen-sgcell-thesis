"""Harmony integration, neighbors, and UMAP."""

from __future__ import annotations

import scanpy as sc


def integrate(adata, config: dict):
    batch_key = config["batch_key"]
    if batch_key not in adata.obs:
        raise KeyError(f"Missing Harmony batch key: {batch_key}")
    sc.external.pp.harmony_integrate(
        adata,
        key=batch_key,
        basis=config["harmony_basis"],
        adjusted_basis=config["harmony_adjusted_basis"],
    )
    sc.pp.neighbors(
        adata,
        use_rep=config["harmony_adjusted_basis"],
        n_neighbors=int(config["n_neighbors"]),
        random_state=int(config["random_state"]),
    )
    sc.tl.umap(adata, random_state=int(config["random_state"]))
    return adata
