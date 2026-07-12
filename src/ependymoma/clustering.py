"""Final Leiden clustering reconstructed from thesis and notebook evidence."""

from __future__ import annotations

import scanpy as sc


def cluster(adata, config: dict):
    sc.tl.leiden(
        adata,
        resolution=float(config["resolution"]),
        key_added=config["key"],
        flavor=config["flavor"],
        n_iterations=int(config["n_iterations"]),
        random_state=int(config["random_state"]),
    )
    return adata
