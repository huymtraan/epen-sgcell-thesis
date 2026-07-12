"""GO Biological Process enrichment for cluster DEGs."""

from __future__ import annotations

from pathlib import Path

import gseapy as gp
import pandas as pd
import scanpy as sc


def enrich_clusters(adata, cluster_key: str, deg_key: str, config: dict, output_dir: str | Path):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, pd.DataFrame] = {}
    for cluster in adata.obs[cluster_key].cat.categories:
        genes = sc.get.rank_genes_groups_df(
            adata,
            group=str(cluster),
            key=deg_key,
            gene_symbols="gene_symbols",
        ).head(int(config["top_n_genes"]))["gene_symbols"].dropna().astype(str).tolist()
        enriched = gp.enrichr(
            gene_list=genes,
            gene_sets=config["gene_set"],
            organism="Human",
            outdir=None,
        ).results
        enriched = enriched[
            enriched["Adjusted P-value"] <= float(config["adjusted_pvalue_max"])
        ].copy()
        enriched.to_csv(output_dir / f"cluster_{cluster}_go_bp.tsv", sep="\t", index=False)
        results[str(cluster)] = enriched
    return results
