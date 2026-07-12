"""Command-line stages for the reconstructed downstream workflow."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import scanpy as sc

from .annotation import apply_annotation
from .clustering import cluster
from .config import load_config
from .copykat import run_copykat_by_sample
from .differential_expression import rank_cluster_genes
from .enrichment import enrich_clusters
from .integration import integrate
from .io import validate_input, write_tsv
from .normalization import (
    pca_on_selected_genes,
    score_cell_cycle,
    select_deviant_genes,
    shifted_log_normalize,
)
from .qc import basic_filter, remove_doublets_by_sample


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["qc", "normalize", "integrate", "cluster-deg", "copykat", "annotate"])
    parser.add_argument("--config", default="config/analysis.yaml")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report-dir", default="reports/run")
    parser.add_argument("--s-genes", default="metadata/cell_cycle_s_genes.txt")
    parser.add_argument("--g2m-genes", default="metadata/cell_cycle_g2m_genes.txt")
    parser.add_argument("--annotation-map", default="metadata/annotation_mapping.tsv")
    parser.add_argument("--log-level", default="INFO")
    return parser


def _write(adata, output: str) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")


def main() -> None:
    args = _parser().parse_args()
    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s %(message)s")
    cfg = load_config(args.config)
    adata = sc.read_h5ad(args.input)
    reports = Path(args.report_dir)

    if args.stage == "qc":
        icfg, qcfg = cfg["input"], cfg["qc"]
        validate_input(adata, icfg["sample_key"], icfg["project_key"], icfg["gene_symbols_key"])
        if icfg["counts_layer"] not in adata.layers:
            adata.layers[icfg["counts_layer"]] = adata.X.copy()
        start = adata.n_obs
        adata = basic_filter(adata, qcfg, icfg["sample_key"], icfg["gene_symbols_key"])
        after_basic = adata.n_obs
        adata, failures = remove_doublets_by_sample(
            adata,
            sample_key=icfg["sample_key"],
            requested_n_pcs=int(qcfg["scrublet_n_prin_comps"]),
            random_state=int(qcfg["scrublet_random_state"]),
        )
        write_tsv([
            {"step": "load", "cells_retained": start, "removed": 0},
            {"step": "basic_qc_mad", "cells_retained": after_basic, "removed": start - after_basic},
            {"step": "remove_doublets", "cells_retained": adata.n_obs, "removed": after_basic - adata.n_obs},
        ], reports / "qc_report.tsv")
        write_tsv(failures, reports / "scrublet_failures.tsv")

    elif args.stage == "normalize":
        ncfg, dcfg = cfg["normalization"], cfg["dimensionality"]
        adata = shifted_log_normalize(adata, ncfg["log_layer"], ncfg["target_sum"])
        adata = select_deviant_genes(
            adata,
            int(ncfg["n_top_genes"]),
            dcfg["batch_key"] if ncfg["batch_aware_deviance"] else None,
        )
        covariates = list(ncfg["regress_covariates"])
        if any(key in covariates for key in ("S_score", "G2M_score")):
            adata = score_cell_cycle(adata, ncfg["log_layer"], args.s_genes, args.g2m_genes)
        adata = pca_on_selected_genes(
            adata,
            ncfg["log_layer"],
            int(dcfg["n_pcs"]),
            dcfg["pca_solver"],
            covariates,
        )

    elif args.stage == "integrate":
        adata = integrate(adata, cfg["dimensionality"])

    elif args.stage == "cluster-deg":
        adata = cluster(adata, cfg["clustering"])
        adata = rank_cluster_genes(
            adata,
            cfg["clustering"]["key"],
            cfg["differential_expression"],
        )
        enrich_clusters(
            adata,
            cfg["clustering"]["key"],
            cfg["differential_expression"]["key"],
            cfg["enrichment"],
            reports / "go_enrichment",
        )

    elif args.stage == "copykat":
        adata, _ = run_copykat_by_sample(adata, cfg["copykat"], reports / "copykat")

    elif args.stage == "annotate":
        acfg = cfg["annotation"]
        adata = apply_annotation(
            adata,
            args.annotation_map,
            acfg["cluster_key"],
            acfg["label_key"],
        )

    _write(adata, args.output)


if __name__ == "__main__":
    main()
