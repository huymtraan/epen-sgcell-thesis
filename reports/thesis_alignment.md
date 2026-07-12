# Thesis alignment

| Thesis statement | Staging implementation | Evidence/status |
|---|---|---|
| Four studies, 46 samples | `metadata/studies.tsv`, `samples.tsv` | Confirmed in final H5AD metadata |
| min_cells 3 | `config/analysis.yaml` | Exact |
| min_genes 100 | `config/analysis.yaml` | Exact |
| mitochondrial percentage ≤40% | QC module | Exact |
| MAD 3 over total counts, genes, top-20 fraction | QC module | Exact; mt is not included in MAD |
| Scrublet | per-sample QC module | Method exact; dynamic PC rank fixes observed small-batch crash |
| shifted-log normalization | normalization module | Exact |
| 2,000 deviance-selected genes | normalization module | Exact; R/scry dependency |
| PCA 50 | config/normalization module | Exact |
| Harmony | integration module by sample | Exact |
| neighbor graph in Harmony space | integration module, 15 neighbors | 15 from successful report |
| Leiden 0.4 | clustering module | Exact |
| Scanpy Wilcoxon DEG | DEG module | Supported by final H5AD and scripts |
| top 1,000 DEGs, GO:BP 2021 via gseapy | enrichment module | Exact thesis description |
| CopyKAT default parameters, internal baseline | CopyKAT module | Defaults retained except operational jobs/sample split |
| marker + DEG + GO + CNV annotation | deterministic mapping + notebooks | Marker panels and cluster mapping included |

Regression of `S_score`, `G2M_score`, and `pct_counts_mt` is explicit in the
analysis configuration.
