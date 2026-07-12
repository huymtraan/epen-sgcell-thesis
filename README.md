# Ependymoma single-cell downstream analysis

This repository provides a downstream workflow for integrating and annotating
public single-cell RNA-seq datasets of ependymoma. It covers Scanpy quality
control, normalization, batch correction, clustering, differential analysis,
CopyKAT, and evidence-guided cell annotation.

## Research scope

The analysis combines four public studies comprising 46 samples and 305,210
cells after quality control. It aims to describe major cellular compartments
and broad candidate malignant expression programs using canonical markers,
differential expression, Gene Ontology enrichment, and inferred copy-number
profiles.

One pattern of interest is a candidate **myeloid-mimicry malignant state**: a
cluster with immune/myeloid-associated expression, predominantly aneuploid
CopyKAT predictions, and motile-cilia-related enrichment. This is treated as a
working annotation rather than a new molecular subtype.

## Data processing

Public accessions and raw reads were retrieved from GEO/SRA with
[`nf-core/fetchngs`](https://nf-co.re/fetchngs) 1.12.0 and supplementary
SRA/Aspera download steps. FASTQ files were processed with a locally modified
[`nf-core/scrnaseq`](https://nf-co.re/scrnaseq) 4.1.0 workflow, using Cell
Ranger 9.0.1 and the 10x Genomics GRCh38 2024-A reference.

The resulting gene-by-cell matrices were converted to AnnData and concatenated
before entering this workflow. The expected input is described in
[`metadata/input_contract.md`](metadata/input_contract.md).

## Workflow

```text
raw-count AnnData
  → QC and Scrublet
  → shifted-log normalization
  → deviance feature selection
  → cell-cycle and mitochondrial regression
  → PCA and Harmony
  → neighbors, UMAP, and Leiden clustering
  → differential expression and GO enrichment
  → CopyKAT
  → annotation and figures
```

Core defaults are defined in [`config/analysis.yaml`](config/analysis.yaml):

- `min_genes = 100`, `min_cells = 3`;
- mitochondrial percentage ≤40%;
- per-sample MAD threshold 3;
- 2,000 deviance-selected genes and 50 PCs;
- Harmony by sample, 15 neighbors, Leiden resolution 0.4;
- regression of S score, G2/M score, and mitochondrial percentage.

## Usage

Install the package and run each stage in order:

```bash
pip install -e .

ependymoma-analysis qc          --input data/concatenated.h5ad --output outputs/01_qc.h5ad
ependymoma-analysis normalize   --input outputs/01_qc.h5ad --output outputs/02_pca.h5ad
ependymoma-analysis integrate   --input outputs/02_pca.h5ad --output outputs/03_harmony.h5ad
ependymoma-analysis cluster-deg --input outputs/03_harmony.h5ad --output outputs/04_clustered.h5ad
ependymoma-analysis copykat     --input outputs/04_clustered.h5ad --output outputs/05_copykat.h5ad
ependymoma-analysis annotate    --input outputs/05_copykat.h5ad --output outputs/06_annotated.h5ad
```

All commands use `config/analysis.yaml` by default. Deviance selection requires
R packages `SingleCellExperiment` and `scry`; CopyKAT requires `copykat`,
`stringr`, `rpy2`, and `infercnvpy`.

## Repository contents

- `src/ependymoma`: analysis modules and CLI;
- `config`: workflow parameters;
- `metadata`: study, sample, marker, and annotation tables;
- `notebooks`: annotation, CopyKAT visualization, and report figures;
- `reports`: workflow/figure provenance and reference run summaries;
- `tests`: metadata and configuration checks.

## Upstream workflow citations

- **nf-core/fetchngs 1.12.0** — [10.5281/zenodo.5070524](https://doi.org/10.5281/zenodo.5070524)
- **nf-core/scrnaseq 4.1.0**, locally modified to use Cell Ranger 9.0.1 — [10.5281/zenodo.3568187](https://doi.org/10.5281/zenodo.3568187)
- Ewels PA, et al. *The nf-core framework for community-curated bioinformatics pipelines.* Nature Biotechnology, 2020. [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)
- Di Tommaso P, et al. *Nextflow enables reproducible computational workflows.* Nature Biotechnology, 2017. [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
