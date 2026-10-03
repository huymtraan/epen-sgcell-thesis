# Integrated single-cell profiling of public ependymoma datasets

**Original undergraduate thesis title**

> **Characterizing the Intra-Tumoral Heterogeneity of Ependymomas Using
> Single-Cell RNA-seq on Multi-center Datasets**

## Scope

The thesis was designed to characterize intra-tumoral heterogeneity in
ependymoma across multiple public datasets. The completed analysis integrated
four single-cell RNA-seq studies and included quality control, normalization,
Harmony integration, clustering, differential expression, Gene Ontology
enrichment, CopyKAT-based CNV inference, and cell annotation.

The integrated dataset supports descriptive profiling of cellular compartments
and transcriptional programs. It does not establish that the reported cell
states are reproducible across cohorts or, by itself, support the broader claim
of characterizing intra-tumoral heterogeneity across studies implied by the
original thesis title. Study-specific replication, leave-one-study-out analysis,
formal cross-study comparison, and systematic re-evaluation of the source
studies were not performed because of the time and methodological constraints
of the thesis project.

## Dataset

The final integrated object contains 305,210 cells from 46 samples across four
public studies.

| Study | BioProject | Samples |
|---|---|---:|
| GSE125969 | PRJNA518131 | 26 |
| GSE163686 | PRJNA687170 | 15 |
| GSE189939 | PRJNA785181 | 4 |
| GSE245674 | PRJNA1029546 | 1 |

Only the ependymoma sample SAMN37873687 was included from GSE245674.

## Analysis workflow

```text
Raw-count AnnData
  → quality control and Scrublet
  → shifted-log normalization
  → deviance-based feature selection
  → cell-cycle and mitochondrial regression
  → PCA and Harmony integration
  → neighbors, UMAP, and Leiden clustering
  → differential expression and GO enrichment
  → CopyKAT CNV inference
  → cell annotation and figures
```

Final analysis parameters: per-sample MAD threshold 3, 2,000
deviance-selected genes, 50 principal components, 15 neighbors, and Leiden
resolution 0.4.

## Main outputs

- Integrated UMAP after quality control and Harmony correction.
- Annotation of immune, vascular, stromal, and resident neural populations.
- Candidate ependymal-like, neuronal-like, cycling, stress-associated, and
  translation-high malignant programs.
- Candidate myeloid-mimicry malignant population with myeloid-associated
  expression and predominantly aneuploid CopyKAT predictions.
- CopyKAT prediction on UMAP and a chromosome-level CNV heatmap.

Cell labels are working annotations based on canonical markers, cluster-level
differential expression, GO enrichment, and inferred CNV profiles.

## Repository contents

| Path | Content |
|---|---|
| `notebooks/01_preprocessing.ipynb` | Data concatenation, QC, normalization, feature selection, PCA, and Harmony |
| `notebooks/02_clustering.ipynb` | Neighbor graph, UMAP, Leiden clustering, and differential expression |
| `notebooks/03_copykat_cnv.ipynb` | Per-sample CopyKAT results, CNV matrix assembly, and chromosome heatmap |
| `notebooks/04_annotation_results.ipynb` | Marker-, DEG-, GO-, and CNV-based annotation and final figures |
| `scripts/` | Long-running Python and Slurm jobs |
| `resources/` | Reference reports and small analysis outputs |

## Data paths and checkpoints

| Variable | Purpose | Default |
|---|---|---|
| `EPEN_SOURCE_DATA` | Existing H5AD checkpoints | Parent workspace `analysis/anndata` |
| `EPEN_OUTPUT_DIR` | Newly written checkpoints | Repository `outputs/` |

```text
00_concatenated.h5ad
01_qc.h5ad
02_norm_pca.h5ad
03_harmony.h5ad
04_clustered.h5ad
05_copykat_merged.h5ad
06_annotated.h5ad
```

## Batch execution

```bash
export EPEN_INPUT_H5AD=/path/to/concatenated.h5ad
export EPEN_OUTPUT_DIR=/path/to/outputs
bash scripts/run_processing.sh
```

Slurm entry points:

```text
scripts/run_processing.sbatch
scripts/run_copykat.sbatch
```

## Dependencies

The preprocessing, clustering, and annotation notebooks use Scanpy/AnnData,
Scrublet, Harmony, igraph/Leiden, pyclustree, rpy2/anndata2ri, R
SingleCellExperiment/scry, and gseapy.

CopyKAT analysis uses a separate environment with infercnvpy and CopyKAT. The
notebook metadata records the original kernels `sgcell-scanpy` and
`infercnvpy`.

## Data availability

Large expression matrices and intermediate H5AD files are not included. Public
source data are identified by GEO and BioProject accession in the dataset table.
