# AnnData input contract

The workflow starts from one concatenated H5AD containing unnormalized UMI
counts for the 46 samples listed in `samples.tsv`.

Required:

- unique `obs_names` representing cells;
- unique `var_names` representing stable gene identifiers;
- `adata.X` containing non-negative integer-like raw counts;
- `adata.obs["sample"]` containing BioSample accessions;
- `adata.obs["project"]` containing one of the four included BioProjects;
- `adata.var["gene_symbols"]` containing human gene symbols.

QC copies `X` to `layers["counts"]` before filtering. CopyKAT uses the counts
layer and requires gene symbols as `var_names`; the CopyKAT stage therefore
creates a per-sample temporary view with gene-symbol indices without mutating
the master object's stable variable index.

The input file is not included in Git. Its creation through download,
Cell Ranger/STAR, nf-core/scrnaseq, matrix conversion, and concatenation is
outside the executable scope of this repository.
