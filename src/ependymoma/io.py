"""AnnData contracts and small report helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse


def validate_input(adata, sample_key: str, project_key: str, gene_symbols_key: str) -> None:
    """Fail early when the concatenated raw-count AnnData violates the contract."""
    missing_obs = [key for key in (sample_key, project_key) if key not in adata.obs]
    if missing_obs:
        raise KeyError(f"Missing required adata.obs columns: {missing_obs}")
    if gene_symbols_key not in adata.var:
        raise KeyError(f"Missing adata.var[{gene_symbols_key!r}]")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("AnnData observation and variable indices must be unique")
    matrix = adata.X
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if values.size and (not np.isfinite(values).all() or (values < 0).any()):
        raise ValueError("Raw counts must be finite and non-negative")
    if values.size and not np.allclose(values, np.rint(values)):
        raise ValueError("Input X is not integer-like; expected raw UMI counts")


def write_tsv(rows: Iterable[dict], path: str | Path) -> None:
    """Write a report table, creating its parent directory."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
