"""Configuration loading and validation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path) -> dict[str, Any]:
    """Load the analysis YAML and validate required top-level sections."""
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    required = {
        "input",
        "qc",
        "normalization",
        "dimensionality",
        "clustering",
        "differential_expression",
        "enrichment",
        "copykat",
        "annotation",
    }
    missing = required.difference(config or {})
    if missing:
        raise ValueError(f"Missing configuration sections: {sorted(missing)}")
    return config
