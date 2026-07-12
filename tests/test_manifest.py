from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).parents[1]


def test_sample_manifest_matches_reference_counts():
    samples = pd.read_csv(ROOT / "metadata/samples.tsv", sep="\t")
    assert len(samples) == 46
    assert samples["project"].nunique() == 4
    assert samples["cells_after_qc"].sum() == 305_210


def test_thesis_parameters_are_explicit():
    cfg = yaml.safe_load((ROOT / "config/analysis.yaml").read_text())
    assert cfg["qc"]["min_cells"] == 3
    assert cfg["qc"]["min_genes"] == 100
    assert cfg["qc"]["max_pct_mito"] == 40.0
    assert cfg["qc"]["nmads"] == 3.0
    assert cfg["normalization"]["n_top_genes"] == 2000
    assert cfg["dimensionality"]["n_pcs"] == 50
    assert cfg["clustering"]["resolution"] == 0.4
