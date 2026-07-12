from pathlib import Path

import pandas as pd


ROOT = Path(__file__).parents[1]


def test_final_mapping_has_all_22_clusters_and_14_labels():
    mapping = pd.read_csv(ROOT / "metadata/annotation_mapping.tsv", sep="\t")
    assert set(mapping["leiden_0_4"].astype(int)) == set(range(22))
    assert mapping["cell_type"].nunique() == 14
