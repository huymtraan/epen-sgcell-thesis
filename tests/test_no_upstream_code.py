from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_excluded_upstream_directories_are_absent():
    for name in ("fetchngs", "scrnaseq-modified", "run_cmd", "cellranger_ref", "work"):
        assert not (ROOT / name).exists()
