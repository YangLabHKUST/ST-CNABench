from pathlib import Path

from st_cnabench.dataset import run_data_preparation

from tests.helpers import get_demo_test_output_root, write_demo_data_config


def test_demo_prep_creates_standardized_bundle(tmp_path: Path) -> None:
    data_cfg = write_demo_data_config(tmp_path)

    success = run_data_preparation(
        config_path=data_cfg,
        dataset_ids=["P6_vis_rep1"],
        overwrite=True,
    )
    assert success is True

    out_root = get_demo_test_output_root(tmp_path) / "P6_vis_rep1"
    assert (out_root / "filtered_feature_bc_matrix" / "matrix.mtx.gz").exists()
    assert (out_root / "filtered_feature_bc_matrix" / "barcodes.tsv.gz").exists()
    assert (out_root / "filtered_feature_bc_matrix" / "features.tsv.gz").exists()
    assert (out_root / "filtered_feature_bc_matrix.h5ad").exists()
    assert (out_root / "spatial" / "tissue_positions.csv").exists()
    assert (out_root / "spatial" / "scalefactors_json.json").exists()
    assert (out_root / "metadata_P6_vis_rep1_tumor_normal.tsv").exists()
    assert (out_root / "possorted_genome_bam.bam").exists()
    assert (out_root / "possorted_genome_bam.bam.bai").exists()
    assert (out_root / "wgs_wes_tumor_data.bedg").exists()
    assert (out_root / "wgs_wes_normal_data.bedg").exists()
