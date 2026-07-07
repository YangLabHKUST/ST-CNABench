from __future__ import annotations

from pathlib import Path

import yaml


def get_repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def get_demo_test_output_root(tmp_path: Path) -> Path:
    return get_repo_root() / ".test_artifacts" / tmp_path.name / "processed_data"


def write_yaml(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def write_demo_data_config(tmp_path: Path) -> Path:
    repo_root = get_repo_root()
    dataset_root = repo_root / "demo_data" / "cscc_demo"
    output_root = get_demo_test_output_root(tmp_path)
    data = {
        "project_settings": {
            "dataset_root": str(dataset_root),
            "output_root": str(output_root),
        },
        "P6_vis_rep1": {
            "dataset_id": "P6_vis_rep1",
            "platform": "Visium",
            "format": "SpaceRanger",
            "genome": "hg38",
            "species": "human",
            "ref_norm": True,
            "tumor_normal_mode": "off",
            "raw": {
                "root": str(dataset_root / "raw" / "P6_vis" / "P6_vis_rep1"),
                "tumor_normal": str(
                    dataset_root
                    / "raw"
                    / "P6_vis"
                    / "P6_vis_rep1"
                    / "metadata_P6_vis_rep1_tumor_normal.tsv"
                ),
                "tumor_normal_gt": str(
                    dataset_root
                    / "raw"
                    / "P6_vis"
                    / "P6_vis_rep1"
                    / "metadata_P6_vis_rep1_tumor_normal.tsv"
                ),
                "cna_gt": str(
                    dataset_root
                    / "raw"
                    / "P6_vis"
                    / "GT"
                    / "WES_CNV_FACETS"
                    / "P6_tumor_vs_P6_normal.vcf"
                ),
                "wgs_wes_tumor_bedg": str(
                    dataset_root
                    / "raw"
                    / "P6_vis"
                    / "GT"
                    / "wes_100kb_counts"
                    / "P6_tumor.hg38.100000b.windows.counts.bedg"
                ),
                "wgs_wes_normal_bedg": str(
                    dataset_root
                    / "raw"
                    / "P6_vis"
                    / "GT"
                    / "wes_100kb_counts"
                    / "P6_normal.hg38.100000b.windows.counts.bedg"
                ),
                "bam": str(
                    dataset_root / "raw" / "P6_vis" / "P6_vis_rep1" / "possorted_genome_bam.bam"
                ),
                "bai": str(
                    dataset_root
                    / "raw"
                    / "P6_vis"
                    / "P6_vis_rep1"
                    / "possorted_genome_bam.bam.bai"
                ),
                "scalefactors": str(
                    dataset_root
                    / "raw"
                    / "P6_vis"
                    / "P6_vis_rep1"
                    / "spatial"
                    / "scalefactors_json.json"
                ),
            },
            "output": {
                "root": str(output_root / "P6_vis_rep1"),
            },
        },
    }
    return write_yaml(tmp_path / "data.yaml", data)


def create_support_tree(tmp_path: Path) -> dict[str, Path]:
    repo_root = get_repo_root()
    refs = tmp_path / "refs"
    exts = tmp_path / "exts"

    files = {
        refs / "hg38_genome_info" / "hg38_genes_simple.txt": (
            repo_root / "refs" / "hg38_genome_info" / "hg38_genes_simple.txt"
        ).read_text(),
        refs / "hg38_genome_info" / "hg38_genes_annot.txt": (
            repo_root / "refs" / "hg38_genome_info" / "hg38_genes_annot.txt"
        ).read_text(),
        refs / "hg38_genome_info" / "hg38.list": "chr17\nchr8\nchr7\nchr13\nchr10\n",
        refs / "hg38_genome_info" / "cytoBand.txt.gz": "placeholder cytoband\n",
        refs
        / "population_phasing"
        / "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz": "placeholder vcf\n",
        refs / "population_phasing" / "1000G_hg38" / "README.txt": "placeholder panel\n",
        exts / "Eagle_v2.4.1" / "eagle": "#!/bin/sh\nexit 0\n",
        exts / "Eagle_v2.4.1" / "tables" / "genetic_map_hg38_withX.txt.gz": "placeholder map\n",
        exts / "numbat" / "inst" / "bin" / "pileup_and_phase.R": "#!/usr/bin/env Rscript\n",
        exts / "BAFExtract" / "bin" / "BAFExtract": "#!/bin/sh\nexit 0\n",
        exts / "clonalscope" / "data-raw" / "README.txt": "placeholder clonalscope data\n",
        exts / "CalicoST" / "README.txt": "placeholder calicost checkout\n",
    }

    for path, content in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)

    return {"refs": refs, "exts": exts}


def write_full_model_config(tmp_path: Path, results_dir: Path, support: dict[str, Path]) -> Path:
    refs = support["refs"]
    exts = support["exts"]
    config = {
        "project_settings": {
            "root_dir": str(get_repo_root()),
            "results_dir": str(results_dir),
            "refs": str(refs),
            "exts": str(exts),
            "default_exec_mode": "docker",
        },
        "CalicoST": {
            "enabled": True,
            "docker_image": "calicost:test",
            "model_name": "CalicoST",
            "calicost_dir": str(exts / "CalicoST"),
            "eagle_dir": str(exts / "Eagle_v2.4.1"),
            "region_vcf": str(
                refs
                / "population_phasing"
                / "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz"
            ),
            "phasing_panel": str(refs / "population_phasing" / "1000G_hg38"),
            "use_tumor_purity": False,
            "n_clones": 2,
            "n_threads": 2,
            "UMItag": "Auto",
            "cellTAG": "CB",
        },
        "CopyKAT": {
            "enabled": True,
            "docker_image": "copykat:test",
            "model_name": "CopyKAT",
            "genome": "hg20",
            "win_size": 25,
            "n_cores": 2,
            "n_clones": 2,
        },
        "InferCNV": {
            "enabled": True,
            "docker_image": "infercnv:test",
            "model_name": "InferCNV",
            "gene_order_file": str(refs / "hg38_genome_info" / "hg38_genes_simple.txt"),
            "cutoff": 0.1,
            "n_threads": 2,
            "k_obs_groups": 2,
        },
        "Clonalscope_NoWGS": {
            "enabled": True,
            "docker_image": "clonalscope:test",
            "model_name": "Clonalscope_NoWGS",
            "gene_coords_file": str(refs / "hg38_genome_info" / "hg38_genes_simple.txt"),
            "aux_data_dir": str(exts / "clonalscope" / "data-raw"),
            "mincell": 2,
        },
        "Clonalscope_WGS": {
            "enabled": True,
            "docker_image": "clonalscope:test",
            "model_name": "Clonalscope_WGS",
            "gene_coords_file": str(refs / "hg38_genome_info" / "hg38_genes_simple.txt"),
            "aux_data_dir": str(exts / "clonalscope" / "data-raw"),
            "hmm_states": [0.8, 1.0, 1.2],
            "mincell": 2,
        },
        "Numbat": {
            "enabled": True,
            "docker_image": "numbat:test",
            "model_name": "Numbat",
            "n_threads": 2,
            "genome_version": "hg38",
            "pileup_script": str(exts / "numbat" / "inst" / "bin" / "pileup_and_phase.R"),
            "eagle_path": str(exts / "Eagle_v2.4.1" / "eagle"),
            "genetic_map": str(exts / "Eagle_v2.4.1" / "tables" / "genetic_map_hg38_withX.txt.gz"),
            "snp_vcf_path": str(
                refs
                / "population_phasing"
                / "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz"
            ),
            "panel_dir": str(refs / "population_phasing" / "1000G_hg38"),
            "UMItag": "Auto",
            "cellTAG": "CB",
            "n_clones": 2,
        },
        "Xclone": {
            "enabled": True,
            "docker_image": "xclone:test",
            "model_name": "Xclone",
            "n_threads": 2,
            "snp_vcf": str(
                refs
                / "population_phasing"
                / "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz"
            ),
            "gene_region": str(refs / "hg38_genome_info" / "hg38_genes_annot.txt"),
            "eagle_path": str(exts / "Eagle_v2.4.1" / "eagle"),
            "genetic_map": str(exts / "Eagle_v2.4.1" / "tables" / "genetic_map_hg38_withX.txt.gz"),
            "panel_dir": str(refs / "population_phasing" / "1000G_hg38"),
            "UMItag": "Auto",
            "cellTAG": "CB",
            "minCOUNT": 2,
            "minMAF": 0.0,
            "n_clusters": 2,
        },
        "SCEVAN": {
            "enabled": True,
            "docker_image": "scevan:test",
            "model_name": "SCEVAN",
            "n_threads": 2,
        },
        "STARCH": {
            "enabled": True,
            "docker_image": "starch:test",
            "model_name": "STARCH",
            "gene_mapping_file": str(refs / "hg38_genome_info" / "hg38_genes_simple.txt"),
            "n_clusters": 2,
            "beta_spot": 2,
            "platform": "Visium",
            "returnnormal": 1,
        },
    }
    return write_yaml(tmp_path / "models.yaml", config)


def write_smoke_model_config(tmp_path: Path, results_dir: Path) -> Path:
    config = {
        "project_settings": {
            "root_dir": str(get_repo_root()),
            "results_dir": str(results_dir),
            "refs": str(get_repo_root() / "demo_data" / "cscc_demo" / "resources"),
            "exts": str(tmp_path / "unused_exts"),
            "default_exec_mode": "conda",
        },
        "CopyKAT": {
            "enabled": True,
            "env_name": "CopyKAT",
            "docker_image": "copykat:test",
            "model_name": "CopyKAT",
            "genome": "hg20",
            "win_size": 25,
            "n_cores": 2,
            "n_clones": 2,
        },
    }
    return write_yaml(tmp_path / "smoke_models.yaml", config)


def write_smoke_eval_config(tmp_path: Path, results_dir: Path, eval_dir: Path) -> Path:
    config = {
        "project_settings": {
            "root_dir": str(get_repo_root()),
            "results_dir": str(results_dir),
            "eval_dir": str(eval_dir),
        },
        "global_params": {
            "bin_size": 100000,
            "genome_version": "hg38",
            "gene_annot_path": str(
                get_repo_root() / "refs" / "hg38_genome_info" / "hg38_genes_annot.txt"
            ),
        },
        "eval_list": {
            "CopyKAT": {
                "eval_name": ["CopyKAT"],
            }
        },
    }
    return write_yaml(tmp_path / "eval.yaml", config)
