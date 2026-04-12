# External Tool Setup

Some model wrappers in ST-CNVBench require external tools under this directory.

## Which Models Need External Tools

| Model               | External tools needed      |
| ------------------- | -------------------------- |
| `CalicoST`          | `CalicoST`, `Eagle_v2.4.1` |
| `Clonalscope_NoWGS` | `clonalscope`              |
| `Clonalscope_WGS`   | `clonalscope`              |
| `Numbat`            | `numbat`, `Eagle_v2.4.1`   |
| `Xclone`            | `Eagle_v2.4.1`             |

## Install

Run these commands from the repository root:

```bash
mkdir -p external_tools
cd external_tools

wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar -xzf Eagle_v2.4.1.tar.gz
rm Eagle_v2.4.1.tar.gz

git clone https://github.com/raphael-group/CalicoST.git CalicoST
git clone https://github.com/seasoncloud/Clonalscope.git clonalscope
git clone https://github.com/kharchenkolab/numbat.git numbat
```

After installation, the directory should look like this:

```text
external_tools/
├── CalicoST/
├── Eagle_v2.4.1/
├── clonalscope/
└── numbat/
```

## Config

If you keep the default public layout, you usually do not need to change the external tool paths.

In `configs/templates/models.template.yaml`, the default root is:

```yaml
project_settings:
  exts: "${root_dir}/external_tools"
```

With that setting, the default paths already match this directory layout:

```yaml
CalicoST:
  calicost_dir: "${exts}/CalicoST"
  eagle_dir: "${exts}/Eagle_v2.4.1"

Clonalscope_NoWGS:
  aux_data_dir: "${exts}/clonalscope/data-raw"

Clonalscope_WGS:
  aux_data_dir: "${exts}/clonalscope/data-raw"

Numbat:
  pileup_script: "${exts}/numbat/inst/bin/pileup_and_phase.R"
  eagle_path: "${exts}/Eagle_v2.4.1/eagle"
  genetic_map: "${exts}/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"

Xclone:
  eagle_path: "${exts}/Eagle_v2.4.1/eagle"
  genetic_map: "${exts}/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
```

If you install these tools somewhere else, update `project_settings.exts` or the corresponding model-specific paths in your `models.yaml`.

## Note

For allele-based wrappers such as `Numbat`, check the BAM tag convention used by your data and update `UMItag` and `cellTAG` in `models.yaml` if needed. For example, some STpipeline outputs use `UMItag: "B3"` and `cellTAG: "B0"`, while some WARP outputs for Slide-seq use `UMItag: "UR"` and `cellTAG: "CB"`.
