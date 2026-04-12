# External Tools And Runtime Notes

Some wrappers depend on external source trees or helper binaries outside the Python package itself.

## Wrappers That Need External Tools

| Model | External tools needed |
| --- | --- |
| `CalicoST` | `CalicoST`, `Eagle_v2.4.1` |
| `Clonalscope_NoWGS` | `clonalscope` |
| `Clonalscope_WGS` | `clonalscope` |
| `Numbat` | `numbat`, `Eagle_v2.4.1` |
| `Xclone` | `Eagle_v2.4.1` |

## Install Commands

From the repository root, run:

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

Expected layout:

```text
external_tools/
├── CalicoST/
├── Eagle_v2.4.1/
├── clonalscope/
└── numbat/
```

## Default Config Mapping

If you keep the default public layout, you usually do not need to change the external tool paths.
The public model templates use:

```yaml
project_settings:
  exts: "${root_dir}/external_tools"
```

With that setting, the default paths resolve as follows:

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

## Runtime Notes

- `conda`: easiest for local debugging and direct inspection
- `docker`: easiest for reproducible containerized runs
- `apptainer`: usually the best fit on HPC systems

## UMI / Cell Tag Reminder

For allele-aware wrappers such as `Numbat` and `Xclone`, make sure `UMItag` and `cellTAG` match your BAM tagging convention.

Examples:

- some STpipeline outputs use `UMItag: B3` and `cellTAG: B0`
- some Slide-seq WARP outputs use `UMItag: UR` and `cellTAG: CB`

See [Model Run](model_run.md) for the model-side configuration fields.
