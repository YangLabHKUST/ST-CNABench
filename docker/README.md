# Docker Runtime Images

This directory contains Dockerfile definitions for the ST-CNVBench method runtimes.

## Images Used By The Public Template

| Model | Docker image |
| --- | --- |
| `CalicoST` | `hans0410/cnv-benchmark-calicost:1.0.0` |
| `CopyKAT` | `hans0410/cnv-benchmark-copykat:1.1.0` |
| `InferCNV` | `trinityctat/infercnv:latest` |
| `Clonalscope_NoWGS` | `hans0410/cnv-benchmark-clonalscope:1.1.0` |
| `Clonalscope_WGS` | `hans0410/cnv-benchmark-clonalscope:1.1.0` |
| `Numbat` | `hans0410/cnv-benchmark-numbat:1.5.1` |
| `Xclone` | `hans0410/cnv-benchmark-xclone:0.4.0` |
| `SCEVAN` | `hans0410/cnv-benchmark-scevan:1.0.3` |
| `STARCH` | `hans0410/cnv-benchmark-starch:1.0.0` |

These names match `docker_image` in [../configs/templates/models.template.yaml](../configs/templates/models.template.yaml).

## Pull Images

Pull only the images for models you enable. Example:

```bash
docker pull hans0410/cnv-benchmark-copykat:1.1.0
```

## Rebuild Images Locally

Build from the repository root and keep the same image tag if you want to use the default template unchanged.

```bash
docker build -t hans0410/cnv-benchmark-copykat:1.1.0 docker/copykat
```

If you use a different local tag, update the corresponding `docker_image` field in your `models.yaml`.
