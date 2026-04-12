# Use Your Own Model

This tutorial covers two different cases.

## Case 1: Enable A Supported Public Wrapper

If the method is already supported in ST-CNVBench, you only need to configure and enable it.

1. Start from `configs/templates/models.template.yaml`
2. Copy the relevant model section into your own `models.yaml`
3. Set `enabled: true`
4. Fill in runtime-specific paths or image names
5. Run with `--models <ModelName>`

See [Model Run](../model_run.md) for the detailed model configuration contract.

## Case 2: Add A New Wrapper

If your method is not yet supported, the public extension path is:

1. Add a wrapper under `src/st_cnvbench/model/tools_scripts/<model_name>/`
2. Register it in `src/st_cnvbench/model/tools_scripts/__init__.py`
3. Add a public config section to `configs/templates/models.template.yaml`
4. If the method should participate in evaluation, add or update the corresponding evaluation loader mappings

## Recommended First Check

Before wiring a full execution path, make sure the new wrapper can at least:

- locate the standardized dataset bundle
- prepare its own input files
- build a valid command line
- fail loudly when required paths or files are missing

## Try Next

- For the packaged cSCC demo, go to [Quickstart Demo And Expected Outputs](quickstart_demo.md)
- For the CNV profile task example, go to [CNV Profile Task Example](cnv_profile_hcc2t.md)
- For the tumor-normal task example, go to [Tumor-Normal Classification Task Example](tumor_normal_hcc2t.md)
- For the subclone task example, go to [Subclone Identification Task Example](subclone_identification_slidednaseq.md)
- To prepare your own input dataset, go to [Use Your Own Dataset](use_your_own_dataset.md)
