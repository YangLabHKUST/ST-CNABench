# Reference Data

This directory stores reference files used by ST-CNABench.

## Included In Git

Small hg38 genome annotation files are bundled under `refs/hg38_genome_info/`:

- `hg38_genes_simple.txt`
- `hg38_genes_annot.txt`
- `cytoBand.txt.gz`
- `hg38.list`

## Download Separately

Large population phasing resources are not bundled in git. They are required only for allele-aware wrappers:

- `CalicoST`
- `Numbat`
- `Xclone`

Download link: [Google Drive](https://drive.google.com/file/d/12-hEUoDdaXTdap-Ro4Ekx4XCSlbYGP_X/view?usp=drive_link)

After download, place the extracted `population_phasing/` directory under `refs/`:

```text
refs/
└── population_phasing/
```

Required files and layout are documented in [population_phasing/README.md](population_phasing/README.md).
