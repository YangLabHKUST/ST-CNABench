# Reference Data

ST-CNABench uses two kinds of reference data.

## Bundled Small References

These files are already included in git under `refs/hg38_genome_info/`:

- `hg38_genes_simple.txt`
- `hg38_genes_annot.txt`
- `cytoBand.txt.gz`
- `hg38.list`

## Population Phasing Bundle

Large population phasing resources are required only for allele-aware wrappers:

- `CalicoST`
- `Numbat`
- `Xclone`

Download the bundle from:

- [Population phasing reference bundle (Google Drive)](https://drive.google.com/file/d/12-hEUoDdaXTdap-Ro4Ekx4XCSlbYGP_X/view?usp=drive_link)

After download, extract the bundle under:

```text
refs/
└── population_phasing/
```

## Required Files

The extracted `population_phasing/` directory must contain:

- `1000G_hg38/`
- `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz`
- `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz.tbi`

Expected layout:

```text
refs/
└── population_phasing/
    ├── 1000G_hg38/
    ├── genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz
    └── genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz.tbi
```

Do not rename the files or the directory.

## Practical Rule

If you are only running the shipped demo with `CopyKAT`, you do not need the population phasing bundle.
