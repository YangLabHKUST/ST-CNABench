# Population Phasing References

This directory contains the population phasing resources used by allele-aware wrappers.

## Used By

- `CalicoST`
- `Numbat`
- `Xclone`

## Required Files

- `1000G_hg38/`
- `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz`
- `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz.tbi`

## Expected Layout

```text
refs/
└── population_phasing/
    ├── 1000G_hg38/
    ├── genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz
    └── genome1K.phase3.SNP_AF5e2.chr1toX.hg38.ensemble_style.sorted.vcf.gz.tbi
```

Do not rename the files or directory.

## Notes

- `1000G_hg38/` is a directory, not a single file.
- `CalicoST`, `Numbat`, and `Xclone` use the common SNP VCF together with the `1000G_hg38` panel directory.
- Other public methods do not require this population phasing bundle.
