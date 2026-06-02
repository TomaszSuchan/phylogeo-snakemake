# Methods documentation

A manuscript-ready methods section is **automatically generated** for each project by the pipeline.

## Output

```
results/<project>/methods.md
```

This file is produced by the `generate_methods` Snakemake rule
(`workflow/rules/methods.smk`) via the script
`workflow/scripts/generate_methods.py`.

## What it contains

- One section per **enabled** analysis (tools set to `false` in `config["projects"][…]["analyses"]` are silently omitted).
- All configuration parameters filled in automatically from `config/config.yaml`
  and `data/mainparams` (STRUCTURE BURNIN / NUMREPS).
- Sample count (`N_SAMPLES`) read from the preprocessing output
  `results/<project>/filtered_data/<project>.samples_to_keep.txt`.
- A dataset-naming note that maps human-readable terms to actual file names:

| Human-readable name | File |
|---|---|
| Unlinked biallelic SNP dataset | `biallelic_snps_thinned.vcf.gz` |
| LD-pruned biallelic SNP dataset | `biallelic_snps_ld_pruned.vcf.gz` |
| All-biallelic SNP dataset | `biallelic_snps_all.vcf.gz` |
| All-sites dataset | `allsites.vcf.gz` |

## Remaining placeholders

All values are filled in automatically. The only fall-through placeholders that can appear
if a value is genuinely unavailable at rule execution time are:

| Placeholder | Condition |
|---|---|
| `[N_SNPS]` | VCF stats file missing or unparseable |
| `[VERSION]` | Package not pinned to an exact version in any `workflow/envs/*.yaml` |

## Re-running

The file is regenerated automatically whenever the config changes or the
preprocessed sample list is updated:

```bash
snakemake results/<project>/methods.md
```
