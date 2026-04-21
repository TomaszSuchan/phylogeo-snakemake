#!/usr/bin/env python3
# pyright: reportMissingImports=false, reportUndefinedVariable=false
"""
Calculate pairwise genetic diffs using mean-allele-frequency imputation.

Implements the bed2diffs_v2-style method:
1) Replace missing genotype calls with per-SNP mean dosage values.
2) Compute similarity matrix as G G^T / n_sites.
3) Convert similarities to pairwise diffs:
   d_ij = s_ii + s_jj - 2 * s_ij
"""

import sys

import numpy as np
import pandas as pd
from bed_reader import open_bed


bed_path = snakemake.input["bed"]
bim_path = snakemake.input["bim"]
fam_path = snakemake.input["fam"]
output_path = snakemake.output["dist"]
log_path = snakemake.log[0]

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

print("=" * 80)
print("Mean-Imputed Genetic Diffs")
print("=" * 80)
print(f"Input BED: {bed_path}")
print(f"Input BIM: {bim_path}")
print(f"Input FAM: {fam_path}")
print(f"Output diffs matrix: {output_path}")
print()

print("Reading PLINK genotype matrix...")
bed = open_bed(bed_path)
genotypes = bed.read(dtype=np.float32)
n_indiv, n_sites = genotypes.shape
print(f"Genotype matrix shape: {n_indiv} samples x {n_sites} SNPs")

sample_ids = bed.iid
if sample_ids.ndim == 2 and sample_ids.shape[1] >= 2:
    sample_ids = sample_ids[:, 1]
sample_ids = [str(sample_id) for sample_id in sample_ids]

missing = np.isnan(genotypes)
missing_count = int(missing.sum())
print(f"Missing genotype calls: {missing_count}")

if missing_count:
    print("Imputing missing genotypes with per-SNP mean dosage...")
    snp_means = np.nanmean(genotypes, axis=0)
    snp_means = np.where(np.isnan(snp_means), 0.0, snp_means)
    missing_rows, missing_cols = np.where(missing)
    genotypes[missing_rows, missing_cols] = snp_means[missing_cols]

print("Computing similarities and diffs matrix...")
similarities = (genotypes @ genotypes.T) / float(n_sites)
self_similarities = np.diag(similarities)
diffs = (
    self_similarities[:, np.newaxis]
    + self_similarities[np.newaxis, :]
    - (2.0 * similarities)
)

# Numerical precision can produce tiny negatives; clamp for stability.
diffs = np.maximum(diffs, 0.0)
print(f"Diffs matrix shape: {diffs.shape[0]} x {diffs.shape[1]}")

diffs_df = pd.DataFrame(diffs, index=sample_ids, columns=sample_ids)

print("Writing diffs matrix...")
diffs_df.to_csv(output_path, sep="\t", index=True, header=True)

print("Done.")
print("=" * 80)
