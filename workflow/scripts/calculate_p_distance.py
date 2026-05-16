#!/usr/bin/env python3
# pyright: reportMissingImports=false, reportUndefinedVariable=false
"""
Calculate pairwise p-distance from PLINK dosage genotypes.

For diploid biallelic SNPs encoded as alternate-allele dosages (0, 1, 2),
the per-locus difference between two individuals is:

    |g_i - g_j| / 2

Pairwise p-distance is the mean of per-locus differences over loci where both
individuals are non-missing (no cross-sample imputation).
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
print("p-distance (biallelic dosage)")
print("=" * 80)
print(f"Input BED: {bed_path}")
print(f"Input BIM: {bim_path}")
print(f"Input FAM: {fam_path}")
print(f"Output: {output_path}")
print()

print("Reading PLINK genotype matrix...")
bed = open_bed(bed_path)
genotypes = bed.read(dtype=np.float32)
n_samples, n_snps = genotypes.shape
print(f"Genotype matrix: {n_samples} samples x {n_snps} SNPs")

sample_ids = bed.iid
if sample_ids.ndim == 2 and sample_ids.shape[1] >= 2:
    sample_ids = sample_ids[:, 1]
sample_ids = [str(sample_id) for sample_id in sample_ids]

missing = np.isnan(genotypes)
print(f"Missing genotype calls: {int(missing.sum())}")
print("Computing pairwise p-distances (using pairwise-complete loci)...")

out = np.zeros((n_samples, n_samples), dtype=np.float64)
for i in range(n_samples):
    sub = genotypes[i:, :]
    ni = n_samples - i
    gi = genotypes[i, :]
    gi_b = np.broadcast_to(gi, (ni, n_snps))
    valid = ~missing[i, :] & ~missing[i:, :]

    d = np.full((ni, n_snps), np.nan, dtype=np.float32)
    d[valid] = np.abs(gi_b[valid] - sub[valid]) / 2.0

    row = np.nanmean(d, axis=1).astype(np.float64)
    out[i, i:] = row
    out[i:, i] = row

np.fill_diagonal(out, 0.0)

if np.isnan(out).any():
    nan_pairs = np.argwhere(np.isnan(out))
    n_nan = nan_pairs.shape[0]
    preview = ", ".join(
        f"({sample_ids[a]}, {sample_ids[b]})" for a, b in nan_pairs[:5]
    )
    raise ValueError(
        "p-distance contains NaN values (pairs with no overlapping non-missing loci). "
        f"NaN entries: {n_nan}. Example pairs: {preview}"
    )

print(f"Distance matrix shape: {out.shape[0]} x {out.shape[1]}")
print("Writing TSV...")
distance_df = pd.DataFrame(out, index=sample_ids, columns=sample_ids)
distance_df.to_csv(output_path, sep="\t", index=True, header=True)

print("Done.")
print("=" * 80)
