#!/usr/bin/env python3
# pyright: reportMissingImports=false, reportUndefinedVariable=false
"""
Pairwise Kosman–Leonard (2005) genetic distance from PLINK dosage data.

For diploid biallelic SNPs encoded as alternate-allele dosages (0, 1, 2), the
per-locus dissimilarity is 0 (identical genotypes), 1 (homozygotes for
different alleles), or 0.5 (one heterozygote and one homozygote, or two
heterozygotes that are not identical — the latter is only 0 vs 1, 1 vs 2).

The multilocus distance for a pair is the mean of per-locus values over loci
where both individuals are non-missing (no mean-imputation across samples).

See Kosman & Leonard, Molecular Ecology 2005; PopGenReport::gd.kosman for R.
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
print("Kosman–Leonard genetic distance (biallelic dosages)")
print("=" * 80)
print(f"Input BED: {bed_path}")
print(f"Input BIM: {bim_path}")
print(f"Input FAM: {fam_path}")
print(f"Output: {output_path}")
print()

print("Reading PLINK genotype matrix...")
bed = open_bed(bed_path)
genotypes = bed.read(dtype=np.float32)
n, L = genotypes.shape
print(f"Genotype matrix: {n} samples x {L} SNPs")

sample_ids = bed.iid
if sample_ids.ndim == 2 and sample_ids.shape[1] >= 2:
    sample_ids = sample_ids[:, 1]
sample_ids = [str(s) for s in sample_ids]

missing = np.isnan(genotypes)
print(f"Missing genotype calls: {int(missing.sum())}")
print("Computing pairwise Kosman distances (loci with both calls only)...")

out = np.zeros((n, n), dtype=np.float64)
for i in range(n):
    sub = genotypes[i:, :]
    ni = n - i
    gi = genotypes[i, :]
    gi_b = np.broadcast_to(gi, (ni, L))
    both = ~missing[i, :] & ~missing[i:, :]

    same = (gi_b == sub) & both
    hom_diff = (((gi_b == 0) & (sub == 2)) | ((gi_b == 2) & (sub == 0))) & both

    d = np.full((ni, L), np.nan, dtype=np.float32)
    d[same] = 0.0
    d[hom_diff] = 1.0
    other = both & ~same & ~hom_diff
    d[other] = 0.5

    row = np.nanmean(d, axis=1).astype(np.float64)
    out[i, i:] = row
    out[i:, i] = row

print(f"Distance matrix shape: {out.shape[0]} x {out.shape[1]}")
print("Writing TSV...")
distance_df = pd.DataFrame(out, index=sample_ids, columns=sample_ids)
distance_df.to_csv(output_path, sep="\t", index=True, header=True)

print("Done.")
print("=" * 80)
