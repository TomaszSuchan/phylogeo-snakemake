#!/usr/bin/env python3
# pyright: reportMissingImports=false, reportUndefinedVariable=false
"""
Calculate dosage-based Euclidean genetic distance from PLINK bed files.

The genotype matrix is read from an existing PLINK dataset and interpreted as
per-SNP alternate-allele dosage values (0, 1, 2). Missing genotypes are mean-
imputed per SNP before pairwise Euclidean distances are calculated.
"""

import sys

import numpy as np
import pandas as pd
from bed_reader import open_bed
from scipy.spatial.distance import pdist, squareform


bed_path = snakemake.input["bed"]
bim_path = snakemake.input["bim"]
fam_path = snakemake.input["fam"]
output_path = snakemake.output["dist"]
log_path = snakemake.log[0]

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

print("=" * 80)
print("Euclidean Genetic Distance")
print("=" * 80)
print(f"Input BED: {bed_path}")
print(f"Input BIM: {bim_path}")
print(f"Input FAM: {fam_path}")
print(f"Output distance matrix: {output_path}")
print()

print("Reading PLINK genotype matrix...")
bed = open_bed(bed_path)
genotypes = bed.read(dtype=np.float32)
print(f"Genotype matrix shape: {genotypes.shape[0]} samples x {genotypes.shape[1]} SNPs")

sample_ids = bed.iid
if sample_ids.ndim == 2 and sample_ids.shape[1] >= 2:
    sample_ids = sample_ids[:, 1]
sample_ids = [str(sample_id) for sample_id in sample_ids]

missing_count = int(np.isnan(genotypes).sum())
print(f"Missing genotype calls: {missing_count}")

if missing_count:
    print("Mean-imputing missing genotypes per SNP...")
    snp_means = np.nanmean(genotypes, axis=0)
    snp_means = np.where(np.isnan(snp_means), 0.0, snp_means)
    missing_rows, missing_cols = np.where(np.isnan(genotypes))
    genotypes[missing_rows, missing_cols] = snp_means[missing_cols]

print("Calculating Euclidean distance matrix...")
distance_matrix = squareform(pdist(genotypes, metric="euclidean"))
print(f"Distance matrix shape: {distance_matrix.shape[0]} x {distance_matrix.shape[1]}")

distance_df = pd.DataFrame(distance_matrix, index=sample_ids, columns=sample_ids)

print("Writing distance matrix...")
distance_df.to_csv(output_path, sep="\t", index=True, header=True)

print("Done.")
print("=" * 80)
