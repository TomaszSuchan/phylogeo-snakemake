#!/usr/bin/env python3
# pyright: reportMissingImports=false, reportUndefinedVariable=false
"""
Calculate a true per-site p-distance from the all-sites (variant + invariant) VCF.

For diploid genotypes coded as alternate-allele dosage (0, 1, 2), the per-site
difference between two individuals is |g_i - g_j| / 2. The p-distance is the mean
of this difference over all sites where both individuals are non-missing, counting
BOTH variant and invariant sites. Invariant sites contribute 0 to the numerator and
1 to the denominator, which normalises the SNP-only dissimilarity into a proper
per-site p-distance whose magnitude does not depend on SNP ascertainment
(MAC/MAF/thinning).

Input is the reconstructed all-sites VCF (the same one pixy uses), read with vcfpy.
Multiallelic sites are skipped (no well-defined biallelic dosage); biallelic and
invariant sites are used. No cross-sample imputation (pairwise-complete sites only).
"""

import re
import sys

import numpy as np
import pandas as pd
import vcfpy


vcf_path = snakemake.input["vcf"]
output_path = snakemake.output["dist"]
log_path = snakemake.log[0]

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

print("=" * 80)
print("p-distance (true per-site, all-sites VCF)")
print("=" * 80)
print(f"Input VCF: {vcf_path}")
print(f"Output: {output_path}")
print()


def parse_dosage(gt):
    """Alternate-allele dosage from a GT string, or NaN if missing/unparseable."""
    if gt is None:
        return np.nan
    alleles = re.split(r"[/|]", str(gt))
    if not alleles or any(a == "." for a in alleles):
        return np.nan
    return float(sum(1 for a in alleles if a != "0"))


print("Reading all-sites VCF and building dosage matrix...")
reader = vcfpy.Reader.from_path(vcf_path)
sample_ids = list(reader.header.samples.names)
n_samples = len(sample_ids)
sample_index = {s: i for i, s in enumerate(sample_ids)}
print(f"Samples: {n_samples}")

dosage_cols = []          # one (n_samples,) array per retained site
n_total = 0
n_multiallelic = 0
for record in reader:
    n_total += 1
    # Keep invariant (no ALT) and biallelic (one ALT); skip multiallelic.
    if len(record.ALT) > 1:
        n_multiallelic += 1
        continue
    col = np.full(n_samples, np.nan, dtype=np.float32)
    for call in record.calls:
        j = sample_index.get(call.sample)
        if j is None:
            continue
        col[j] = parse_dosage(call.data.get("GT"))
    dosage_cols.append(col)
reader.close()

n_sites = len(dosage_cols)
print(f"Total records: {n_total}")
print(f"Retained sites (biallelic + invariant): {n_sites}")
print(f"Skipped multiallelic sites: {n_multiallelic}")
if n_sites == 0:
    raise ValueError("No usable (biallelic or invariant) sites found in the all-sites VCF.")

# genotypes: n_samples x n_sites, NaN for missing.
genotypes = np.ascontiguousarray(np.array(dosage_cols, dtype=np.float32).T)
del dosage_cols
missing = np.isnan(genotypes)
print(f"Missing genotype calls: {int(missing.sum())}")
print("Computing pairwise p-distances (pairwise-complete sites)...")

out = np.zeros((n_samples, n_samples), dtype=np.float64)
for i in range(n_samples):
    sub = genotypes[i:, :]
    ni = n_samples - i
    gi = genotypes[i, :]
    gi_b = np.broadcast_to(gi, (ni, n_sites))
    valid = ~missing[i, :] & ~missing[i:, :]

    d = np.full((ni, n_sites), np.nan, dtype=np.float32)
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
        "p-distance contains NaN values (pairs with no overlapping non-missing sites). "
        f"NaN entries: {n_nan}. Example pairs: {preview}"
    )

print(f"Distance matrix shape: {out.shape[0]} x {out.shape[1]}")
print("Writing TSV...")
distance_df = pd.DataFrame(out, index=sample_ids, columns=sample_ids)
distance_df.to_csv(output_path, sep="\t", index=True, header=True)

print("Done.")
print("=" * 80)
