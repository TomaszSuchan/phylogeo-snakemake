#!/usr/bin/env python
"""
Create popmap file for genome scan with two populations.
"""

import pandas as pd
import sys

# Snakemake automatically passes input/output/log via the 'snakemake' object
samples_file = snakemake.input["samples_file"]
indpopdata_file = snakemake.input["indpopdata"]
output_file = snakemake.output["popmap"]
pop_column = snakemake.params["pop_column"]
pop1 = snakemake.params["pop1"]
pop2 = snakemake.params["pop2"]
log_file = snakemake.log[0]

# Redirect stdout and stderr to log file
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print("=" * 80)
print("Genome Scan Popmap Creation")
print("=" * 80)
print(f"Samples file: {samples_file}")
print(f"Indpopdata file: {indpopdata_file}")
print(f"Population column: {pop_column}")
print(f"Population 1: {pop1}")
print(f"Population 2: {pop2}")
print(f"Output popmap: {output_file}")
print()

# Read indpopdata
print("Reading indpopdata file...")
indpopdata = pd.read_csv(indpopdata_file, sep="\t", header=0)
print(f"  Loaded {len(indpopdata)} individuals")
print()

# Read samples from file
print("Reading samples from file...")
try:
    with open(samples_file, 'r') as f:
        vcf_samples = [s.strip() for s in f.readlines() if s.strip()]
except Exception as e:
    print(f"ERROR: Failed to read samples from file")
    print(f"Error: {e}")
    sys.exit(1)

print(f"  Found {len(vcf_samples)} samples in VCF")
print()

# Filter indpopdata to only include samples in VCF and the two populations
print(f"Filtering indpopdata for populations {pop1} and {pop2}...")
filtered = indpopdata[
    (indpopdata["Ind"].isin(vcf_samples)) &
    (indpopdata[pop_column].isin([pop1, pop2]))
].copy()

print(f"  Found {len(filtered)} individuals in filtered data")
print(f"  Population {pop1}: {len(filtered[filtered[pop_column] == pop1])} individuals")
print(f"  Population {pop2}: {len(filtered[filtered[pop_column] == pop2])} individuals")
print()

# Create popmap (Ind, Site format where Site is the population from pop_column)
popmap = filtered[["Ind", pop_column]].copy()
popmap.columns = ["Ind", "Site"]

# Write popmap
print(f"Writing popmap to {output_file}...")
popmap.to_csv(output_file, sep="\t", header=False, index=False)
print(f"  Wrote {len(popmap)} entries")
print(f"  Population {pop1}: {len(popmap[popmap['Site'] == pop1])} entries")
print(f"  Population {pop2}: {len(popmap[popmap['Site'] == pop2])} entries")

