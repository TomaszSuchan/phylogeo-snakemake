#!/usr/bin/env python
"""
Extract sample names for two populations from indpopdata based on a specified column.
"""

import pandas as pd
import sys

# Snakemake automatically passes input/output/log via the 'snakemake' object
indpopdata_file = snakemake.input["indpopdata"]
output_file = snakemake.output["samples_file"]
pop_column = snakemake.params["pop_column"]
pop1 = snakemake.params["pop1"]
pop2 = snakemake.params["pop2"]
log_file = snakemake.log[0]

# Redirect stdout and stderr to log file
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print("=" * 80)
print("Genome Scan Sample Extraction")
print("=" * 80)
print(f"Input indpopdata: {indpopdata_file}")
print(f"Population column: {pop_column}")
print(f"Population 1: {pop1}")
print(f"Population 2: {pop2}")
print(f"Output samples file: {output_file}")
print()

# Read indpopdata
print("Reading indpopdata file...")
indpopdata = pd.read_csv(indpopdata_file, sep="\t", header=0)
print(f"  Loaded {len(indpopdata)} individuals")
print(f"  Columns: {', '.join(indpopdata.columns)}")
print()

# Check required columns
if "Ind" not in indpopdata.columns:
    print("ERROR: 'Ind' column not found in indpopdata")
    print(f"Available columns: {', '.join(indpopdata.columns)}")
    sys.exit(1)

if pop_column not in indpopdata.columns:
    print(f"ERROR: Population column '{pop_column}' not found in indpopdata")
    print(f"Available columns: {', '.join(indpopdata.columns)}")
    sys.exit(1)

# Extract samples for each population
print(f"Extracting samples from column '{pop_column}'...")
pop1_samples = indpopdata[indpopdata[pop_column] == pop1]["Ind"].tolist()
pop2_samples = indpopdata[indpopdata[pop_column] == pop2]["Ind"].tolist()

print(f"  Population 1 ({pop1}): {len(pop1_samples)} samples")
print(f"  Population 2 ({pop2}): {len(pop2_samples)} samples")
print()

# Check if populations were found
if len(pop1_samples) == 0:
    print(f"ERROR: No samples found for population '{pop1}'")
    print(f"Available values in column '{pop_column}': {', '.join(sorted(indpopdata[pop_column].dropna().unique()))}")
    sys.exit(1)

if len(pop2_samples) == 0:
    print(f"ERROR: No samples found for population '{pop2}'")
    print(f"Available values in column '{pop_column}': {', '.join(sorted(indpopdata[pop_column].dropna().unique()))}")
    sys.exit(1)

# Combine samples
all_samples = pop1_samples + pop2_samples

# Write samples to file (one per line)
with open(output_file, 'w') as f:
    for sample in all_samples:
        f.write(f"{sample}\n")

print(f"Successfully wrote {len(all_samples)} samples to {output_file}")
print(f"  Population 1: {len(pop1_samples)} samples")
print(f"  Population 2: {len(pop2_samples)} samples")

