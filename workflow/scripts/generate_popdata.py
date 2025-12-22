#!/usr/bin/env python

import pandas as pd
import sys

# Snakemake automatically passes input/output/log via the 'snakemake' object
popdata_path = snakemake.params.popdata
popmap_path = snakemake.input.popmap
output_path = snakemake.output.indpopdata
log_file = snakemake.log[0]

# Redirect stdout and stderr to log file
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Reading popdata file: {popdata_path}")
# Read the provided popdata file
popdata_df = pd.read_csv(popdata_path, sep="\t", header=0)
print(f"  Found {len(popdata_df)} unique sites in popdata")

print(f"Reading popmap file: {popmap_path}")
# Read the generated popmap file
popmap_df = pd.read_csv(popmap_path, sep="\t", header=None, names=["Ind", "Site"])
print(f"  Found {len(popmap_df)} individuals in popmap")
print(f"  Found {popmap_df['Site'].nunique()} unique sites in popmap")

# Check if all sites in popmap exist in popdata
unique_sites_popmap = set(popmap_df["Site"].unique())
unique_sites_popdata = set(popdata_df["Site"].unique())
missing_sites = unique_sites_popmap - unique_sites_popdata

if missing_sites:
    print("\n" + "=" * 80)
    print("ERROR: Sites found in popmap but missing from popdata file")
    print("=" * 80)
    print(f"\nPopdata file: {popdata_path}")
    print(f"Popmap file:  {popmap_path}")
    print(f"\nMissing {len(missing_sites)} site(s):\n")

    for site in sorted(missing_sites):
        # Find all individuals with this site
        individuals = popmap_df[popmap_df["Site"] == site]["Ind"].tolist()
        print(f"  Site: '{site}'")
        print(f"    Affected individuals ({len(individuals)}):")
        for ind in individuals[:5]:  # Show first 5
            print(f"      - {ind}")
        if len(individuals) > 5:
            print(f"      ... and {len(individuals) - 5} more")
        print()

    print("=" * 80)
    print("Please add these sites to the popdata file or update the popmap.")
    print("=" * 80)
    sys.exit(1)

# Merge to ensure consistency
merged_df = pd.merge(popmap_df, popdata_df, left_on="Site", right_on="Site", how="left")

# Rename the 'Individual_x' column back to 'Individual' (if exists)
merged_df = merged_df.rename(columns={"Ind_x": "Ind"})

# Select relevant columns: individual, population, latitude, longitude, and any additional columns
output_columns = ["Ind", "Site"] + list(merged_df.columns[2:])
final_popdata_df = merged_df[output_columns]

# Save the final popdata file
final_popdata_df.to_csv(output_path, sep="\t", header=True, index=False)