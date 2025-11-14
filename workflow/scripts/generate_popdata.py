#!/usr/bin/env python

import pandas as pd
import sys

# Snakemake automatically passes input/output/log via the 'snakemake' object
popdata_path = snakemake.params.popdata
popmap_path = snakemake.input.popmap
output_path = snakemake.output.indpopdata

# Read the provided popdata file
popdata_df = pd.read_csv(popdata_path, sep="\t", header=0)

# Read the generated popmap file
popmap_df = pd.read_csv(popmap_path, sep="\t", header=None, names=["Ind", "Site"])

# Merge to ensure consistency
merged_df = pd.merge(popmap_df, popdata_df, left_on="Site", right_on="Site", how="left")

# Rename the 'Individual_x' column back to 'Individual' (if exists)
merged_df = merged_df.rename(columns={"Ind_x": "Ind"})

# Select relevant columns: individual, population, latitude, longitude, and any additional columns
output_columns = ["Ind", "Site"] + list(merged_df.columns[2:])
final_popdata_df = merged_df[output_columns]

# Save the final popdata file
final_popdata_df.to_csv(output_path, sep="\t", header=True, index=False)