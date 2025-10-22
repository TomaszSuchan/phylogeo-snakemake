#!/usr/bin/env python

import pandas as pd
import sys

# Snakemake automatically passes input/output/log via the 'snakemake' object
popdata_path = snakemake.input.popdata
popmap_path = snakemake.input.popmap
output_path = snakemake.output.popdata

# Read the provided popdata file
popdata_df = pd.read_csv(popdata_path, sep="\t", header=0)

# Read the generated popmap file
popmap_df = pd.read_csv(popmap_path, sep="\t", header=None, names=["Individual", "Population"])

# Merge to ensure consistency
merged_df = pd.merge(popmap_df, popdata_df, left_on="Population", right_on="Population", how="left")

# Select relevant columns: individual, population, latitude, longitude, and any additional columns
output_columns = ["Individual", "Population"] + list(merged_df.columns[3:])
final_popdata_df = merged_df[output_columns]

# Rename the 'Individual_x' column back to 'Individual' (if exists)
final_popdata_df = final_popdata_df.rename(columns={"Individual_x": "Individual"})

# Save the final popdata file
final_popdata_df.to_csv(output_path, sep="\t", header=True, index=False)