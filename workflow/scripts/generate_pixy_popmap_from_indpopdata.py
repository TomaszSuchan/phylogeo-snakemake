#!/usr/bin/env python

import os
import sys

import pandas as pd

# Snakemake inputs/outputs/params
indpopdata_path = snakemake.input["indpopdata"]
output_path = snakemake.output["popmap"]
log_path = snakemake.log[0]
group_by = snakemake.params.get("group_by")

# Redirect stdout/stderr to log
sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

print(f"Reading indpopdata file: {indpopdata_path}")
if not os.path.exists(indpopdata_path):
    print(f"ERROR: indpopdata file does not exist: {indpopdata_path}")
    sys.exit(1)

indpopdata = pd.read_csv(indpopdata_path, sep="\t", header=0)
print(f"  Loaded {len(indpopdata)} individuals")
print(f"  Columns: {', '.join(indpopdata.columns)}")

# Determine individual ID column
if "Ind" in indpopdata.columns:
    ind_col = "Ind"
elif "Sample" in indpopdata.columns:
    ind_col = "Sample"
else:
    print("ERROR: Could not find individual ID column ('Ind' or 'Sample') in indpopdata.")
    sys.exit(1)

# Determine grouping column
if isinstance(group_by, (list, tuple)):
    group_by_cols = [g for g in group_by if g and g != "none"]
    group_col = group_by_cols[0] if group_by_cols else "Site"
else:
    group_col = group_by or "Site"

print(f"\nUsing '{ind_col}' as individual column")
print(f"Using '{group_col}' as grouping column for pixy populations")

if group_col not in indpopdata.columns:
    print(f"ERROR: Grouping column '{group_col}' not found in indpopdata.")
    print(f"Available columns: {', '.join(indpopdata.columns)}")
    sys.exit(1)

# Build popmap: two columns, no header: <individual>\t<group>
popmap_df = indpopdata[[ind_col, group_col]].copy()
popmap_df.columns = ["Ind", "Pop"]  # internal names; pixy only cares about order

unique_groups = popmap_df["Pop"].nunique(dropna=True)
print(f"\nGenerated pixy popmap with {len(popmap_df)} individuals and {unique_groups} groups")

os.makedirs(os.path.dirname(output_path), exist_ok=True)
popmap_df.to_csv(output_path, sep="\t", header=False, index=False)
print(f"\nSaved pixy popmap to: {output_path}")

