#!/usr/bin/env python
"""
Create a summary table of populations with number of individuals per site.
Includes all columns from popdata (Lat, Lon, Region, etc.) plus a count of individuals.
"""

import pandas as pd
import sys

# Snakemake automatically passes input/output/log via the 'snakemake' object
indpopdata_path = snakemake.input.indpopdata
output_path = snakemake.output.summary
log_file = snakemake.log[0]

# Redirect stdout and stderr to log file
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Reading indpopdata file: {indpopdata_path}")
# Read the individual-level popdata
indpopdata_df = pd.read_csv(indpopdata_path, sep="\t", header=0)
print(f"  Found {len(indpopdata_df)} individuals")
print(f"  Found {indpopdata_df['Site'].nunique()} unique sites")
print(f"  Columns: {', '.join(indpopdata_df.columns)}")

# Count individuals per site
print("\nCounting individuals per site...")
counts_df = indpopdata_df.groupby('Site').size().reset_index(name='no_Inds')
print(f"  Created counts for {len(counts_df)} sites")

# Get all unique values per site for other columns (Lat, Lon, Region, etc.)
# For columns other than Site and Ind, take the first value (should be same for all individuals in a site)
print("\nAggregating site-level data...")
site_cols = ['Site']
other_cols = [col for col in indpopdata_df.columns if col not in ['Site', 'Ind']]

# Group by Site and take first value for each column (all should be same within a site)
site_data = indpopdata_df.groupby('Site')[other_cols].first().reset_index()

print(f"  Aggregated data for columns: {', '.join(other_cols)}")

# Merge counts with site data
summary_df = pd.merge(counts_df, site_data, on='Site', how='left')

# Reorder columns: Site, no_Inds, then all other columns
column_order = ['Site', 'no_Inds'] + other_cols
summary_df = summary_df[column_order]

print(f"\nGenerated summary table with {len(summary_df)} sites")
print(f"  Columns: {', '.join(summary_df.columns)}")
print(f"\nSummary statistics:")
print(f"  Total individuals: {summary_df['no_Inds'].sum()}")
print(f"  Mean individuals per site: {summary_df['no_Inds'].mean():.2f}")
print(f"  Min individuals per site: {summary_df['no_Inds'].min()}")
print(f"  Max individuals per site: {summary_df['no_Inds'].max()}")

# Save the summary table
summary_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"\nSaved population summary to: {output_path}")

# Close log file
sys.stdout.close()

