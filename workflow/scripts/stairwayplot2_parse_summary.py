#!/usr/bin/env python3
"""Parse Stairway Plot 2 *.final.summary into a tidy Ne trajectory TSV."""

import sys
from pathlib import Path

import pandas as pd

summary_in = Path(snakemake.input["summary"])
tsv_out = Path(snakemake.output["tsv"])

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

# Columns written by Stairway_output_summary_plot2, renamed to the Ne table schema
# shared with the GONE2 and currentNe2 plots.
rename = {
    "year": "year",
    "Ne_median": "ne_median",
    "Ne_12.5%": "ne_lower_75",
    "Ne_87.5%": "ne_upper_75",
    "Ne_2.5%": "ne_lower_95",
    "Ne_97.5%": "ne_upper_95",
    "mutation_per_site": "mutation_per_site",
    "n_estimation": "n_estimation",
    "theta_per_site_median": "theta_per_site_median",
}
out = pd.read_csv(summary_in, sep="\t")[list(rename)].rename(columns=rename)
out.insert(0, "generation", out["year"] / float(snakemake.params["year_per_generation"]))
out.to_csv(tsv_out, sep="\t", index=False)
print(f"Wrote {tsv_out} ({len(out)} rows)")
