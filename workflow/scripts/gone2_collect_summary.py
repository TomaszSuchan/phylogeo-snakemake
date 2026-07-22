#!/usr/bin/env python3
"""Summarise per-population GONE2 outputs."""

import sys
from pathlib import Path

import pandas as pd

project = snakemake.wildcards.project
populations_path = Path(snakemake.input["populations"])
prep_dir = Path(snakemake.params["prep_dir"])
out_dir = Path(snakemake.params["out_dir"])
summary_path = Path(snakemake.output["summary"])
log_path = Path(snakemake.log[0])

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

pop_df = pd.read_csv(populations_path, sep="\t", dtype=str)
rows = []
for _, row in pop_df.iterrows():
    pop = row["pop"]
    vcf_path = prep_dir / f"{project}.{pop}.vcf"
    ne_path = out_dir / f"{project}.{pop}_GONE2_Ne"
    if not ne_path.is_file():
        legacy = list(out_dir.glob(f"{project}.{pop}_GONE_Ne"))
        if not legacy:
            legacy = list((out_dir / "vcf").glob(f"{project}.{pop}_GONE2_Ne"))
        ne_path = legacy[0] if legacy else ne_path
    d2 = list(out_dir.glob(f"{project}.{pop}_GONE2_d2"))
    if not d2:
        d2 = list((out_dir / "vcf").glob(f"{project}.{pop}_GONE2_d2"))
    stats = list(out_dir.glob(f"{project}.{pop}_GONE2_STATS"))
    if not stats:
        stats = list((out_dir / "vcf").glob(f"{project}.{pop}_GONE2_STATS"))
    ne_summary = out_dir / f"{project}.{pop}_GONE2_Ne_summary.tsv"
    ne_seeds = out_dir / f"{project}.{pop}_GONE2_Ne_seeds.tsv"
    chrom_filter = prep_dir / f"{project}.{pop}.chrom_filter.tsv"
    rows.append(
        {
            "population": row["population"],
            "pop": pop,
            "n_individuals": row.get("n_individuals", ""),
            "vcf_file": str(vcf_path) if vcf_path.is_file() else "",
            "ne_file": str(ne_path) if Path(ne_path).is_file() else "",
            "ne_summary_file": str(ne_summary) if ne_summary.is_file() else "",
            "ne_seeds_file": str(ne_seeds) if ne_seeds.is_file() else "",
            "d2_file": str(d2[0]) if d2 else "",
            "stats_file": str(stats[0]) if stats else "",
            "chrom_filter_file": str(chrom_filter) if chrom_filter.is_file() else "",
            "status": "ok" if Path(ne_path).is_file() else "missing_ne_output",
        }
    )

out_df = pd.DataFrame(rows)
out_df.to_csv(summary_path, sep="\t", index=False)
print(out_df.to_string(index=False))
