#!/usr/bin/env python3
"""Summarise per-population GONE2 outputs."""

import sys
from pathlib import Path

import pandas as pd

project = snakemake.wildcards.project
populations_path = Path(snakemake.input["populations"])
vcf_dir = Path(snakemake.params["vcf_dir"])
summary_path = Path(snakemake.output["summary"])
log_path = Path(snakemake.log[0])

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

pop_df = pd.read_csv(populations_path, sep="\t", dtype=str)
rows = []
for _, row in pop_df.iterrows():
    pop = row["pop"]
    vcf_path = vcf_dir / f"{project}.{pop}.vcf"
    ne_path = vcf_dir / f"{project}.{pop}_GONE2_Ne"
    if not ne_path.is_file():
        legacy = list(vcf_dir.glob(f"{project}.{pop}_GONE_Ne"))
        ne_path = legacy[0] if legacy else ne_path
    d2 = list(vcf_dir.glob(f"{project}.{pop}_GONE2_d2"))
    stats = list(vcf_dir.glob(f"{project}.{pop}_GONE2_STATS"))
    rows.append(
        {
            "population": row["population"],
            "pop": pop,
            "n_individuals": row.get("n_individuals", ""),
            "vcf_file": str(vcf_path) if vcf_path.is_file() else "",
            "ne_file": str(ne_path) if ne_path.is_file() else "",
            "d2_file": str(d2[0]) if d2 else "",
            "stats_file": str(stats[0]) if stats else "",
            "status": "ok" if ne_path.is_file() else "missing_ne_output",
        }
    )

out_df = pd.DataFrame(rows)
out_df.to_csv(summary_path, sep="\t", index=False)
print(out_df.to_string(index=False))
