#!/usr/bin/env python3
"""Write per-population sample lists and a populations table from indpopdata."""

import re
import sys
from pathlib import Path

import pandas as pd

indpopdata = pd.read_csv(snakemake.input["indpopdata"], sep="\t", dtype=str)
pop_col = snakemake.params["population_column"]
min_individuals = int(snakemake.params["min_individuals"])
project = snakemake.wildcards.project
samples_dir = Path(snakemake.output["samples_dir"])
populations_path = Path(snakemake.output["populations"])
log_path = Path(snakemake.log[0])

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

samples_dir.mkdir(parents=True, exist_ok=True)
populations_path.parent.mkdir(parents=True, exist_ok=True)


def pop_safe(label: str) -> str:
    s = re.sub(r"[^\w.\-]+", "_", str(label).strip())
    return (s[:120] or "group")


rows = []
for pop_val, grp in indpopdata.groupby(pop_col, dropna=True):
    label = str(pop_val).strip()
    if not label or label.lower() == "nan":
        continue
    samples = grp["Ind"].dropna().astype(str).tolist()
    if len(samples) < min_individuals:
        print(f"Skipping {label}: {len(samples)} individuals < min_individuals={min_individuals}")
        continue
    safe = pop_safe(label)
    sample_file = samples_dir / f"{project}.{safe}.samples.txt"
    sample_file.write_text("\n".join(samples) + "\n")
    rows.append({"population": label, "pop": safe, "n_individuals": len(samples)})
    print(f"{label} -> {safe}: {len(samples)} samples")

if not rows:
    print("ERROR: no populations met min_individuals")
    sys.exit(1)

pd.DataFrame(rows).to_csv(populations_path, sep="\t", index=False)
print(f"Wrote {len(rows)} populations to {populations_path}")
