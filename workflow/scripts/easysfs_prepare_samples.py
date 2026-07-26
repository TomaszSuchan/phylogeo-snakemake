#!/usr/bin/env python3
"""Write per-population sample lists, easySFS popmaps and projections, plus a populations table."""

import re
import sys
from pathlib import Path

import pandas as pd

indpopdata = pd.read_csv(snakemake.input["indpopdata"], sep="\t", dtype=str)
pop_col = snakemake.params["population_column"]
min_individuals = int(snakemake.params["min_individuals"])
n_diploids = snakemake.params["project_to_n_diploids"]
max_diploids = int(snakemake.params["max_project_diploids"])
project = snakemake.wildcards.project
samples_dir = Path(snakemake.output["samples_dir"])
populations_path = Path(snakemake.output["populations"])
log_path = Path(snakemake.log[0])

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

samples_dir.mkdir(parents=True, exist_ok=True)


def pop_safe(label):
    s = re.sub(r"[^\w.\-]+", "_", str(label).strip())
    return s[:120] or "group"


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
    diploids = min(len(samples), int(n_diploids) if n_diploids else max_diploids)
    proj = 2 * diploids  # easySFS --proj counts gene copies

    (samples_dir / f"{project}.{safe}.samples.txt").write_text("\n".join(samples) + "\n")
    (samples_dir / f"{project}.{safe}.popmap.txt").write_text(
        "".join(f"{s}\t{safe}\n" for s in samples)
    )
    (samples_dir / f"{project}.{safe}.proj.txt").write_text(f"{proj}\n")

    rows.append(
        {
            "population": label,
            "pop": safe,
            "n_individuals": len(samples),
            "proj_haploids": proj,
        }
    )
    print(f"{label} -> {safe}: {len(samples)} samples, projected to {proj} gene copies")

if not rows:
    sys.exit("ERROR: no populations met min_individuals")

pd.DataFrame(rows).to_csv(populations_path, sep="\t", index=False)
print(f"Wrote {len(rows)} populations to {populations_path}")
