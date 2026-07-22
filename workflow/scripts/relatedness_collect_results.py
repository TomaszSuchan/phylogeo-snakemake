#!/usr/bin/env python3
"""Collect per-stratum relatedness outputs into one summary table."""

import csv
import sys
from pathlib import Path

import pandas as pd

KIND_CONFIG = {
    "related_coancestry": {
        "suffix": ".related_coancestry.tsv",
        "sep": "\t",
        "header": 0,
    },
    "relatedness": {
        "suffix": ".relatedness",
        "sep": r"\s+",
        "header": 0,
        "engine": "python",
    },
    "relatedness2": {
        "suffix": ".relatedness2",
        "sep": r"\s+",
        "header": 0,
        "engine": "python",
    },
    "genome": {
        "suffix": ".genome",
        "sep": r"\s+",
        "header": 0,
        "engine": "python",
    },
}

kind = snakemake.params["kind"]
if kind not in KIND_CONFIG:
    raise ValueError(f"Unsupported collect kind: {kind}")

cfg = KIND_CONFIG[kind]
project = snakemake.wildcards.project
populations_path = Path(snakemake.input["populations"])
output_path = Path(snakemake.output["summary"])
log_path = Path(snakemake.log[0])

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout

with populations_path.open(newline="") as handle:
    populations = list(csv.DictReader(handle, delimiter="\t"))

results_dir = populations_path.parent  # results/{project}/relatedness_by_pop
frames = []
for row in populations:
    stratum = row["pop"]
    path = results_dir / f"{project}.{stratum}{cfg['suffix']}"
    if not path.exists():
        print(f"WARNING: missing {path}")
        continue

    read_kwargs = {
        "sep": cfg["sep"],
        "header": cfg["header"],
        "comment": "#",
    }
    if "engine" in cfg:
        read_kwargs["engine"] = cfg["engine"]

    df = pd.read_csv(path, **read_kwargs)
    df.columns = [str(col).strip() for col in df.columns]
    if "population" not in df.columns:
        df.insert(0, "population", row["population"])
    if "stratum" not in df.columns:
        df.insert(1, "stratum", stratum)
    frames.append(df)
    print(f"Loaded {len(df)} rows from {path.name}")

if not frames:
    print(f"ERROR: no {kind} outputs found to collect")
    sys.exit(1)

summary = pd.concat(frames, ignore_index=True)
output_path.parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(output_path, sep="\t", index=False)
print(f"Wrote {len(summary)} rows to {output_path}")
