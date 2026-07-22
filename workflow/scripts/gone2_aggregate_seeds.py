#!/usr/bin/env python3
"""Aggregate multi-seed GONE2 Ne trajectories into mean ± SD summaries."""

import io
import re
import sys
from pathlib import Path

import pandas as pd

seed_files = list(snakemake.input)
ne_mean_path = Path(snakemake.output["ne"])
summary_path = Path(snakemake.output["summary"])
seeds_path = Path(snakemake.output["seeds"])
log_path = Path(snakemake.log[0])

sys.stdout = open(log_path, "w")
sys.stderr = sys.stdout


def read_gone2_ne(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="replace").splitlines()
    data_lines = [ln for ln in lines if not ln.startswith("#")]
    if not data_lines:
        raise ValueError(f"No data lines in {path}")
    df = pd.read_csv(io.StringIO("\n".join(data_lines)), sep="\t")
    gen_col = next((c for c in df.columns if c.lower() == "generation"), None)
    ne_col = next((c for c in df.columns if c.startswith("Ne_")), None)
    if gen_col is None or ne_col is None:
        raise ValueError(f"Missing Generation / Ne_* columns in {path}: {list(df.columns)}")
    out = pd.DataFrame(
        {
            "generation": pd.to_numeric(df[gen_col], errors="coerce"),
            "ne": pd.to_numeric(df[ne_col], errors="coerce"),
        }
    )
    return out[out["generation"].notna() & out["ne"].notna() & (out["ne"] > 0)].copy()


seed_re = re.compile(r"\.seed(\d+)_GONE2_Ne$")
frames = []
for path in sorted(seed_files, key=str):
    path = Path(path)
    m = seed_re.search(path.name)
    if not m:
        raise ValueError(f"Could not parse seed from filename: {path.name}")
    seed = int(m.group(1))
    df = read_gone2_ne(path)
    df["seed"] = seed
    frames.append(df)
    print(f"Read seed {seed}: {len(df)} generations from {path}")

if not frames:
    raise SystemExit("ERROR: no seed Ne files to aggregate")

all_seeds = pd.concat(frames, ignore_index=True)
all_seeds = all_seeds.sort_values(["seed", "generation"])
all_seeds.to_csv(seeds_path, sep="\t", index=False)

summary = (
    all_seeds.groupby("generation", as_index=False)
    .agg(ne_mean=("ne", "mean"), ne_sd=("ne", "std"), n_seeds=("ne", "count"))
    .sort_values("generation")
)
summary["ne_sd"] = summary["ne_sd"].fillna(0.0)
summary.to_csv(summary_path, sep="\t", index=False)

mean_out = pd.DataFrame(
    {
        "Generation": summary["generation"].astype(int),
        "Ne_diploids": summary["ne_mean"],
    }
)
ne_mean_path.parent.mkdir(parents=True, exist_ok=True)
mean_out.to_csv(ne_mean_path, sep="\t", index=False)
print(f"Wrote mean Ne ({len(summary)} generations) to {ne_mean_path}")
print(f"Wrote seed table to {seeds_path}")
print(f"Wrote summary to {summary_path}")
