#!/usr/bin/env python3
"""Write a Stairway Plot 2 blueprint from the folded easySFS spectrum."""

import sys
from pathlib import Path

import pandas as pd

sfs_dir = Path(snakemake.input["sfs_dir"])
blueprint_path = Path(snakemake.output["blueprint"])
stats_path = Path(snakemake.output["stats"])
params = snakemake.params
project = snakemake.wildcards.project
stratum = snakemake.wildcards.stratum
popid = f"{project}_{stratum}"

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

L = int(Path(snakemake.input["L"]).read_text().strip())
nseq = int(Path(params["proj_file"]).read_text().strip())
dadi_sfs = sfs_dir / "dadi" / f"{stratum}-{nseq}.sfs"
print(f"Reading {dadi_sfs} (nseq={nseq}, L={L})")

header, counts = dadi_sfs.read_text().splitlines()[:2]
if "unfolded" in header.lower():
    sys.exit(f"ERROR: {dadi_sfs} is not folded")
bins = [float(x) for x in counts.split()]
if len(bins) != nseq + 1:
    sys.exit(f"ERROR: expected {nseq + 1} dadi bins in {dadi_sfs}, got {len(bins)}")

# dadi bin 0 holds the monomorphic sites (from --total-length); folded bins are 1..nseq/2.
monomorphic = bins[0]
sfs = bins[1 : nseq // 2 + 1]
if sum(sfs) <= 0:
    sys.exit(f"ERROR: no polymorphic sites in {dadi_sfs}")
print(f"polymorphic={sum(sfs):g}; monomorphic={monomorphic:g}")

nrand = params["nrand"]
if not nrand:
    # Stairway Plot 2 default: quarters of (nseq - 2) break points.
    base = max(nseq - 2, 4)
    nrand = sorted({max(1, base * k // 4) for k in (1, 2, 3, 4)})
elif isinstance(nrand, str):
    nrand = [int(x) for x in nrand.replace(",", " ").split()]

settings = {
    "popid": popid,
    "nseq": nseq,
    "L": L,
    "whether_folded": "true",
    "SFS": "\t".join(f"{x:.10g}" for x in sfs),
    "smallest_size_of_SFS_bin_used_for_estimation": 2 if params["exclude_singletons"] else 1,
    "largest_size_of_SFS_bin_used_for_estimation": nseq // 2,
    "pct_training": params["pct_training"],
    "nrand": "\t".join(str(int(x)) for x in nrand),
    "project_dir": (blueprint_path.parent / "run").resolve(),
    "stairway_plot_dir": params["stairway_plot_dir"],
    "ninput": params["ninput"],
    "random_seed": params["random_seed"],
    "mu": params["mu"],
    "year_per_generation": params["year_per_generation"],
    # Stairpainter reads the plot settings below; "0,0" lets it pick the axis ranges.
    "plot_title": popid,
    "xrange": "0,0",
    "yrange": "0,0",
    "xspacing": 2,
    "yspacing": 2,
    "fontsize": 12,
}
if settings["random_seed"] is None:
    del settings["random_seed"]

blueprint_path.write_text("".join(f"{k}: {v}\n" for k, v in settings.items()))

pd.DataFrame(
    [
        {
            "population": stratum,
            "popid": popid,
            "nseq": nseq,
            "n_diploids_projected": nseq // 2,
            "L": L,
            "polymorphic_mass": sum(sfs),
            "monomorphic_mass": monomorphic,
            "dadi_sfs": str(dadi_sfs),
            "blueprint": str(blueprint_path),
        }
    ]
).to_csv(stats_path, sep="\t", index=False)

print(f"Wrote {blueprint_path} and {stats_path}")
