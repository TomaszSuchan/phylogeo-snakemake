#!/usr/bin/env python3
"""Write easySFS sample lists, popmaps and projections.

Modes (params.mode):
  column — group by indpopdata[population_column] (Stairway Plot 2).
  demes  — group by explicit site lists from params.demes / deme_order (moments).

Both modes write per-group files under samples_dir:
  {project}.{pop}.samples.txt
  {project}.{pop}.popmap.txt
  {project}.{pop}.proj.txt

Demes mode also writes joint inputs for a multi-population SFS:
  {project}.joint.samples.txt
  {project}.joint.popmap.txt
  {project}.joint.proj.txt
"""

import re
import sys
from pathlib import Path

import pandas as pd

indpopdata = pd.read_csv(snakemake.input["indpopdata"], sep="\t", dtype=str)
mode = str(snakemake.params.get("mode", "column")).strip().lower()
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


def project_haploids(n_samples):
    diploids = min(n_samples, int(n_diploids) if n_diploids else max_diploids)
    return 2 * diploids


def write_group_files(safe, samples, proj):
    (samples_dir / f"{project}.{safe}.samples.txt").write_text("\n".join(samples) + "\n")
    (samples_dir / f"{project}.{safe}.popmap.txt").write_text(
        "".join(f"{s}\t{safe}\n" for s in samples)
    )
    (samples_dir / f"{project}.{safe}.proj.txt").write_text(f"{proj}\n")


def deme_sites(deme_name, deme_def):
    """Accept a plain site list, or a legacy {sites: [...]} mapping."""
    if isinstance(deme_def, list):
        sites = deme_def
    elif isinstance(deme_def, dict):
        sites = deme_def.get("sites") or []
        if deme_def.get("massifs") or deme_def.get("exclude_sites"):
            sys.exit(
                f"ERROR: deme '{deme_name}' still uses massifs/exclude_sites; "
                "list sites explicitly instead"
            )
    else:
        sys.exit(f"ERROR: deme '{deme_name}' must be a list of Site names")
    sites = [str(s).strip() for s in sites if str(s).strip()]
    if not sites:
        sys.exit(f"ERROR: deme '{deme_name}' has an empty site list")
    return sites


rows = []

if mode == "column":
    pop_col = snakemake.params["population_column"]
    for pop_val, grp in indpopdata.groupby(pop_col, dropna=True):
        label = str(pop_val).strip()
        if not label or label.lower() == "nan":
            continue
        samples = grp["Ind"].dropna().astype(str).tolist()
        if len(samples) < min_individuals:
            print(
                f"Skipping {label}: {len(samples)} individuals "
                f"< min_individuals={min_individuals}"
            )
            continue
        safe = pop_safe(label)
        proj = project_haploids(len(samples))
        write_group_files(safe, samples, proj)
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

elif mode == "demes":
    deme_order = [str(x) for x in snakemake.params["deme_order"]]
    demes_cfg = snakemake.params["demes"]
    outgroup = snakemake.params.get("outgroup")
    outgroup = str(outgroup) if outgroup not in (None, "", "null", "NULL") else None
    require_n = snakemake.params.get("require_n_demes")
    if require_n is not None and len(deme_order) != int(require_n):
        sys.exit(
            f"ERROR: moments model requires exactly {require_n} demes, got {deme_order}"
        )
    if set(deme_order) != set(demes_cfg.keys()):
        sys.exit(
            f"ERROR: deme_order {deme_order} must match demes keys {list(demes_cfg.keys())}"
        )
    if outgroup is not None and outgroup not in deme_order:
        sys.exit(f"ERROR: outgroup '{outgroup}' not in deme_order {deme_order}")

    all_samples = []
    popmap_lines = []
    for deme in deme_order:
        sites = deme_sites(deme, demes_cfg[deme])
        missing = [s for s in sites if s not in set(indpopdata["Site"])]
        if missing:
            sys.exit(f"ERROR: deme '{deme}' sites not in indpopdata: {missing}")
        sub = indpopdata[indpopdata["Site"].isin(sites)]
        samples = sub["Ind"].dropna().astype(str).tolist()
        if len(samples) < min_individuals:
            sys.exit(
                f"ERROR: deme '{deme}' has {len(samples)} individuals "
                f"< min_individuals={min_individuals}"
            )
        safe = pop_safe(deme)
        proj = project_haploids(len(samples))
        write_group_files(safe, samples, proj)
        for s in samples:
            popmap_lines.append(f"{s}\t{safe}")
            all_samples.append(s)
        rows.append(
            {
                "population": deme,
                "pop": safe,
                "n_individuals": len(samples),
                "proj_haploids": proj,
                "n_sites": len(sites),
                "sites": ",".join(sites),
                "is_outgroup": (outgroup is not None and deme == outgroup),
            }
        )
        print(
            f"{deme} -> {safe}: {len(samples)} samples from {len(sites)} sites, "
            f"proj={proj} haploids; outgroup={outgroup is not None and deme == outgroup}"
        )

    seen = set()
    unique_samples = []
    for s in all_samples:
        if s in seen:
            sys.exit(f"ERROR: sample {s} assigned to more than one deme")
        seen.add(s)
        unique_samples.append(s)

    joint_samples = Path(snakemake.output["joint_samples"])
    joint_popmap = Path(snakemake.output["joint_popmap"])
    joint_proj = Path(snakemake.output["joint_proj"])
    joint_samples.write_text("\n".join(unique_samples) + "\n")
    joint_popmap.write_text("\n".join(popmap_lines) + "\n")
    joint_proj.write_text(",".join(str(r["proj_haploids"]) for r in rows) + "\n")
    print(
        f"Wrote joint inputs; outgroup={outgroup}; "
        f"proj={joint_proj.read_text().strip()}"
    )

    # moments.LD.Parsing pop file: header "sample pop", deme names (not safe labels)
    ld_pop_file = getattr(snakemake.output, "ld_pop_file", None)
    if ld_pop_file is not None:
        ld_lines = ["sample\tpop"]
        for deme in deme_order:
            sites = deme_sites(deme, demes_cfg[deme])
            for s in indpopdata.loc[indpopdata["Site"].isin(sites), "Ind"].dropna().astype(str):
                ld_lines.append(f"{s}\t{deme}")
        Path(ld_pop_file).write_text("\n".join(ld_lines) + "\n")
        print(f"Wrote LD pop file: {ld_pop_file}")

else:
    sys.exit(f"ERROR: unknown easysfs prepare mode '{mode}' (use column or demes)")

pd.DataFrame(rows).to_csv(populations_path, sep="\t", index=False)
print(f"Wrote {len(rows)} populations to {populations_path}")
