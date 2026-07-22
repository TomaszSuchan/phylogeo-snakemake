#!/usr/bin/env python3
"""Summarise per-population currentNe2 outputs."""

import re
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


def _parse_output(path: Path) -> dict:
    """Parse currentNe2 text output for Ne point estimate and CIs.

    Prefer the whole-genome integration block when both estimates are present.
    """
    out = {
        "ne": "",
        "ne_ci50_low": "",
        "ne_ci50_high": "",
        "ne_ci90_low": "",
        "ne_ci90_high": "",
        "estimate_type": "",
        "status": "missing_output",
    }
    if not path.is_file():
        return out

    text = path.read_text(errors="replace")
    sections = re.split(r"(?=# Ne estimation)", text)
    chosen = text
    estimate_type = "unknown"
    for section in sections:
        if "integration over the whole genome" in section:
            chosen = section
            estimate_type = "whole_genome"
            break
        if "LD between chromosomes" in section:
            chosen = section
            estimate_type = "between_chromosomes"
    else:
        if "Ne point estimate" in text:
            estimate_type = "primary"

    def _value_after(label_re: str):
        m = re.search(label_re + r"\n([^\n#]+)", chosen)
        if not m:
            return ""
        val = m.group(1).strip()
        try:
            return f"{float(val):.6g}"
        except ValueError:
            return val

    ne = _value_after(r"# Ne point estimate:")
    if not ne:
        ne = _value_after(r"# N_T of the metapopulation[^:]*:")
        if ne:
            estimate_type = "metapopulation_Nt"
            out["ne_ci50_low"] = _value_after(r"# Lower 50% limit of the N_T estimate:")
            out["ne_ci50_high"] = _value_after(r"# Upper 50% limit of the N_T estimate:")
            out["ne_ci90_low"] = _value_after(r"# Lower 90% limit of the N_T estimate:")
            out["ne_ci90_high"] = _value_after(r"# Upper 90% limit of the N_T estimate:")
            out["ne"] = ne
            out["estimate_type"] = estimate_type
            out["status"] = "ok"
            return out

    out["ne"] = ne
    out["ne_ci50_low"] = _value_after(r"# Lower limit of 50% CI:")
    out["ne_ci50_high"] = _value_after(r"# Upper (?:bound|limit) of 50% CI:")
    out["ne_ci90_low"] = _value_after(r"# Lower limit of 90% CI:")
    out["ne_ci90_high"] = _value_after(r"# Upper limit of 90% CI:")
    out["estimate_type"] = estimate_type if ne else "no_ne_estimate"
    out["status"] = "ok" if ne else "missing_ne_estimate"
    return out


pop_df = pd.read_csv(populations_path, sep="\t", dtype=str)
rows = []
for _, row in pop_df.iterrows():
    pop = row["pop"]
    vcf_path = prep_dir / f"{project}.{pop}.vcf"
    out_path = out_dir / f"{project}.{pop}_currentNe2_OUTPUT.txt"
    if not out_path.is_file():
        mix = out_dir / f"{project}.{pop}_currentNe2_mix_OUTPUT.txt"
        if mix.is_file():
            out_path = mix
    parsed = _parse_output(out_path)
    rows.append(
        {
            "population": row["population"],
            "pop": pop,
            "n_individuals": row.get("n_individuals", ""),
            "vcf_file": str(vcf_path) if vcf_path.is_file() else "",
            "output_file": str(out_path) if out_path.is_file() else "",
            "estimate_type": parsed["estimate_type"],
            "ne": parsed["ne"],
            "ne_ci50_low": parsed["ne_ci50_low"],
            "ne_ci50_high": parsed["ne_ci50_high"],
            "ne_ci90_low": parsed["ne_ci90_low"],
            "ne_ci90_high": parsed["ne_ci90_high"],
            "status": parsed["status"] if out_path.is_file() else "missing_output",
        }
    )

out_df = pd.DataFrame(rows)
out_df.to_csv(summary_path, sep="\t", index=False)
print(out_df.to_string(index=False))
