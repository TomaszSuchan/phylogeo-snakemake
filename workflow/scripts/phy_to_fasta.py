#!/usr/bin/env python3
# pyright: reportMissingImports=false, reportUndefinedVariable=false
"""Convert ipyrad PHYLIP alignment to FASTA for rapidNJ."""

from Bio import SeqIO

phy_path = snakemake.input["phy"]
fasta_path = snakemake.output["fasta"]

for fmt in ("phylip", "phylip-relaxed"):
    try:
        records = list(SeqIO.parse(phy_path, fmt))
        if records:
            SeqIO.write(records, fasta_path, "fasta")
            break
    except Exception:
        records = []
else:
    raise ValueError(f"Could not parse PHYLIP alignment: {phy_path}")
