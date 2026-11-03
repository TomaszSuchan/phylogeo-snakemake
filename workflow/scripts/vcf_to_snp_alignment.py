#!/usr/bin/env python3
"""Build FASTA/PHYLIP alignments from biallelic SNP genotypes in a VCF.

Diploid genotypes are encoded as single characters for RapidNJ:
  0/0 -> REF, 1/1 -> ALT, het -> IUPAC ambiguity, missing -> N.
"""

from pathlib import Path

import vcfpy

IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("A", "C")): "M",
    frozenset(("G", "T")): "K",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
}


def _encode_gt(call, ref, alt):
    if call is None or call.gt_type is None or call.gt_bases is None:
        return "N"
    bases = [b for b in call.gt_bases if b not in (None, ".")]
    if not bases:
        return "N"
    if call.is_het:
        return IUPAC.get(frozenset(bases[:2]), "N")
    b = bases[0]
    if b == ref:
        return ref
    if b == alt:
        return alt
    return "N"


def main(vcf_path, fasta_path, phy_path, summary_path):
    # Contig header warnings from ipyrad VCFs are noisy and harmless here.
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, module="vcfpy")

    reader = vcfpy.Reader.from_path(str(vcf_path))
    samples = list(reader.header.samples)
    if not samples:
        raise SystemExit(f"No samples in VCF: {vcf_path}")

    seqs = {sample: [] for sample in samples}
    n_snps = 0
    n_skipped = 0

    for record in reader:
        if not record.is_snp:
            n_skipped += 1
            continue
        if len(record.ALT) != 1 or record.ALT[0] is None:
            n_skipped += 1
            continue
        ref = record.REF
        alt = record.ALT[0].value
        if len(ref) != 1 or len(alt) != 1:
            n_skipped += 1
            continue
        calls = record.call_for_sample
        for sample in samples:
            seqs[sample].append(_encode_gt(calls.get(sample), ref, alt))
        n_snps += 1

    if n_snps == 0:
        raise SystemExit(f"No biallelic SNPs retained from {vcf_path}")

    fasta_path = Path(fasta_path)
    phy_path = Path(phy_path)
    summary_path = Path(summary_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    # Single-line FASTA (required by rapidNJ -i fa).
    with fasta_path.open("w", encoding="utf-8") as out:
        for sample in samples:
            out.write(f">{sample}\n{''.join(seqs[sample])}\n")

    # Relaxed PHYLIP.
    with phy_path.open("w", encoding="utf-8") as out:
        out.write(f"{len(samples)} {n_snps}\n")
        for sample in samples:
            out.write(f"{sample}  {''.join(seqs[sample])}\n")

    n_missing = sum(s.count("N") for s in ("".join(seqs[sample]) for sample in samples))
    total_cells = len(samples) * n_snps
    with summary_path.open("w", encoding="utf-8") as out:
        out.write("n_samples\tn_snps\tn_skipped_records\tn_missing_cells\tmissing_fraction\n")
        out.write(
            f"{len(samples)}\t{n_snps}\t{n_skipped}\t{n_missing}\t"
            f"{(n_missing / total_cells) if total_cells else 0:.6f}\n"
        )


if __name__ == "__main__":
    main(
        snakemake.input.vcf,
        snakemake.output.fasta,
        snakemake.output.phy,
        snakemake.output.summary,
    )
