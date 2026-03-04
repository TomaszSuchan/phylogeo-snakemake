#!/usr/bin/env python3
"""
Extract invariant and variable sites from ipyrad .loci and generate VCF.

Loci file formats (both supported):
- Reference-mapped / genome alignment: locus block includes
  |N:CHROM:START-END| (e.g. |0:OY992832.1:6295-7119|), often on the line
  with //. START is 0-based (BED-like). VCF uses that CHROM and 1-based
  POS = START + 1 + offset so positions match ipyrad VCF.
- De novo assembly: locus block has only |N| (e.g. |0|, |1|). VCF uses
  synthetic CHROM (RAD_0, RAD_1, ...) and 1-based position within the locus.

Definitions:
- Valid nucleotides: A, C, G, T
- Ambiguous/IUPAC (e.g., R, Y, S...) and N are treated as missing by default.
- Variable site: among valid nucleotides (A,C,G,T) there are >= 2 distinct alleles.
- Invariant site: among valid nucleotides there is exactly 1 allele.

Options:
- --mode: output 'all' (default), only 'invariant', or only 'polymorphic'.
- --use-iupac: decode IUPAC ambiguity codes into heterozygous GT if possible.
- --ref-choice: 'majority' (default) or 'first-sample' to determine REF allele.
- --samples-file: file containing sample names to include (one per line)

VCF encoding:
- REF = chosen allele per --ref-choice
- ALT = comma-separated remaining observed alleles
- GT per sample:
    - A/C/G/T -> 0/0 if REF, k/k if matches ALT[k]
    - Ambiguous (IUPAC) -> optionally heterozygous if both alleles are in {REF}+ALT
    - N or missing -> ./.
"""

import sys
import re
import argparse
from datetime import datetime
from collections import defaultdict

# Pattern for locus reference line: |0:OY992832.1:6295-7119|
LOCUS_REF_PATTERN = re.compile(r"\|\d+:([^:]+):(\d+)-(\d+)\|")

VALID = set("ACGT")

IUPAC_TO_ALLELES = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": set(),  # treat as missing
}

def parse_loci_file(filename):
    """
    Parse ipyrad .loci file and return a list of dicts, each with:
    - 'sequences': {sample_id -> sequence}
    - 'chrom': chromosome/contig name from |N:CHROM:START-END| line, or None
    - 'start': genomic start (0-based), or None
    - 'end': genomic end, or None
    """
    loci_data = []
    current_sequences = {}
    current_meta = {"chrom": None, "start": None, "end": None}
    with open(filename, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("//") or "//" in line:
                # Reference may be on same line as // (e.g. "// ... |0:OY992832.1:6296-7119|")
                m = LOCUS_REF_PATTERN.search(line)
                if m:
                    current_meta["chrom"] = m.group(1)
                    current_meta["start"] = int(m.group(2))
                    current_meta["end"] = int(m.group(3))
                # End of current locus
                if current_sequences:
                    loci_data.append({
                        "sequences": dict(current_sequences),
                        "chrom": current_meta["chrom"],
                        "start": current_meta["start"],
                        "end": current_meta["end"],
                    })
                    current_sequences = {}
                    current_meta = {"chrom": None, "start": None, "end": None}
                continue

            # Locus reference line: |N:CHROM:START-END| (reference-mapped) or |N| (de novo)
            if line.startswith("|"):
                m = LOCUS_REF_PATTERN.search(line)
                if m:
                    current_meta["chrom"] = m.group(1)
                    current_meta["start"] = int(m.group(2))
                    current_meta["end"] = int(m.group(3))
                # If no match (e.g. |0| only), chrom/start/end stay None → VCF uses RAD_* and locus position
                continue

            # Sample line: "sample_id ACTGNN..."
            parts = line.split()
            if len(parts) >= 2:
                sample_id, sequence = parts[0], parts[1].upper()
                if any(c in "ACGTNRYSWKMBVDH" for c in sequence):
                    current_sequences[sample_id] = sequence

    # Capture last locus if file doesn't end with //
    if current_sequences:
        loci_data.append({
            "sequences": dict(current_sequences),
            "chrom": current_meta["chrom"],
            "start": current_meta["start"],
            "end": current_meta["end"],
        })
    return loci_data

def read_samples_file(filename):
    """Read sample names from a file (one per line)"""
    samples = set()
    with open(filename, 'r') as f:
        for line in f:
            sample = line.strip()
            if sample:  # Skip empty lines
                samples.add(sample)
    return samples

def choose_ref_allele(counts, sequences, pos, ref_choice, samples_order):
    """
    Choose REF allele:
    - 'majority': allele with highest count
    - 'first-sample': the base present in the first sample with a valid nucleotide
    """
    if not counts:
        return None
    if ref_choice == "majority":
        return max(counts.items(), key=lambda kv: kv[1])[0]
    elif ref_choice == "first-sample":
        for s in samples_order:
            seq = sequences.get(s)
            if seq and pos < len(seq):
                base = seq[pos].upper()
                if base in VALID:
                    return base
        # fallback to majority if none valid in first sample
        return max(counts.items(), key=lambda kv: kv[1])[0]
    else:
        raise ValueError("ref_choice must be 'majority' or 'first-sample'")

def summarize_site(sequences, pos, ref_choice, samples_order):
    """
    Return per-position summary:
    {
      'pos': int,
      'ref': str,
      'alt': [str,...],
      'counts': dict(base->count),
      'ns_valid': int,
      'is_invariant': bool
    }
    Only considers A/C/G/T as valid; ambiguous/N are ignored for allele counting.
    However, if any IUPAC ambiguity codes are present, the site is marked as variable.
    """
    counts = defaultdict(int)
    ns_valid = 0
    has_ambiguous = False
    
    # IUPAC ambiguity codes: any of these makes a site variable
    # R=purine, Y=pyrimidine, S=strong, W=weak, K=keto, M=amino,
    # B=not A, D=not C, H=not G, V=not T, N=any
    ambiguity_codes = {'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'}
    
    for seq in sequences.values():
        if pos < len(seq):
            b = seq[pos].upper()
            if b in VALID:
                counts[b] += 1
                ns_valid += 1
            elif b in ambiguity_codes:
                # Any ambiguity code indicates uncertainty and makes site variable
                has_ambiguous = True
            # else: missing/other characters ignored for counts

    if not counts:
        return None

    alleles = sorted(counts.keys())
    # Position is invariant if:
    # 1. All valid nucleotides are the same (len(alleles) == 1)
    # 2. No ambiguity codes present (any ambiguity code indicates uncertainty and makes site variable)
    is_invariant = (len(alleles) == 1) and not has_ambiguous
    ref = choose_ref_allele(counts, sequences, pos, ref_choice, samples_order)
    alt = [a for a in alleles if a != ref]
    return {
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "counts": counts,
        "ns_valid": ns_valid,
        "is_invariant": is_invariant
    }

def encode_gt_for_sample(base, ref, alt, use_iupac):
    """
    Map a per-sample base to VCF GT:
    - If base == ref -> 0/0
    - If base == ALT[k] -> k+1/k+1
    - If IUPAC and use_iupac:
        - If exactly two alleles and both in {ref} U alt, encode heterozygote a/b with indices.
        - Else -> ./.
    - Else (N/unknown) -> ./.
    """
    if base in VALID:
        if base == ref:
            return "0/0"
        try:
            idx = alt.index(base)
            return f"{idx+1}/{idx+1}"
        except ValueError:
            # Base not in ALT set (shouldn't happen if ALT == all non-REF valid alleles)
            return "./."
    if use_iupac and base in IUPAC_TO_ALLELES and len(IUPAC_TO_ALLELES[base]) == 2:
        a1, a2 = sorted(IUPAC_TO_ALLELES[base])
        # Map each to allele index if present
        def allele_to_index(a):
            if a == ref:
                return 0
            if a in alt:
                return alt.index(a) + 1
            return None

        i1 = allele_to_index(a1)
        i2 = allele_to_index(a2)
        if i1 is not None and i2 is not None:
            # sort to keep canonical representation
            i_low, i_high = min(i1, i2), max(i1, i2)
            return f"{i_low}/{i_high}"
        else:
            return "./."
    # N or unsupported ambiguous code
    return "./."

def generate_vcf(loci_data, output_file, mode="all", use_iupac=False, ref_choice="majority", samples_filter=None):
    """
    Write a VCF including invariant and/or polymorphic sites based on 'mode'.
    
    Args:
        loci_data: List of dicts per locus, each with 'sequences', 'chrom', 'start', 'end'
        output_file: Output VCF file path
        mode: 'all', 'invariant', or 'polymorphic'
        use_iupac: Whether to decode IUPAC ambiguity codes
        ref_choice: 'majority' or 'first-sample'
        samples_filter: Optional set of sample names to include
    """
    # Global sample list across all loci (stable order)
    all_samples = set()
    for locus in loci_data:
        all_samples.update(locus["sequences"].keys())
    
    # Apply filter if provided
    if samples_filter is not None:
        all_samples = all_samples.intersection(samples_filter)
    
    samples_order = sorted(all_samples)
    
    # Filter loci_data to only include samples in filter
    if samples_filter is not None:
        filtered_loci_data = []
        for locus in loci_data:
            filtered_seqs = {s: seq for s, seq in locus["sequences"].items() if s in samples_filter}
            if filtered_seqs:
                filtered_loci_data.append({
                    "sequences": filtered_seqs,
                    "chrom": locus.get("chrom"),
                    "start": locus.get("start"),
                    "end": locus.get("end"),
                })
        loci_data = filtered_loci_data

    with open(output_file, "w") as vcf:
        # Header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        vcf.write("##source=extract_invariant_vcf.py\n")
        vcf.write("##INFO=<ID=INVARIANT,Number=0,Type=Flag,Description=\"Invariant site among valid A/C/G/T\">\n")
        vcf.write("##INFO=<ID=POLYMORPHIC,Number=0,Type=Flag,Description=\"Variable site among valid A/C/G/T\">\n")
        vcf.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with valid A/C/G/T at site\">\n")
        vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

        # Column header
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in samples_order:
            vcf.write(f"\t{s}")
        vcf.write("\n")

        # Emit records locus by locus
        for locus_idx, locus in enumerate(loci_data):
            sequences = locus["sequences"]
            # Reference-mapped: use locus chrom; de novo: use synthetic RAD_0, RAD_1, ...
            chrom = locus.get("chrom") or f"RAD_{locus_idx}"
            locus_start = locus.get("start")  # 0-based genomic start (reference-mapped) or None (de novo)
            # Determine maximum sequence length in this locus
            max_len = max((len(seq) for seq in sequences.values()), default=0)
            for pos in range(max_len):
                summary = summarize_site(sequences, pos, ref_choice, samples_order)
                if summary is None:
                    continue  # no valid alleles at this position

                is_inv = summary["is_invariant"]
                ref = summary["ref"]
                alt = summary["alt"]
                if mode == "invariant" and not is_inv:
                    continue
                if mode == "polymorphic" and (is_inv or not alt):
                    # Skip invariant and sites with only IUPAC ambiguity (no actual alternate allele)
                    continue

                # Reference-mapped: START is 0-based (BED-like), so 1-based VCF POS = start+1+pos
                # de novo: 1-based position within locus
                vcf_pos = (locus_start + pos + 1) if locus_start is not None else (pos + 1)
                var_id = f"loc{locus_idx}_pos{pos}"

                alt_field = ",".join(alt) if alt else "."
                info_tags = []
                info_tags.append("INVARIANT" if is_inv else "POLYMORPHIC")
                info_tags.append(f"NS={summary['ns_valid']}")
                info = ";".join(info_tags)

                # Write the row
                vcf.write(f"{chrom}\t{vcf_pos}\t{var_id}\t{ref}\t{alt_field}\t.\tPASS\t{info}\tGT")

                # Per-sample genotypes
                for s in samples_order:
                    seq = sequences.get(s)
                    if not seq or pos >= len(seq):
                        vcf.write("\t./.")
                        continue
                    base = seq[pos].upper()
                    gt = encode_gt_for_sample(base, ref, alt, use_iupac=use_iupac)
                    vcf.write(f"\t{gt}")
                vcf.write("\n")

def main():
    parser = argparse.ArgumentParser(
        description=("Extract invariant and variable sites from ipyrad .loci to VCF.\n"
                     "Valid bases are A/C/G/T. Ambiguous/N are missing unless --use-iupac.")
    )
    parser.add_argument("input_file", help="Input ipyrad .loci file")
    parser.add_argument("-o", "--output", default="sites.vcf", help="Output VCF path (default: sites.vcf)")
    parser.add_argument("--mode", choices=["all", "invariant", "polymorphic"], default="all",
                        help="Which sites to include (default: all)")
    parser.add_argument("--use-iupac", action="store_true",
                        help="Decode IUPAC ambiguity codes into heterozygous GT if possible")
    parser.add_argument("--ref-choice", choices=["majority", "first-sample"], default="majority",
                        help="How to choose REF allele per site (default: majority)")
    parser.add_argument("--samples-file", help="File containing sample names to include (one per line)")

    args = parser.parse_args()

    # Read samples file if provided
    samples_filter = None
    if args.samples_file:
        samples_filter = read_samples_file(args.samples_file)
        if not samples_filter:
            print(f"ERROR: No samples found in filter file: {args.samples_file}", file=sys.stderr)
            sys.exit(1)

    print(f"Reading loci file: {args.input_file}")
    loci_data = parse_loci_file(args.input_file)
    print(f"Found {len(loci_data)} loci")

    # Check for critical errors
    if not loci_data:
        print(f"ERROR: No loci found in input file: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    # Collect all samples in loci
    all_samples_in_loci = set()
    for locus in loci_data:
        all_samples_in_loci.update(locus["sequences"].keys())

    if not all_samples_in_loci:
        print(f"ERROR: No samples found in loci file: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    # Check if samples in filter file are missing from loci
    if samples_filter:
        samples_in_filter_not_in_loci = samples_filter - all_samples_in_loci
        if samples_in_filter_not_in_loci:
            print(f"ERROR: {len(samples_in_filter_not_in_loci)} sample(s) in filter file not found in loci file:", file=sys.stderr)
            for sample in sorted(samples_in_filter_not_in_loci):
                print(f"  {sample}", file=sys.stderr)
            sys.exit(1)

        # Check if any samples will be included after filtering
        filtered_samples = all_samples_in_loci.intersection(samples_filter)
        if not filtered_samples:
            print(f"ERROR: No samples match between filter file and loci file", file=sys.stderr)
            sys.exit(1)

    generate_vcf(
        loci_data,
        output_file=args.output,
        mode=args.mode,
        use_iupac=args.use_iupac,
        ref_choice=args.ref_choice,
        samples_filter=samples_filter,
    )
    print(f"VCF written: {args.output}")

if __name__ == "__main__":
    main()
