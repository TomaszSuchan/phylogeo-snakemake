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
- Ambiguous/IUPAC (e.g., R, Y, S...) and N are treated as missing.
- Variable site: among valid nucleotides (A,C,G,T) there are >= 2 distinct alleles.
- Invariant site: among valid nucleotides there is exactly 1 allele.

Options:
- --samples-file: optional; if omitted, all samples in the loci file are included.
  If given, one sample name per line (only these samples are written to the VCF).
- --template-vcf: ipyrad VCF(.gz) used as the exact template for polymorphic
  records (CHROM/POS/REF/ALT). This guarantees polymorphic site congruence.

VCF encoding:
- Polymorphic records: REF/ALT are taken exactly from --template-vcf
- Invariant records: REF is chosen as majority valid allele; ALT is "."
- GT per sample:
    - A/C/G/T -> 0/0 if REF, k/k if matches ALT[k]
    - Ambiguous (IUPAC) or N -> ./.
    - N or missing -> ./.
"""

import sys
import re
import argparse
import gzip
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

def parse_template_vcf(vcf_path):
    """
    Parse template VCF and return:
    - variants: {(chrom, pos): {"ref": REF, "alt": [ALT...]}}
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    variants = {}
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _vid, ref, alt = parts[:5]
            alt_list = [a for a in alt.split(",") if a and a != "."]
            variants[(chrom, int(pos))] = {"ref": ref, "alt": alt_list}
    return variants

def choose_ref_allele(counts):
    """
    Choose REF allele:
    Allele with highest count.
    """
    if not counts:
        return None
    return max(counts.items(), key=lambda kv: kv[1])[0]

def summarize_site(sequences, pos):
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
    """
    counts = defaultdict(int)
    ns_valid = 0
    for seq in sequences.values():
        if pos < len(seq):
            b = seq[pos].upper()
            if b in VALID:
                counts[b] += 1
                ns_valid += 1
            # else: missing/other characters ignored for counts

    if not counts:
        return None

    alleles = sorted(counts.keys())
    # Position is invariant if all valid nucleotides are the same.
    is_invariant = len(alleles) == 1
    ref = choose_ref_allele(counts)
    alt = [a for a in alleles if a != ref]
    return {
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "counts": counts,
        "ns_valid": ns_valid,
        "is_invariant": is_invariant,
    }

def encode_gt_for_sample(base, ref, alt):
    """
    Map a per-sample base to VCF GT:
    - If base == ref -> 0/0
    - If base == ALT[k] -> k+1/k+1
    - Any non-ACGT (including IUPAC and N) -> ./.
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
    # N/IUPAC/unknown
    return "./."

def generate_vcf(loci_data, output_file, template_variants, samples_filter=None):
    """
    Write a VCF including invariant and/or polymorphic sites based on 'mode'.
    
    Args:
        loci_data: List of dicts per locus, each with 'sequences', 'chrom', 'start', 'end'
        output_file: Output VCF file path
        template_variants: Dict[(chrom, pos) -> {"ref": str, "alt": [str,...]}]
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

    emitted_records = 0
    emitted_invariant = 0
    emitted_polymorphic = 0
    skipped_no_summary = 0
    loci_polymorphic_not_in_template = 0
    seen_template_positions = set()
    template_only_rows_written = 0

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
                summary = summarize_site(sequences, pos)
                if summary is None:
                    skipped_no_summary += 1
                    continue  # no valid alleles at this position

                # Reference-mapped: START is 0-based (BED-like), so 1-based VCF POS = start+1+pos
                # de novo: 1-based position within locus
                vcf_pos = (locus_start + pos + 1) if locus_start is not None else (pos + 1)
                var_id = f"loc{locus_idx}_pos{pos}"
                key = (chrom, vcf_pos)

                # Polymorphic sites are taken exactly from template VCF
                if key in template_variants:
                    is_inv = False
                    ref = template_variants[key]["ref"]
                    alt = list(template_variants[key]["alt"])
                    seen_template_positions.add(key)
                else:
                    # If loci says polymorphic but template does not, keep template congruence:
                    # do not emit as polymorphic.
                    if not summary["is_invariant"]:
                        loci_polymorphic_not_in_template += 1
                        continue
                    is_inv = True
                    ref = summary["ref"]
                    alt = []

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
                    gt = encode_gt_for_sample(base, ref, alt)
                    vcf.write(f"\t{gt}")
                vcf.write("\n")
                emitted_records += 1
                if is_inv:
                    emitted_invariant += 1
                else:
                    emitted_polymorphic += 1

        # Emit template variants that were not encountered in loci iteration
        # (rare edge-case). This preserves exact polymorphic concordance.
        missing_template_positions = set(template_variants.keys()) - seen_template_positions
        for chrom, vcf_pos in sorted(missing_template_positions):
            tmpl = template_variants[(chrom, vcf_pos)]
            var_id = f"template_{chrom}_{vcf_pos}"
            ref = tmpl["ref"]
            alt = list(tmpl["alt"])
            alt_field = ",".join(alt) if alt else "."
            info = "POLYMORPHIC;NS=0"
            vcf.write(f"{chrom}\t{vcf_pos}\t{var_id}\t{ref}\t{alt_field}\t.\tPASS\t{info}\tGT")
            for _s in samples_order:
                vcf.write("\t./.")
            vcf.write("\n")
            emitted_records += 1
            emitted_polymorphic += 1
            template_only_rows_written += 1

    return {
        "emitted_records": emitted_records,
        "emitted_invariant": emitted_invariant,
        "emitted_polymorphic": emitted_polymorphic,
        "skipped_no_summary": skipped_no_summary,
        "template_variants_total": len(template_variants),
        "template_variants_seen_in_loci": len(seen_template_positions),
        "template_variants_missing_from_loci": len(set(template_variants.keys()) - seen_template_positions),
        "loci_polymorphic_not_in_template": loci_polymorphic_not_in_template,
        "template_only_rows_written": template_only_rows_written,
        "loci_after_filter": len(loci_data),
        "samples_after_filter": len(samples_order),
    }

def main():
    parser = argparse.ArgumentParser(
        description=("Extract invariant and variable sites from ipyrad .loci to VCF.\n"
                     "Polymorphic records are taken exactly from --template-vcf for congruence.")
    )
    parser.add_argument("input_file", help="Input ipyrad .loci file")
    parser.add_argument("-o", "--output", default="sites.vcf", help="Output VCF path (default: sites.vcf)")
    parser.add_argument(
        "--samples-file",
        default=None,
        metavar="PATH",
        help="Optional. Restrict output to these samples (one name per line). "
        "Default: include every sample present in the loci file.",
    )
    parser.add_argument("--template-vcf", required=True, help="Input ipyrad VCF(.gz) used as polymorphic template")

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
    template_variants = parse_template_vcf(args.template_vcf)
    generate_vcf(
        loci_data,
        output_file=args.output,
        template_variants=template_variants,
        samples_filter=samples_filter,
    )
    print(f"VCF written: {args.output}")

if __name__ == "__main__":
    main()
