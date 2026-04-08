#!/usr/bin/env python3
"""
Extract invariant and variable sites from ipyrad .loci and write one VCF.

Purpose
    Walk the .loci alignment in order. Only columns where at least one sample has A, C,
    G, or T are considered; columns with no A/C/G/T (e.g. only N, gaps, or B/D/H/V) are
    skipped with no VCF row. For each such column, classify using all samples in the locus
    (see Definitions): if invariant, emit an INVARIANT row (REF = that consensus base,
    ALT = “.”, GT from loci); if variable, emit a POLYMORPHIC row from --template-vcf
    (REF/ALT/GT as in ipyrad). Before writing, require that the set of variable columns’
    (CHROM, POS) equals the template’s variant rows. The .loci file should list the same
    samples as the VCF, or SNPs can look monomorphic in .loci and the check fails.

Algorithm (high level)
    1. Stream .loci in two passes (one locus in memory at a time) to avoid RAM blow-up on
       large alignments.
    2. Parse --template-vcf: for each row, store CHROM, POS, REF, ALT, and per-sample GT
       (GT subfield via FORMAT).
    3. Optionally restrict to --samples-file (output columns only). Classify invariant vs
       variable using all samples present in each locus (not the filter).
    4. Pass 1: union of (CHROM, POS) for every column classified as variable in .loci.
       Compare to template variant keys; if unequal, abort.
    5. Pass 2: Walk .loci columns again in the same order and write the VCF:
       a. No sample has A/C/G/T at this column → skip (no row).
       b. Classify per Definitions; invariant → INVARIANT; REF/ALT/GT from loci (code).
       c. Variable → POLYMORPHIC; REF/ALT/GT from the template row (exists after step 4).

Definitions
    Invariant site: at least one A/C/G/T in the column, all counted A/C/G/T are the same
    allele, and no sample has R, Y, S, W, K, or M at that position.
    Variable site: two or more distinct A/C/G/T among samples, or any R/Y/S/W/K/M at the
    column (IUPAC ambiguity counts as variable).

Loci formats
    Reference-mapped: |N:CHROM:START-END| on the // line; START is the VCF POS of alignment
    column 0 (matches ipyrad .vcf), so POS = START + column_index.
    De novo: |N| only → CHROM = RAD_0, RAD_1, …; POS = column_index + 1 (1-based along locus).

CLI
    --samples-file: optional; one sample name per line. Default: all samples in .loci.
    --template-vcf: required ipyrad variant VCF(.gz); must match variable sites in .loci
    exactly (bidirectional set equality on CHROM/POS).
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

# Two-base IUPAC codes: column is treated as variable (not INVARIANT), even if all A/C/G/T agree.
IUPAC_VARIABLE = frozenset("RYSWKM")

DISALLOWED_IUPAC = frozenset("BDHV")


def vcf_position(locus_start, col_index):
    """
    VCF POS (1-based) for alignment column col_index (0-based).

    Reference-mapped: START from |N:CHROM:START-END| is ipyrad's POS for column 0.
    De novo: locus_start is None → POS counts from 1 along the locus.
    """
    if locus_start is not None:
        return locus_start + col_index
    return col_index + 1


def iter_loci_file(filename):
    """
    Stream ipyrad .loci file: yield one locus dict at a time (low memory).

    Each dict has:
    - 'sequences': {sample_id -> sequence}
    - 'chrom': chromosome/contig from |N:CHROM:START-END| on the // line, or None
    - 'start': genomic start (0-based), or None
    - 'end': genomic end, or None
    """
    current_sequences = {}
    current_meta = {"chrom": None, "start": None, "end": None}
    with open(filename, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("//") or "//" in line:
                m = LOCUS_REF_PATTERN.search(line)
                if m:
                    current_meta["chrom"] = m.group(1)
                    current_meta["start"] = int(m.group(2))
                    current_meta["end"] = int(m.group(3))
                if current_sequences:
                    yield {
                        "sequences": dict(current_sequences),
                        "chrom": current_meta["chrom"],
                        "start": current_meta["start"],
                        "end": current_meta["end"],
                    }
                    current_sequences = {}
                    current_meta = {"chrom": None, "start": None, "end": None}
                continue

            if line.startswith("|"):
                m = LOCUS_REF_PATTERN.search(line)
                if m:
                    current_meta["chrom"] = m.group(1)
                    current_meta["start"] = int(m.group(2))
                    current_meta["end"] = int(m.group(3))
                continue

            parts = line.split()
            if len(parts) >= 2:
                sample_id, sequence = parts[0], parts[1].upper()
                if any(c in "ACGTNRYSWKMBVDH" for c in sequence):
                    current_sequences[sample_id] = sequence

    if current_sequences:
        yield {
            "sequences": dict(current_sequences),
            "chrom": current_meta["chrom"],
            "start": current_meta["start"],
            "end": current_meta["end"],
        }

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
    variants: {(chrom, pos): {
        "ref": REF,
        "alt": [ALT...],
        "gts": {sample_name: GT string from template},
    }}
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    variants = {}
    template_samples = []
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                template_samples = line.rstrip().split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, pos, _vid, ref, alt = parts[:5]
            fmt = parts[8]
            fmt_keys = fmt.split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                gt_idx = 0
            alt_list = [a for a in alt.split(",") if a and a != "."]
            gts = {}
            for i, sname in enumerate(template_samples):
                col = parts[9 + i] if 9 + i < len(parts) else "."
                sub = col.split(":")
                gt_val = sub[gt_idx] if gt_idx < len(sub) else "./."
                gts[sname] = gt_val
            variants[(chrom, int(pos))] = {"ref": ref, "alt": alt_list, "gts": gts}
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
    A/C/G/T are counted for allele multiplicity. N, gap, B/D/H/V are ignored for counts.
    R/Y/S/W/K/M at any sample force variable (not invariant).
    """
    counts = defaultdict(int)
    ns_valid = 0
    has_iupac_variable = False
    for seq in sequences.values():
        if pos < len(seq):
            b = seq[pos].upper()
            if b in IUPAC_VARIABLE:
                has_iupac_variable = True
            if b in VALID:
                counts[b] += 1
                ns_valid += 1

    if not counts:
        return None

    alleles = sorted(counts.keys())
    # Invariant only if a single A/C/G/T allele among counted bases and no RYSWKM in column.
    is_invariant = (len(alleles) == 1) and not has_iupac_variable
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


def _iter_locus_columns(locus, locus_idx):
    """Yield (locus_idx, chrom, locus_start, pos, sequences) for each column in one locus."""
    sequences = locus["sequences"]
    chrom = locus.get("chrom") or f"RAD_{locus_idx}"
    locus_start = locus.get("start")
    max_len = max((len(seq) for seq in sequences.values()), default=0)
    for pos in range(max_len):
        yield locus_idx, chrom, locus_start, pos, sequences


def _variable_keys_one_locus(locus, locus_idx):
    """Set of (CHROM, POS) for variable columns in a single locus."""
    keys = set()
    for _li, chrom, locus_start, pos, sequences in _iter_locus_columns(locus, locus_idx):
        summary = summarize_site(sequences, pos)
        if summary is None or summary["is_invariant"]:
            continue
        vcf_pos = vcf_position(locus_start, pos)
        keys.add((chrom, vcf_pos))
    return keys


def _format_key(k):
    chrom, pos = k
    return f"{chrom}:{pos}"


def _validate_loci_template_variable_sets(loci_var, template_variants):
    """
    Require exact equality: variable (CHROM, POS) from .loci vs variant rows in template VCF.
    """
    tmpl = set(template_variants.keys())
    only_loci = sorted(loci_var - tmpl, key=lambda x: (x[0], x[1]))
    only_vcf = sorted(tmpl - loci_var, key=lambda x: (x[0], x[1]))
    if only_loci:
        sample = ", ".join(_format_key(k) for k in only_loci[:15])
        more = f" … (+{len(only_loci) - 15} more)" if len(only_loci) > 15 else ""
        print(
            "ERROR: variable in .loci but no matching row in --template-vcf "
            f"({len(only_loci)} site(s)): {sample}{more}",
            file=sys.stderr,
        )
        sys.exit(1)
    if only_vcf:
        sample = ", ".join(_format_key(k) for k in only_vcf[:15])
        more = f" … (+{len(only_vcf) - 15} more)" if len(only_vcf) > 15 else ""
        print(
            "ERROR: variant in --template-vcf but no variable column in .loci at that "
            f"(CHROM, POS) ({len(only_vcf)} site(s)): {sample}{more}\n"
            "  Check CHROM/POS mapping or whether --samples-file made a site monomorphic.",
            file=sys.stderr,
        )
        sys.exit(1)


def encode_gt_for_sample_invariant(base, ref, warned_bdhv):
    """GT for INVARIANT rows (ALT is '.')."""
    if base in DISALLOWED_IUPAC:
        if warned_bdhv is not None and not warned_bdhv[0]:
            print(
                "WARNING: loci contain B/D/H/V; treating as missing (./.).",
                file=sys.stderr,
            )
            warned_bdhv[0] = True
        return "./."
    if base in ("N", "-", "."):
        return "./."
    if base in VALID:
        return "0/0" if base == ref else "./."
    return "./."


def collect_loci_pass1(loci_path):
    """
    Single streamed pass: variable-site keys (for template check) and union of sample IDs.
    Returns (loci_var, all_samples, n_loci).
    """
    loci_var = set()
    all_samples = set()
    n_loci = 0
    for locus_idx, locus in enumerate(iter_loci_file(loci_path)):
        loci_var.update(_variable_keys_one_locus(locus, locus_idx))
        all_samples.update(locus["sequences"].keys())
        n_loci += 1
    return loci_var, all_samples, n_loci


def write_loci_pass2(loci_path, output_file, template_variants, samples_order, n_loci):
    """Second streamed pass: emit VCF rows (caller must have validated keys vs template)."""
    emitted_records = 0
    emitted_invariant = 0
    emitted_polymorphic = 0
    skipped_no_summary = 0
    warned_bdhv = [False]

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

        for locus_idx, locus in enumerate(iter_loci_file(loci_path)):
            for _li, chrom, locus_start, pos, sequences in _iter_locus_columns(locus, locus_idx):
                summary = summarize_site(sequences, pos)
                if summary is None:
                    skipped_no_summary += 1
                    continue

                vcf_pos = vcf_position(locus_start, pos)
                var_id = f"loc{locus_idx}_pos{pos}"
                key = (chrom, vcf_pos)

                if summary["is_invariant"]:
                    ref = summary["ref"]
                    alt_field = "."
                    info = f"INVARIANT;NS={summary['ns_valid']}"
                    vcf.write(f"{chrom}\t{vcf_pos}\t{var_id}\t{ref}\t{alt_field}\t.\tPASS\t{info}\tGT")
                    for s in samples_order:
                        seq = sequences.get(s)
                        if not seq or pos >= len(seq):
                            vcf.write("\t./.")
                            continue
                        base = seq[pos].upper()
                        gt = encode_gt_for_sample_invariant(base, ref, warned_bdhv)
                        vcf.write(f"\t{gt}")
                    vcf.write("\n")
                    emitted_records += 1
                    emitted_invariant += 1
                    continue

                tmpl = template_variants[key]
                ref = tmpl["ref"]
                alt = list(tmpl["alt"])
                alt_field = ",".join(alt) if alt else "."
                info = f"POLYMORPHIC;NS={summary['ns_valid']}"

                vcf.write(f"{chrom}\t{vcf_pos}\t{var_id}\t{ref}\t{alt_field}\t.\tPASS\t{info}\tGT")

                tmpl_gts = tmpl.get("gts", {})
                for s in samples_order:
                    vcf.write(f"\t{tmpl_gts.get(s, './.')}")
                vcf.write("\n")
                emitted_records += 1
                emitted_polymorphic += 1

    return {
        "emitted_records": emitted_records,
        "emitted_invariant": emitted_invariant,
        "emitted_polymorphic": emitted_polymorphic,
        "skipped_no_summary": skipped_no_summary,
        "template_variants_total": len(template_variants),
        "template_variants_seen_in_loci": emitted_polymorphic,
        "template_variants_missing_from_loci": 0,
        "loci_polymorphic_not_in_template": 0,
        "template_only_rows_written": 0,
        "loci_after_filter": n_loci,
        "samples_after_filter": len(samples_order),
    }

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Walk ipyrad .loci: emit INVARIANT rows from the alignment, POLYMORPHIC rows from --template-vcf.\n"
            "Variable sites in .loci and variant rows in the VCF must match exactly (same CHROM/POS set)."
        )
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
    parser.add_argument(
        "--template-vcf",
        required=True,
        help="ipyrad variant VCF(.gz). Set of (CHROM,POS) must equal variable columns in .loci; REF/ALT/GT copied per variable site.",
    )

    args = parser.parse_args()

    # Read samples file if provided
    samples_filter = None
    if args.samples_file:
        samples_filter = read_samples_file(args.samples_file)
        if not samples_filter:
            print(f"ERROR: No samples found in filter file: {args.samples_file}", file=sys.stderr)
            sys.exit(1)

    print(f"Reading loci file (pass 1/2): {args.input_file}", flush=True)
    loci_var, all_samples_in_loci, n_loci = collect_loci_pass1(args.input_file)
    print(f"Found {n_loci} loci", flush=True)

    if n_loci == 0:
        print(f"ERROR: No loci found in input file: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    if not all_samples_in_loci:
        print(f"ERROR: No samples found in loci file: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    if samples_filter:
        samples_in_filter_not_in_loci = samples_filter - all_samples_in_loci
        if samples_in_filter_not_in_loci:
            print(f"ERROR: {len(samples_in_filter_not_in_loci)} sample(s) in filter file not found in loci file:", file=sys.stderr)
            for sample in sorted(samples_in_filter_not_in_loci):
                print(f"  {sample}", file=sys.stderr)
            sys.exit(1)

        filtered_samples = all_samples_in_loci.intersection(samples_filter)
        if not filtered_samples:
            print("ERROR: No samples match between filter file and loci file", file=sys.stderr)
            sys.exit(1)

    print(f"Reading template VCF: {args.template_vcf}", flush=True)
    template_variants = parse_template_vcf(args.template_vcf)
    _validate_loci_template_variable_sets(loci_var, template_variants)

    if samples_filter is not None:
        samples_order = sorted(all_samples_in_loci.intersection(samples_filter))
    else:
        samples_order = sorted(all_samples_in_loci)

    print("Writing VCF (pass 2/2)...", flush=True)
    write_loci_pass2(
        args.input_file,
        args.output,
        template_variants,
        samples_order,
        n_loci,
    )
    print(f"VCF written: {args.output}")

if __name__ == "__main__":
    main()
