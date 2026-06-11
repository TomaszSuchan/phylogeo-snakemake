#!/usr/bin/env python3
"""
Diagnose ipyrad .loci vs template SNP VCF for extract_invariant_vcf-style parsing.

Checks (current extract_invariant_vcf.py behavior):
  - Strict layout: sample rows, then one closing "//" line per locus; coordinates only
    from |N:CHROM:START-END| on that closing line (no separate "|" metadata lines).
  - For reference-mapped loci: every alignment column i maps to POS = START + i and
    must satisfy START <= POS <= END when columns are capped by span.
  - Optional: stream template VCF (CHROM, POS) keys and report sites classified as
    variable by summarize_site() in .loci but missing from the VCF (IUPAC pattern).

Streams VCF (CHROM, POS) keys only; does not load per-site genotypes.
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import sys
from collections import Counter
from pathlib import Path


def _load_extract_module():
    here = Path(__file__).resolve().parent
    mod_path = here / "extract_invariant_vcf.py"
    spec = importlib.util.spec_from_file_location("extract_invariant_vcf", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def stream_vcf_keys(vcf_path: str) -> set[tuple[str, int]]:
    keys: set[tuple[str, int]] = set()
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 2)
            if len(parts) < 2:
                continue
            keys.add((parts[0], int(parts[1])))
    return keys


def closing_trailer_snippet(line: str, width: int = 90) -> str:
    s = line.rstrip("\r\n")
    if len(s) <= width:
        return s
    return "…" + s[-width:]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loci", required=True, help="Path to ipyrad .loci file")
    parser.add_argument(
        "--vcf",
        default=None,
        help="Optional template VCF (.vcf or .vcf.gz); used for key set and mismatch scan",
    )
    parser.add_argument(
        "--max-loci",
        type=int,
        default=None,
        metavar="N",
        help="Stop after N loci (default: all)",
    )
    parser.add_argument(
        "--mismatch-examples",
        type=int,
        default=5,
        metavar="K",
        help="Max examples of variable-in-loci / missing-in-VCF to print (default 5)",
    )
    args = parser.parse_args()

    eiv = _load_extract_module()

    tmpl: set[tuple[str, int]] | None = None
    if args.vcf:
        print(f"Streaming VCF keys: {args.vcf}", file=sys.stderr)
        tmpl = stream_vcf_keys(args.vcf)
        print(f"  variant rows: {len(tmpl)}", file=sys.stderr)

    n_loci = 0
    n_refmapped = 0
    n_denovo = 0
    col_oob = 0  # POS > END or POS < START for any column in iteration range
    mismatch_examples: list[tuple[str, int, int, str, Counter[str]]] = []
    mismatch_total = 0

    # Light pass: show trailing part of first few closing "//" lines (coordinates trailer).
    raw_closing_snippets: list[str] = []
    with open(args.loci, "r", encoding="utf-8", errors="replace") as raw_f:
        for raw in raw_f:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("//") and "|" in line and ":" in line:
                raw_closing_snippets.append(closing_trailer_snippet(line))
                if len(raw_closing_snippets) >= 3:
                    break

    # Use the same iterator as production (validates strict layout)
    for locus_idx, locus in enumerate(eiv.iter_loci_file(args.loci)):
        if args.max_loci is not None and locus_idx >= args.max_loci:
            break
        n_loci += 1
        chrom = locus.get("chrom")
        start = locus.get("start")
        end = locus.get("end")
        seqs = locus["sequences"]

        if chrom and start is not None and end is not None:
            n_refmapped += 1
            span = end - start + 1
            max_len = max((len(s) for s in seqs.values()), default=0)
            cap = min(max_len, span)
            for pos in range(cap):
                pos_vcf = eiv.vcf_position(start, pos)
                if pos_vcf < start or pos_vcf > end:
                    col_oob += 1
                if tmpl is None:
                    continue
                sm = eiv.summarize_site(seqs, pos)
                if sm is None or sm["is_invariant"]:
                    continue
                key = (chrom, pos_vcf)
                if key not in tmpl:
                    mismatch_total += 1
                    if len(mismatch_examples) < args.mismatch_examples:
                        bases = [seq[pos].upper() for seq in seqs.values() if pos < len(seq)]
                        mismatch_examples.append(
                            (chrom, pos_vcf, locus_idx, f"col={pos}", Counter(bases))
                        )
        else:
            n_denovo += 1

    print("=== diagnose_loci_template_vcf ===")
    print(f"loci_file: {args.loci}")
    print(f"loci_parsed: {n_loci} (ref_mapped={n_refmapped}, de_novo_or_no_coords={n_denovo})")
    print(f"column_POS_outside_[START,END]_in_scanned_range: {col_oob}")
    if args.vcf:
        print(f"template_vcf: {args.vcf}")
        print(
            f"variable_in_loci_missing_in_vcf: total={mismatch_total} "
            f"(showing up to {args.mismatch_examples} examples)"
        )
        for chrom, pos_vcf, li, col_label, cnt in mismatch_examples:
            top = ", ".join(f"{b}={n}" for b, n in cnt.most_common(8))
            print(f"  {chrom}:{pos_vcf}  locus_idx={li}  {col_label}  [{top}]")

    print("\n--- closing // line trailers (first up to 3 ref-style lines seen) ---")
    if raw_closing_snippets:
        for i, snip in enumerate(raw_closing_snippets, 1):
            print(f"  [{i}] {snip}")
    else:
        print("  (none captured; file may be empty or format differs)")

    print(
        "\nNote: coordinates are taken only from the closing \"//\" line after each "
        "alignment block (see extract_invariant_vcf.iter_loci_file). "
        "Snpstring * / - on that line are not used for POS.",
    )


if __name__ == "__main__":
    main()
