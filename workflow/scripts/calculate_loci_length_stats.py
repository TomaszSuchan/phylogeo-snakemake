#!/usr/bin/env python3
"""
RAD locus length distributions from an ipyrad .loci file.

Within a locus, ipyrad typically pads sequences to a shared alignment width, but
the number of called bases (non-N / non-gap) still varies among individuals.
This script therefore reports:

  - alignment width per locus (len of the padded sequences)
  - effective length per sample×locus (count of bases not in - / N / .)
  - per-locus summaries of effective length (min / median / max across samples
    present at that locus) and the within-locus effective-length range

Optional --samples-file restricts to retained samples (same idea as
loci_to_alignment / pixy).
"""

import os
import sys
from pathlib import Path

try:
    _SCRIPT_DIR = Path(snakemake.scriptdir)
except NameError:
    _SCRIPT_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(_SCRIPT_DIR / "extract-invariant-vcf"))
from extract_invariant_vcf import iter_loci_file  # noqa: E402

_MISSING = frozenset("-N.")


def _read_keep_samples(path):
    samples = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            sample = line.strip()
            if sample:
                samples.add(sample)
    if not samples:
        raise SystemExit(f"No sample names found in {path}")
    return samples


def _effective_length(seq):
    return sum(1 for c in seq if c not in _MISSING)


def _percentile(sorted_vals, p):
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    rank = (len(sorted_vals) - 1) * (p / 100.0)
    lo = int(rank)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = rank - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def _fmt(x, digits=2):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.{digits}f}"
    return str(x)


def _summarize(values):
    if not values:
        return None
    vals = sorted(values)
    n = len(vals)
    return {
        "n": n,
        "min": vals[0],
        "q1": _percentile(vals, 25),
        "median": _percentile(vals, 50),
        "mean": sum(vals) / n,
        "q3": _percentile(vals, 75),
        "max": vals[-1],
    }


def _write_summary(out, title, summary, unit="bp"):
    out.write(f"{title}\n")
    if summary is None:
        out.write("  (no data)\n")
        return
    out.write(f"  n: {summary['n']}\n")
    out.write(
        f"  min / Q1 / median / mean / Q3 / max ({unit}): "
        f"{_fmt(summary['min'], 0)} / {_fmt(summary['q1'])} / "
        f"{_fmt(summary['median'])} / {_fmt(summary['mean'])} / "
        f"{_fmt(summary['q3'])} / {_fmt(summary['max'], 0)}\n"
    )


def _write_bins(out, values, bins):
    if not values:
        return
    n = len(values)
    for label, lo, hi in bins:
        if hi is None:
            count = sum(1 for x in values if x >= lo)
        else:
            count = sum(1 for x in values if lo <= x < hi)
        out.write(f"  {label}: {count} ({100.0 * count / n:.2f}%)\n")


_LENGTH_BINS = (
    ("short (<200 bp)", 0, 200),
    ("medium (200-399 bp)", 200, 400),
    ("long (>=400 bp)", 400, None),
)


def calculate_loci_length_stats(loci_path, output_file, samples_file=None):
    keep = _read_keep_samples(samples_file) if samples_file else None

    alignment_widths = []
    effective_obs = []  # every sample×locus effective length
    locus_eff_min = []
    locus_eff_median = []
    locus_eff_max = []
    locus_eff_range = []
    n_loci_align_heterogeneous = 0
    n_loci_eff_heterogeneous = 0
    n_loci_skipped_no_keep = 0
    n_samples_seen = set()
    n_ref_span = 0
    ref_spans = []

    for locus in iter_loci_file(loci_path):
        sequences = locus["sequences"]
        if keep is not None:
            present = [s for s in sequences if s in keep]
            if not present:
                n_loci_skipped_no_keep += 1
                continue
        else:
            present = list(sequences)

        lengths = [len(sequences[s]) for s in present]
        effs = [_effective_length(sequences[s]) for s in present]
        for s in present:
            n_samples_seen.add(s)

        if len(set(lengths)) > 1:
            n_loci_align_heterogeneous += 1
        if len(set(effs)) > 1:
            n_loci_eff_heterogeneous += 1

        width = max(lengths)
        alignment_widths.append(width)
        effective_obs.extend(effs)

        eff_sorted = sorted(effs)
        locus_eff_min.append(eff_sorted[0])
        locus_eff_median.append(_percentile(eff_sorted, 50))
        locus_eff_max.append(eff_sorted[-1])
        locus_eff_range.append(eff_sorted[-1] - eff_sorted[0])

        start, end = locus.get("start"), locus.get("end")
        if start is not None and end is not None and end >= start:
            n_ref_span += 1
            ref_spans.append(end - start + 1)

    if not alignment_widths:
        raise SystemExit(
            f"No loci retained from {loci_path}"
            + (f" after filtering to {samples_file}" if samples_file else "")
        )

    loci_abspath = os.path.abspath(loci_path)
    path_parts = loci_abspath.split(os.sep)
    if "results" in path_parts:
        results_idx = path_parts.index("results")
        project_root = os.sep.join(path_parts[:results_idx])
        loci_relpath = os.path.relpath(loci_abspath, project_root)
    else:
        loci_relpath = loci_abspath

    n_loci = len(alignment_widths)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as out:
        out.write(f"Loci file name: {os.path.basename(loci_path)}\n")
        out.write(f"Loci file location: {loci_relpath}\n")
        if samples_file:
            out.write(f"Samples file: {samples_file}\n")
            out.write(f"Samples retained (in keep list and seen): {len(n_samples_seen)}\n")
        else:
            out.write(f"Samples seen: {len(n_samples_seen)}\n")
        out.write(f"Number of loci: {n_loci}\n")
        if keep is not None:
            out.write(
                f"Loci skipped (no keep-list sample present): {n_loci_skipped_no_keep}\n"
            )
        out.write(
            f"Loci with heterogeneous alignment width across samples: "
            f"{n_loci_align_heterogeneous} "
            f"({100.0 * n_loci_align_heterogeneous / n_loci:.2f}%)\n"
        )
        out.write(
            f"Loci with heterogeneous effective length across samples: "
            f"{n_loci_eff_heterogeneous} "
            f"({100.0 * n_loci_eff_heterogeneous / n_loci:.2f}%)\n"
        )
        out.write("\n")
        out.write("How lengths are defined\n")
        out.write(
            "  Each locus block in the .loci file lists one sequence per sample that "
            "has data at that locus (after optional keep-list filtering). ipyrad pads "
            "those sequences so that within a locus they usually share the same string "
            "length, with flanking or internal missing sequence coded as N (and gaps as -).\n"
        )
        out.write(
            "  Alignment width for a locus = len(sequence) of that padded string "
            "(max across samples if widths ever differ). One value per locus.\n"
        )
        out.write(
            "  Effective length for one sample at one locus = number of characters in "
            "that sample's sequence that are neither gap (-) nor missing (N or .). "
            "IUPAC ambiguity codes other than N count as present. Because samples differ "
            "in how much of the locus is called vs N-padded, effective length usually "
            "varies among individuals at the same locus even when alignment width does not.\n"
        )
        out.write(
            "  For each locus we therefore have many effective lengths (one per sample "
            "present). We summarise those within the locus (min / median / max / range), "
            "then summarise that one-number-per-locus list across all loci. The "
            "min/Q1/median/mean/Q3/max lines below are always that outer, "
            "across-loci (or across-observations) distribution.\n"
        )
        out.write(
            "  Length bins: short <200 bp, medium 200-399 bp, long >=400 bp "
            "(counts of values falling in each bin).\n"
        )
        out.write("\n")

        out.write("Alignment width (per locus)\n")
        out.write(
            "  Calculation: for each locus, take the padded sequence length "
            "(usually identical for all samples at that locus). The summary below "
            "is the distribution of those one-per-locus widths across loci.\n"
        )
        _write_summary(out, "  Summary", _summarize(alignment_widths))
        _write_bins(out, alignment_widths, _LENGTH_BINS)
        out.write("\n")

        out.write("Effective length (per sample × locus observation)\n")
        out.write(
            "  Calculation: for every retained sample at every locus, compute that "
            "sample's effective length (non-gap / non-N bases). Pool all of those "
            "sample×locus values into one list (so a locus with 200 samples "
            "contributes 200 numbers). The summary below is the distribution of "
            "that pooled list. This is the most direct view of how much sequence "
            "individuals actually contribute, but loci with more samples weigh "
            "more heavily than loci with fewer samples.\n"
        )
        _write_summary(out, "  Summary", _summarize(effective_obs))
        _write_bins(out, effective_obs, _LENGTH_BINS)
        out.write("\n")

        out.write("Per-locus median effective length\n")
        out.write(
            "  Calculation: for each locus separately, collect the effective lengths "
            "of all samples present at that locus and take their median (one number "
            "per locus — a typical individual's called length at that locus). Then "
            "summarise that list of per-locus medians across loci. Unlike the pooled "
            "sample×locus section above, every locus contributes exactly once, "
            "regardless of how many samples it has.\n"
        )
        _write_summary(out, "  Summary", _summarize(locus_eff_median))
        _write_bins(out, locus_eff_median, _LENGTH_BINS)
        out.write("\n")

        out.write("Per-locus minimum effective length\n")
        out.write(
            "  Calculation: for each locus, take the smallest effective length among "
            "samples present (the sample with the least called sequence at that "
            "locus). Then summarise those per-locus minima across loci.\n"
        )
        _write_summary(out, "  Summary", _summarize(locus_eff_min))
        out.write("\n")

        out.write("Per-locus maximum effective length\n")
        out.write(
            "  Calculation: for each locus, take the largest effective length among "
            "samples present (the sample with the most called sequence at that "
            "locus). Then summarise those per-locus maxima across loci. This is often "
            "closer to the alignment width than the median, but can still be shorter "
            "when even the best-covered sample has internal or terminal Ns.\n"
        )
        _write_summary(out, "  Summary", _summarize(locus_eff_max))
        out.write("\n")

        out.write("Within-locus effective-length range (max − min across samples)\n")
        out.write(
            "  Calculation: for each locus, subtract the minimum sample effective "
            "length from the maximum sample effective length at that same locus. "
            "That range is how much individuals differ in called length within the "
            "locus. Then summarise those per-locus ranges across loci. A range of 0 "
            "means every sample at that locus had the same effective length.\n"
        )
        _write_summary(out, "  Summary", _summarize(locus_eff_range))

        if ref_spans:
            out.write("\n")
            out.write("Reference span (per locus)\n")
            out.write(
                f"  Calculation: for reference-mapped loci, the closing // line of "
                f"the .loci block carries |N:CHROM:START-END|. Reference span = "
                f"END − START + 1 (one value per mapped locus; "
                f"{n_ref_span}/{n_loci} loci had this trailer). This is the genomic "
                f"interval the locus is mapped to, not the amount of sequence called "
                f"in any sample. The summary below is the distribution of those spans "
                f"across mapped loci.\n"
            )
            _write_summary(out, "  Summary", _summarize(ref_spans))
    print(
        f"Wrote locus length stats for {n_loci} loci "
        f"({len(n_samples_seen)} samples) to {output_file}",
        flush=True,
    )


if __name__ == "__main__":
    try:
        loci_path = snakemake.input.loci
        samples_file = snakemake.input.samples
        output_file = snakemake.output.stats
    except NameError:
        import argparse

        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("loci", help="ipyrad .loci file")
        parser.add_argument("output", help="output stats text file")
        parser.add_argument(
            "--samples-file",
            default=None,
            help="optional keep list (one sample per line)",
        )
        args = parser.parse_args()
        loci_path = args.loci
        samples_file = args.samples_file
        output_file = args.output

    calculate_loci_length_stats(loci_path, output_file, samples_file=samples_file)
