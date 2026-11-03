#!/usr/bin/env python3
"""
Calculate VCF statistics:
- Number of chromosomes (unique values in column 1)
- Number of RAD fragments (unique first part of column 3, before "_pos")
- Number of variants (total rows)
- Number of samples (from VCF header)
- SNPs-per-RAD-locus distribution
"""

import sys
import gzip
import os
from collections import defaultdict


def extract_rad_fragment(loc_id):
    """Extract RAD fragment ID from location ID (first part before '_pos')."""
    parts = loc_id.split("_pos")
    return parts[0] if parts else loc_id


def _percentile(sorted_vals, p):
    """Linear-interpolated percentile for a non-empty sorted list."""
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
    """Return dict of count / min / quartiles / mean / max for a numeric list."""
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


def _write_summary(out, title, summary):
    out.write(f"{title}\n")
    if summary is None:
        out.write("  (no data)\n")
        return
    out.write(f"  n: {summary['n']}\n")
    out.write(
        f"  min / Q1 / median / mean / Q3 / max: "
        f"{_fmt(summary['min'], 0)} / {_fmt(summary['q1'])} / "
        f"{_fmt(summary['median'])} / {_fmt(summary['mean'])} / "
        f"{_fmt(summary['q3'])} / {_fmt(summary['max'], 0)}\n"
    )


def calculate_stats_vcf(vcf_file, output_file):
    """Calculate and write VCF statistics to output file."""
    chromosomes = set()
    snps_per_locus = defaultdict(int)
    variant_count = 0
    num_samples = 0

    if vcf_file.endswith(".gz"):
        f = gzip.open(vcf_file, "rt")
    else:
        f = open(vcf_file, "r")

    try:
        with f:
            for line in f:
                if line.startswith("#CHROM"):
                    fields = line.strip().split("\t")
                    if len(fields) > 9:
                        num_samples = len(fields) - 9
                    continue

                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 3:
                    continue

                chromosomes.add(fields[0])
                snps_per_locus[extract_rad_fragment(fields[2])] += 1
                variant_count += 1
    except Exception as e:
        print(f"Error reading VCF file: {e}", file=sys.stderr)
        sys.exit(1)

    snp_count_summary = _summarize(list(snps_per_locus.values()))

    vcf_filename = os.path.basename(vcf_file)
    vcf_abspath = os.path.abspath(vcf_file)

    path_parts = vcf_abspath.split(os.sep)
    if "results" in path_parts:
        results_idx = path_parts.index("results")
        project_root = os.sep.join(path_parts[:results_idx])
        vcf_relpath = os.path.relpath(vcf_abspath, project_root)
    else:
        vcf_relpath = vcf_abspath

    n_rad = len(snps_per_locus)
    try:
        with open(output_file, "w") as out:
            out.write(f"VCF file name: {vcf_filename}\n")
            out.write(f"VCF file location: {vcf_relpath}\n")
            out.write(f"Number of chromosomes: {len(chromosomes)}\n")
            out.write(f"Number of RAD fragments: {n_rad}\n")
            out.write(f"Number of variants: {variant_count}\n")
            out.write(f"Number of samples: {num_samples}\n")
            out.write("\n")
            _write_summary(out, "SNPs per RAD locus", snp_count_summary)
            if snp_count_summary is not None:
                n_singleton = sum(1 for c in snps_per_locus.values() if c == 1)
                pct = 100.0 * n_singleton / snp_count_summary["n"]
                out.write(
                    f"  loci with exactly 1 SNP: {n_singleton} ({pct:.2f}%)\n"
                )
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Statistics calculated:")
    print(f"  VCF file: {vcf_filename}")
    print(f"  VCF location: {vcf_relpath}")
    print(f"  Chromosomes: {len(chromosomes)}")
    print(f"  RAD fragments: {n_rad}")
    print(f"  Variants: {variant_count}")
    print(f"  Samples: {num_samples}")


if __name__ == "__main__":
    try:
        vcf_file = snakemake.input.vcf
        output_file = snakemake.output.stats
    except NameError:
        if len(sys.argv) != 3:
            print(
                "Usage: calculate_vcf_stats.py <input.vcf.gz> <output.txt>",
                file=sys.stderr,
            )
            sys.exit(1)
        vcf_file = sys.argv[1]
        output_file = sys.argv[2]

    calculate_stats_vcf(vcf_file, output_file)
