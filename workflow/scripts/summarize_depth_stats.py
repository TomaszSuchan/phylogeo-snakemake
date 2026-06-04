#!/usr/bin/env python3
"""Summarise VCFtools sequencing-depth outputs."""

import math
import statistics
import sys


DEPTH_THRESHOLDS = (3, 5, 10)
HIGH_DEPTH_THRESHOLD = 100


def read_table(path):
    """Read a whitespace-delimited VCFtools table into a list of dictionaries."""
    with open(path) as handle:
        header = None
        rows = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                continue
            rows.append(dict(zip(header, fields)))
    return rows


def as_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def mean(values):
    return statistics.mean(values) if values else float("nan")


def median(values):
    return statistics.median(values) if values else float("nan")


def stdev(values):
    return statistics.stdev(values) if len(values) > 1 else float("nan")


def fmt(value):
    return "NA" if value is None or not math.isfinite(value) else f"{value:.6f}"


def summarize(individual_depth_path, site_depth_path, output_path):
    ind_rows = read_table(individual_depth_path)
    site_rows = read_table(site_depth_path)

    individual_depths = []
    individual_called_sites = []
    for row in ind_rows:
        depth = as_float(row.get("MEAN_DEPTH"))
        if depth is not None:
            individual_depths.append(depth)
            individual_called_sites.append(as_int(row.get("N_SITES")))

    site_depths = []
    for row in site_rows:
        depth = as_float(row.get("MEAN_DEPTH"))
        if depth is not None:
            site_depths.append(depth)

    total_called_genotypes = sum(individual_called_sites)
    if total_called_genotypes:
        weighted_depth_sum = sum(
            depth * n_sites
            for depth, n_sites in zip(individual_depths, individual_called_sites)
        )
        overall_mean_depth = weighted_depth_sum / total_called_genotypes
    else:
        overall_mean_depth = float("nan")

    with open(output_path, "w") as out:
        out.write("Sequencing depth summary\n")
        out.write(f"Individual depth file: {individual_depth_path}\n")
        out.write(f"Site depth file: {site_depth_path}\n")
        out.write(f"Number of individuals: {len(individual_depths)}\n")
        out.write(f"Number of sites: {len(site_depths)}\n")
        out.write(f"Called genotypes with depth: {total_called_genotypes}\n")
        out.write(f"Overall mean depth across called genotypes: {fmt(overall_mean_depth)}\n")
        out.write("\n")

        out.write("Per-individual mean depth\n")
        out.write(f"Mean: {fmt(mean(individual_depths))}\n")
        out.write(f"Median: {fmt(median(individual_depths))}\n")
        out.write(f"SD: {fmt(stdev(individual_depths))}\n")
        out.write(f"Minimum: {fmt(min(individual_depths) if individual_depths else float('nan'))}\n")
        out.write(f"Maximum: {fmt(max(individual_depths) if individual_depths else float('nan'))}\n")
        for threshold in DEPTH_THRESHOLDS:
            count = sum(depth < threshold for depth in individual_depths)
            pct = 100 * count / len(individual_depths) if individual_depths else float("nan")
            out.write(
                f"Individuals with mean depth < {threshold}: {count} ({fmt(pct)}%)\n"
            )
        count = sum(depth > HIGH_DEPTH_THRESHOLD for depth in individual_depths)
        pct = 100 * count / len(individual_depths) if individual_depths else float("nan")
        out.write(
            f"Individuals with mean depth > {HIGH_DEPTH_THRESHOLD}: {count} ({fmt(pct)}%)\n"
        )
        out.write("\n")

        out.write("Per-site mean depth\n")
        out.write(f"Mean: {fmt(mean(site_depths))}\n")
        out.write(f"Median: {fmt(median(site_depths))}\n")
        out.write(f"SD: {fmt(stdev(site_depths))}\n")
        out.write(f"Minimum: {fmt(min(site_depths) if site_depths else float('nan'))}\n")
        out.write(f"Maximum: {fmt(max(site_depths) if site_depths else float('nan'))}\n")
        for threshold in DEPTH_THRESHOLDS:
            count = sum(depth < threshold for depth in site_depths)
            pct = 100 * count / len(site_depths) if site_depths else float("nan")
            out.write(f"Sites with mean depth < {threshold}: {count} ({fmt(pct)}%)\n")
        count = sum(depth > HIGH_DEPTH_THRESHOLD for depth in site_depths)
        pct = 100 * count / len(site_depths) if site_depths else float("nan")
        out.write(f"Sites with mean depth > {HIGH_DEPTH_THRESHOLD}: {count} ({fmt(pct)}%)\n")


if __name__ == "__main__":
    try:
        snakemake_obj = snakemake  # type: ignore[name-defined]
        individual_depth = snakemake_obj.input.individual_depth
        site_depth = snakemake_obj.input.site_depth
        output = snakemake_obj.output.summary
    except NameError:
        if len(sys.argv) != 4:
            print(
                "Usage: summarize_depth_stats.py <input.idepth> <input.ldepth.mean> <output.txt>",
                file=sys.stderr,
            )
            sys.exit(1)
        individual_depth, site_depth, output = sys.argv[1:]

    summarize(individual_depth, site_depth, output)
