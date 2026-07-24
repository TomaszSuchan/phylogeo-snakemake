#!/usr/bin/env python3
import math
import sys
import traceback
from collections import Counter
from functools import lru_cache

import numpy as np
import pandas as pd

# Redirect stdout and stderr to log file first so any error is visible in the rule log
log_file = snakemake.log[0]
log_handle = open(log_file, "w")
sys.stdout = log_handle
sys.stderr = log_handle

try:
    bootstrap_replicates = snakemake.params.bootstrap_replicates
except Exception:
    traceback.print_exc()
    raise


def weighted_mean(values, weights):
    return (values * weights).sum() / weights.sum()


def bootstrap_ci(values, weights, B=1000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(values)
    boot = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        boot[b] = weighted_mean(values[idx], weights[idx])
    return np.percentile(boot, [2.5, 97.5]), boot.std(ddof=1)


def bootstrap_ratio_ci(numerators, denominators, B=1000, seed=42):
    """Bootstrap CI for sum(numerators) / sum(denominators) over resampled windows."""
    rng = np.random.default_rng(seed)
    n = len(numerators)
    boot = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        denom = denominators[idx].sum()
        boot[b] = np.nan if denom == 0 else numerators[idx].sum() / denom
    valid = boot[~np.isnan(boot)]
    if len(valid) == 0:
        return (np.nan, np.nan), np.nan
    return np.percentile(valid, [2.5, 97.5]), float(np.std(valid, ddof=1)) if len(valid) > 1 else 0.0


@lru_cache(maxsize=1024)
def _tajima_constants(n):
    """Return (e1, e2) Tajima 1989 coefficients (mirrors pixy.calc._tajima_constants)."""
    if n < 2:
        return 0.0, 0.0
    i = np.arange(1, n)
    a1 = float(np.sum(1.0 / i))
    a2 = float(np.sum(1.0 / (i ** 2)))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 ** 2))
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)
    return e1, e2


def parse_tajima_d_s_counts(value):
    """Parse pixy tajima_d_s_counts ('n:s,n:s,...') into a Counter."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return Counter()
    text = str(value).strip()
    if text in ("", "NA", "nan"):
        return Counter()
    counts = Counter()
    for item in text.split(","):
        n_str, s_str = item.split(":")
        counts[int(n_str)] += int(float(s_str))
    return counts


def combine_tajima_d_s_counts(values):
    total = Counter()
    for value in values:
        total.update(parse_tajima_d_s_counts(value))
    return total


def calc_tajima_d_stdev(s_counts):
    """Denominator of Tajima's D from summed observed-allele-count classes (Bailey et al. 2025)."""
    d_stdev = 0.0
    for n, s in s_counts.items():
        if n < 2 or s <= 0:
            continue
        e1, e2 = _tajima_constants(int(n))
        variance = (e1 * s) + (e2 * s * (s - 1))
        if variance < 0:
            return float("nan")
        d_stdev += math.sqrt(variance)
    return float(d_stdev)


def aggregate_tajima_d(raw_pi, raw_watterson_theta, s_count_values):
    raw_pi_sum = float(np.nansum(raw_pi))
    raw_theta_sum = float(np.nansum(raw_watterson_theta))
    d_stdev = calc_tajima_d_stdev(combine_tajima_d_s_counts(s_count_values))
    if not np.isfinite(d_stdev) or d_stdev <= 0:
        tajima_d = np.nan
    else:
        tajima_d = (raw_pi_sum - raw_theta_sum) / d_stdev
    return tajima_d, d_stdev, raw_pi_sum, raw_theta_sum


def bootstrap_tajima_d_ci(raw_pi, raw_watterson_theta, s_count_values, B=1000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(raw_pi)
    boot = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        boot[b], _, _, _ = aggregate_tajima_d(
            raw_pi[idx],
            raw_watterson_theta[idx],
            [s_count_values[i] for i in idx],
        )
    valid = boot[~np.isnan(boot)]
    if len(valid) == 0:
        return (np.nan, np.nan), np.nan
    return np.percentile(valid, [2.5, 97.5]), float(np.std(valid, ddof=1)) if len(valid) > 1 else 0.0


def process_pi(file, output, bootstrap_replicates=1000):
    df = pd.read_csv(file, sep="\t")
    df_clean = df[df["avg_pi"].notna()].copy()

    results = []
    for pop, group_df in df_clean.groupby("pop"):
        values = group_df["avg_pi"].to_numpy()
        weights = group_df["no_sites"].to_numpy()

        m = weighted_mean(values, weights)
        (ci_low, ci_high), se = bootstrap_ci(values, weights, B=bootstrap_replicates)

        results.append({
            "population": pop,
            "mean_pi": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_loci": len(group_df),
            "n_sites_total": int(weights.sum()),
        })

    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Pi results written to {output}")


def process_fst(file, output, bootstrap_replicates=1000):
    df = pd.read_csv(file, sep="\t")
    df_clean = df[df["avg_wc_fst"].notna()].copy()

    results = []
    for (pop1, pop2), group_df in df_clean.groupby(["pop1", "pop2"]):
        values = group_df["avg_wc_fst"].to_numpy()
        weights = group_df["no_snps"].to_numpy()

        m = weighted_mean(values, weights)
        (ci_low, ci_high), se = bootstrap_ci(values, weights, B=bootstrap_replicates)

        results.append({
            "pop1": pop1,
            "pop2": pop2,
            "mean_fst": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_loci": len(group_df),
            "n_snps_total": int(weights.sum()),
        })

    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Fst results written to {output}")


def process_dxy(file, output, bootstrap_replicates=1000):
    df = pd.read_csv(file, sep="\t")
    df_clean = df[df["avg_dxy"].notna()].copy()

    results = []
    for (pop1, pop2), group_df in df_clean.groupby(["pop1", "pop2"]):
        values = group_df["avg_dxy"].to_numpy()
        weights = group_df["no_sites"].to_numpy()

        m = weighted_mean(values, weights)
        (ci_low, ci_high), se = bootstrap_ci(values, weights, B=bootstrap_replicates)

        results.append({
            "pop1": pop1,
            "pop2": pop2,
            "mean_dxy": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_loci": len(group_df),
            "n_sites_total": int(weights.sum()),
        })

    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Dxy results written to {output}")


def process_watterson_theta(file, output, bootstrap_replicates=1000):
    df = pd.read_csv(file, sep="\t")
    required = {"avg_watterson_theta", "raw_watterson_theta", "no_sites"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"watterson_theta output missing columns {sorted(missing)}; "
            "need pixy >= 2.0.0 with --stats watterson_theta"
        )
    df_clean = df[df["avg_watterson_theta"].notna()].copy()

    results = []
    for pop, group_df in df_clean.groupby("pop"):
        raw = group_df["raw_watterson_theta"].to_numpy(dtype=float)
        sites = group_df["no_sites"].to_numpy(dtype=float)
        site_sum = float(sites.sum())
        m = np.nan if site_sum == 0 else float(raw.sum() / site_sum)
        (ci_low, ci_high), se = bootstrap_ratio_ci(raw, sites, B=bootstrap_replicates)

        results.append({
            "population": pop,
            "mean_watterson_theta": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_loci": len(group_df),
            "n_sites_total": int(site_sum),
        })

    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Watterson theta results written to {output}")


def process_tajima_d(file, output, bootstrap_replicates=1000):
    df = pd.read_csv(file, sep="\t")
    required = {"tajima_d", "raw_pi", "raw_watterson_theta", "no_sites", "tajima_d_s_counts"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"tajima_d output missing columns {sorted(missing)}; "
            "run pixy with --stats tajima_d --tajima_components"
        )
    df_clean = df[df["tajima_d"].notna()].copy()

    results = []
    for pop, group_df in df_clean.groupby("pop"):
        raw_pi = group_df["raw_pi"].to_numpy(dtype=float)
        raw_theta = group_df["raw_watterson_theta"].to_numpy(dtype=float)
        s_counts = group_df["tajima_d_s_counts"].tolist()
        sites = group_df["no_sites"].to_numpy(dtype=float)

        m, d_stdev, raw_pi_sum, raw_theta_sum = aggregate_tajima_d(raw_pi, raw_theta, s_counts)
        (ci_low, ci_high), se = bootstrap_tajima_d_ci(
            raw_pi, raw_theta, s_counts, B=bootstrap_replicates
        )

        results.append({
            "population": pop,
            "mean_tajima_d": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "tajima_d_stdev": d_stdev,
            "raw_pi_total": raw_pi_sum,
            "raw_watterson_theta_total": raw_theta_sum,
            "n_loci": len(group_df),
            "n_sites_total": int(sites.sum()),
        })

    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Tajima D results written to {output}")


# Which stat to process is set by the rule (params.stat)
try:
    stat = getattr(snakemake.params, "stat", None)
    if stat is None:
        # Fallback: infer from which input key is present
        for key in ("pi", "fst", "dxy", "watterson_theta", "tajima_d"):
            if hasattr(snakemake.input, key):
                stat = key
                break
    if stat is None:
        raise ValueError(
            "params.stat is required; add params: "
            "stat = 'pi'|'fst'|'dxy'|'watterson_theta'|'tajima_d' to the rule"
        )
    infile = getattr(snakemake.input, stat)
    outfile = getattr(snakemake.output, stat)

    if stat == "pi":
        process_pi(infile, outfile, bootstrap_replicates)
    elif stat == "fst":
        process_fst(infile, outfile, bootstrap_replicates)
    elif stat == "dxy":
        process_dxy(infile, outfile, bootstrap_replicates)
    elif stat == "watterson_theta":
        process_watterson_theta(infile, outfile, bootstrap_replicates)
    elif stat == "tajima_d":
        process_tajima_d(infile, outfile, bootstrap_replicates)
    else:
        raise ValueError(f"Unknown stat: {stat}")
except Exception:
    traceback.print_exc()
    raise
