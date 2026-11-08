#!/usr/bin/env python3
import math
import sys
import traceback
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


def ratio_sum(numerators, denominators):
    """sum(numerators) / sum(denominators); nan unless the denominator is positive.

    Matches pixy >= 2.2.5: a negative numerator is a valid estimate, and a
    non-positive denominator is undefined.
    """
    denom = float(np.sum(denominators))
    if not np.isfinite(denom) or not denom > 0:
        return np.nan
    return float(np.sum(numerators) / denom)


def finite_rows(df, columns):
    """Rows whose component columns are all finite. Zeros are kept."""
    keep = pd.Series(True, index=df.index)
    for column in columns:
        values = pd.to_numeric(df[column], errors="coerce")
        keep &= values.notna()
    return df.loc[keep].copy()


def bootstrap_ratio_ci(numerators, denominators, B=1000, seed=42):
    """Bootstrap CI for sum(numerators) / sum(denominators) over resampled windows."""
    rng = np.random.default_rng(seed)
    n = len(numerators)
    boot = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        denom = float(denominators[idx].sum())
        boot[b] = np.nan if not (np.isfinite(denom) and denom > 0) else numerators[idx].sum() / denom
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


def parse_tajima_d_components(value):
    """Parse pixy >= 2.2.3 tajima_d_components ('nsum=<int>,mut=<int>')."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    if text in ("", "NA", "nan"):
        return None
    fields = {}
    for item in text.split(","):
        key, sep, val = item.partition("=")
        if not sep:
            raise ValueError(
                f"Unrecognized tajima_d_components {text!r}. "
                "pixy >= 2.2.3 writes 'nsum=<int>,mut=<int>' "
                "(pixy <= 2.2.1 wrote n:s allele-count classes). Re-run pixy."
            )
        fields[key] = int(float(val))
    if "nsum" not in fields or "mut" not in fields:
        raise ValueError(
            f"tajima_d_components {text!r} is missing nsum or mut"
        )
    return fields["nsum"], fields["mut"]


def calc_tajima_d_stdev(total_allele_count, num_sites, num_mutations):
    """Tajima's D denominator at the pooled mean allele count (pixy >= 2.2.3)."""
    if num_sites <= 0 or num_mutations <= 0:
        return float("nan")
    n = int(np.rint(total_allele_count / num_sites))
    if n < 2:
        return float("nan")
    e1, e2 = _tajima_constants(n)
    s = num_mutations
    variance = (e1 * s) + (e2 * s * (s - 1))
    if variance < 0:
        return float("nan")
    return math.sqrt(variance)


def aggregate_tajima_d(raw_pi, raw_watterson_theta, nsum, mut, no_sites):
    raw_pi_sum = float(np.nansum(raw_pi))
    raw_theta_sum = float(np.nansum(raw_watterson_theta))
    d_stdev = calc_tajima_d_stdev(
        float(np.nansum(nsum)),
        float(np.nansum(no_sites)),
        float(np.nansum(mut)),
    )
    if not np.isfinite(d_stdev) or d_stdev <= 0:
        tajima_d = np.nan
    else:
        tajima_d = (raw_pi_sum - raw_theta_sum) / d_stdev
    return tajima_d, d_stdev, raw_pi_sum, raw_theta_sum


def bootstrap_tajima_d_ci(raw_pi, raw_watterson_theta, nsum, mut, no_sites, B=1000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(raw_pi)
    boot = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        boot[b], _, _, _ = aggregate_tajima_d(
            raw_pi[idx],
            raw_watterson_theta[idx],
            nsum[idx],
            mut[idx],
            no_sites[idx],
        )
    valid = boot[~np.isnan(boot)]
    if len(valid) == 0:
        return (np.nan, np.nan), np.nan
    return np.percentile(valid, [2.5, 97.5]), float(np.std(valid, ddof=1)) if len(valid) > 1 else 0.0


def _require_columns(df, required, stat):
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{stat} output missing columns {sorted(missing)}; "
            "rerun pixy so window components are available for aggregation"
        )


def process_pi(file, output, bootstrap_replicates=1000):
    # sum(count_diffs) / sum(count_comparisons), not a mean of avg_pi.
    df = pd.read_csv(file, sep="\t")
    _require_columns(df, {"pop", "count_diffs", "count_comparisons", "no_sites"}, "pi")

    results = []
    for pop, group_df in df.groupby("pop"):
        used = finite_rows(group_df, ["count_diffs", "count_comparisons"])
        if used.empty:
            continue
        num = used["count_diffs"].to_numpy(dtype=float)
        den = used["count_comparisons"].to_numpy(dtype=float)
        sites = pd.to_numeric(used["no_sites"], errors="coerce").fillna(0).to_numpy(dtype=float)

        m = ratio_sum(num, den)
        (ci_low, ci_high), se = bootstrap_ratio_ci(num, den, B=bootstrap_replicates)

        results.append({
            "population": pop,
            "mean_pi": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_loci": len(used),
            "n_sites_total": int(sites.sum()),
        })

    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Pi results written to {output}")


def _fst_components(df):
    """Per-window numerator and denominator for Weir–Cockerham or Hudson FST."""
    if {"wc_fst_a", "wc_fst_b", "wc_fst_c"} <= set(df.columns):
        used = finite_rows(df, ["wc_fst_a", "wc_fst_b", "wc_fst_c"])
        if used.empty:
            return used, np.array([]), np.array([])
        a = used["wc_fst_a"].to_numpy(dtype=float)
        b = used["wc_fst_b"].to_numpy(dtype=float)
        c = used["wc_fst_c"].to_numpy(dtype=float)
        return used, a, a + b + c
    if {"hudson_fst_num", "hudson_fst_den"} <= set(df.columns):
        used = finite_rows(df, ["hudson_fst_num", "hudson_fst_den"])
        if used.empty:
            return used, np.array([]), np.array([])
        return (
            used,
            used["hudson_fst_num"].to_numpy(dtype=float),
            used["hudson_fst_den"].to_numpy(dtype=float),
        )
    raise ValueError(
        "fst output has no wc_fst_a/b/c or hudson_fst_num/den columns; "
        "rerun pixy with --fst_components"
    )


def process_fst(file, output, bootstrap_replicates=1000):
    # sum(components) recomputed as one ratio, not a mean of per-window FST.
    df = pd.read_csv(file, sep="\t")
    _require_columns(df, {"pop1", "pop2", "no_snps"}, "fst")
    has_wc = {"wc_fst_a", "wc_fst_b", "wc_fst_c"} <= set(df.columns)
    has_hudson = {"hudson_fst_num", "hudson_fst_den"} <= set(df.columns)
    if not has_wc and not has_hudson:
        raise ValueError(
            "fst output has no wc_fst_a/b/c or hudson_fst_num/den columns; "
            "rerun pixy with --fst_components"
        )

    results = []
    for (pop1, pop2), group_df in df.groupby(["pop1", "pop2"]):
        used, num, den = _fst_components(group_df)
        if used.empty:
            continue
        snps = pd.to_numeric(used["no_snps"], errors="coerce").fillna(0).to_numpy(dtype=float)

        m = ratio_sum(num, den)
        (ci_low, ci_high), se = bootstrap_ratio_ci(num, den, B=bootstrap_replicates)

        results.append({
            "pop1": pop1,
            "pop2": pop2,
            "mean_fst": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_loci": len(used),
            "n_snps_total": int(snps.sum()),
        })

    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Fst results written to {output}")


def process_dxy(file, output, bootstrap_replicates=1000):
    # sum(count_diffs) / sum(count_comparisons), not a mean of avg_dxy.
    df = pd.read_csv(file, sep="\t")
    _require_columns(df, {"pop1", "pop2", "count_diffs", "count_comparisons", "no_sites"}, "dxy")

    results = []
    for (pop1, pop2), group_df in df.groupby(["pop1", "pop2"]):
        used = finite_rows(group_df, ["count_diffs", "count_comparisons"])
        if used.empty:
            continue
        num = used["count_diffs"].to_numpy(dtype=float)
        den = used["count_comparisons"].to_numpy(dtype=float)
        sites = pd.to_numeric(used["no_sites"], errors="coerce").fillna(0).to_numpy(dtype=float)

        m = ratio_sum(num, den)
        (ci_low, ci_high), se = bootstrap_ratio_ci(num, den, B=bootstrap_replicates)

        results.append({
            "pop1": pop1,
            "pop2": pop2,
            "mean_dxy": m,
            "bootstrap_se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n_loci": len(used),
            "n_sites_total": int(sites.sum()),
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
    required = {"tajima_d", "raw_pi", "raw_watterson_theta", "no_sites", "tajima_d_components"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"tajima_d output missing columns {sorted(missing)}; "
            "run pixy >= 2.2.3 with --stats tajima_d --tajima_components"
        )
    parsed = df["tajima_d_components"].map(parse_tajima_d_components)
    df = df.loc[parsed.notna()].copy()
    df["_nsum"] = [pair[0] for pair in parsed.loc[parsed.notna()]]
    df["_mut"] = [pair[1] for pair in parsed.loc[parsed.notna()]]

    results = []
    for pop, group_df in df.groupby("pop"):
        raw_pi = pd.to_numeric(group_df["raw_pi"], errors="coerce").to_numpy(dtype=float)
        raw_theta = pd.to_numeric(group_df["raw_watterson_theta"], errors="coerce").to_numpy(dtype=float)
        nsum = group_df["_nsum"].to_numpy(dtype=float)
        mut = group_df["_mut"].to_numpy(dtype=float)
        sites = pd.to_numeric(group_df["no_sites"], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(raw_pi) & np.isfinite(raw_theta) & np.isfinite(sites)
        raw_pi, raw_theta, nsum, mut, sites = (
            raw_pi[ok], raw_theta[ok], nsum[ok], mut[ok], sites[ok]
        )
        if len(raw_pi) == 0:
            continue

        m, d_stdev, raw_pi_sum, raw_theta_sum = aggregate_tajima_d(
            raw_pi, raw_theta, nsum, mut, sites
        )
        (ci_low, ci_high), se = bootstrap_tajima_d_ci(
            raw_pi, raw_theta, nsum, mut, sites, B=bootstrap_replicates
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
            "n_loci": len(raw_pi),
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
