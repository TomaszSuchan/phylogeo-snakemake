#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import traceback

# Redirect stdout and stderr to log file first so any error is visible in the rule log
log_file = snakemake.log[0]
log_handle = open(log_file, 'w')
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
            "n_sites_total": int(weights.sum())
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
            "n_snps_total": int(weights.sum())
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
            "n_sites_total": int(weights.sum())
        })
    
    out = pd.DataFrame(results)
    out.to_csv(output, sep="\t", index=False)
    print(f"Dxy results written to {output}")

# Which stat to process is set by the rule (params.stat)
try:
    stat = getattr(snakemake.params, "stat", None)
    if stat is None:
        # Fallback: infer from which input key is present
        for key in ("pi", "fst", "dxy"):
            if hasattr(snakemake.input, key):
                stat = key
                break
    if stat is None:
        raise ValueError("params.stat is required; add params: stat = 'pi'|'fst'|'dxy' to the rule")
    infile = getattr(snakemake.input, stat)
    outfile = getattr(snakemake.output, stat)

    if stat == "pi":
        process_pi(infile, outfile, bootstrap_replicates)
    elif stat == "fst":
        process_fst(infile, outfile, bootstrap_replicates)
    elif stat == "dxy":
        process_dxy(infile, outfile, bootstrap_replicates)
    else:
        raise ValueError(f"Unknown stat: {stat}")
except Exception:
    traceback.print_exc()
    raise