#!/usr/bin/env python

import os
import sys

import pandas as pd
import vcfpy

# Snakemake automatically passes input/output/log via the 'snakemake' object
vcf_path = snakemake.input.vcf
popdata_path = snakemake.params.popdata
input_popmap_path = snakemake.params.popmap
separator = snakemake.params.separator
output_path = snakemake.output.indpopdata
log_file = snakemake.log[0]

# Redirect stdout and stderr to log file
sys.stdout = open(log_file, "w")
sys.stderr = sys.stdout


def load_vcf_samples(vcf_file):
    """Read samples from a VCF while preserving the header order."""
    try:
        reader = vcfpy.Reader.from_path(vcf_file)
        samples = reader.header.samples.names
        reader.close()
        return samples
    except Exception as exc:
        print(f"ERROR: Failed to read samples from VCF file: {vcf_file}")
        print(f"  {exc}")
        sys.exit(1)


def build_site_mapping(vcf_samples, popmap_file, popseparator):
    """Return a dataframe with Ind/Site in VCF sample order."""
    if popmap_file and popmap_file.strip():
        print(f"Using configured popmap file for Site assignments: {popmap_file}")
        try:
            popmap_df = pd.read_csv(
                popmap_file,
                sep="\t",
                header=None,
                names=["Ind", "Site"],
            )
        except FileNotFoundError:
            print(f"ERROR: Popmap file not found: {popmap_file}")
            sys.exit(1)
        except Exception as exc:
            print(f"ERROR: Failed to read popmap file: {popmap_file}")
            print(f"  {exc}")
            sys.exit(1)

        popmap_dict = dict(zip(popmap_df["Ind"], popmap_df["Site"]))
        missing_samples = [sample for sample in vcf_samples if sample not in popmap_dict]
        if missing_samples:
            print("\n" + "=" * 80)
            print("ERROR: Samples found in VCF but missing from configured popmap")
            print("=" * 80)
            print(f"\nVCF file:    {vcf_path}")
            print(f"Popmap file: {popmap_file}")
            print(f"\nMissing {len(missing_samples)} sample(s):\n")
            for sample in missing_samples:
                print(sample)
            print("\nPlease add these samples to the configured popmap file.")
            print("=" * 80)
            sys.exit(1)

        output_rows = [[sample, popmap_dict[sample]] for sample in vcf_samples]
        return pd.DataFrame(output_rows, columns=["Ind", "Site"])

    print("No configured popmap provided. Deriving Site assignments from sample names.")
    print(f"Using separator: '{popseparator}'")
    output_rows = []
    for sample in vcf_samples:
        parts = sample.split(popseparator)
        population = parts[0] if len(parts) > 1 else "UNKNOWN"
        output_rows.append([sample, population])
    return pd.DataFrame(output_rows, columns=["Ind", "Site"])


print(f"Reading VCF file: {vcf_path}")
vcf_samples = load_vcf_samples(vcf_path)
print(f"  Found {len(vcf_samples)} samples in VCF")

ind_site_df = build_site_mapping(vcf_samples, input_popmap_path, separator)
print(f"  Built Site assignments for {len(ind_site_df)} individuals")
print(f"  Found {ind_site_df['Site'].nunique()} unique sites")

# Check if popdata file is provided and exists
if popdata_path and popdata_path.strip() and os.path.exists(popdata_path):
    print(f"\nReading popdata file: {popdata_path}")
    popdata_df = pd.read_csv(popdata_path, sep="\t", header=0)
    print(f"  Found {len(popdata_df)} rows in popdata (before deduplication)")

    if "Site" not in popdata_df.columns:
        print("ERROR: popdata file must contain a 'Site' column")
        sys.exit(1)

    n_before = len(popdata_df)
    popdata_df = popdata_df.drop_duplicates(subset=["Site"], keep="first")
    n_after = len(popdata_df)
    if n_before > n_after:
        print(f"  Removed {n_before - n_after} duplicate site entries")
    print(f"  Found {n_after} unique sites in popdata")

    unique_sites_ind = set(ind_site_df["Site"].unique())
    unique_sites_popdata = set(popdata_df["Site"].unique())
    missing_sites = unique_sites_ind - unique_sites_popdata

    if missing_sites:
        print("\n" + "=" * 80)
        print("ERROR: Sites found in VCF/population assignments but missing from popdata file")
        print("=" * 80)
        print(f"\nPopdata file: {popdata_path}")
        print(f"VCF file:     {vcf_path}")
        print(f"\nMissing {len(missing_sites)} site(s):\n")
        for site in sorted(missing_sites):
            print(site)
        print("\nPlease add these sites to the popdata file or update the Site assignments.")
        print("=" * 80)
        sys.exit(1)

    merged_df = pd.merge(ind_site_df, popdata_df, on="Site", how="left")
    output_columns = ["Ind", "Site"] + [col for col in merged_df.columns if col not in ["Ind", "Site"]]
    final_popdata_df = merged_df[output_columns]
    print(
        f"\nGenerated indpopdata with {len(final_popdata_df.columns)} columns: "
        f"{', '.join(final_popdata_df.columns)}"
    )
else:
    print("\nNo popdata file provided or file does not exist.")
    print("Generating indpopdata with only Ind and Site columns.")
    final_popdata_df = ind_site_df[["Ind", "Site"]]
    print(f"Generated indpopdata with {len(final_popdata_df.columns)} columns: {', '.join(final_popdata_df.columns)}")

final_popdata_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"\nSaved indpopdata to: {output_path}")