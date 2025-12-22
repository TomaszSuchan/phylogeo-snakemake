#!/usr/bin/env python

import pandas as pd
import sys
import vcfpy

# Snakemake automatically passes input/output/params via the 'snakemake' object
vcf_file = snakemake.input.vcf
output_popmap = snakemake.output.popmap
input_popmap = snakemake.params.popmap
separator = snakemake.params.separator
log_file = snakemake.log[0]

# Redirect stdout and stderr to log file
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Generating popmap from VCF: {vcf_file}")

# Get list of samples from VCF using vcfpy
try:
    reader = vcfpy.Reader.from_path(vcf_file)
    vcf_samples = reader.header.samples.names
    reader.close()
    print(f"  Found {len(vcf_samples)} samples in VCF")
except Exception as e:
    print(f"ERROR: Failed to read samples from VCF file: {vcf_file}")
    print(f"  {e}")
    sys.exit(1)

if input_popmap and input_popmap != "":
    # User provided a popmap file - validate against VCF samples
    print(f"Validating provided popmap: {input_popmap}")

    try:
        # Read the provided popmap file
        popmap_df = pd.read_csv(input_popmap, sep="\t", header=None, names=["Ind", "Site"])
        print(f"  Found {len(popmap_df)} entries in popmap")

        # Create a lookup dictionary: sample -> population
        popmap_dict = dict(zip(popmap_df["Ind"], popmap_df["Site"]))

        # Check which VCF samples are missing from popmap
        missing_samples = [sample for sample in vcf_samples if sample not in popmap_dict]

        if missing_samples:
            print("\n" + "=" * 80)
            print("ERROR: Samples found in VCF but missing from popmap file")
            print("=" * 80)
            print(f"\nVCF file:    {vcf_file}")
            print(f"Popmap file: {input_popmap}")
            print(f"\nMissing {len(missing_samples)} sample(s):\n")
            for sample in missing_samples:
                print(f"  - {sample}")
            print("\n" + "=" * 80)
            print("Please add these samples to the popmap file with their population assignments.")
            print("Example format (tab-separated):")
            print("  sample_name<TAB>population_name")
            print("=" * 80)
            sys.exit(1)

        # Create output popmap with VCF sample order
        output_data = []
        for sample in vcf_samples:
            output_data.append([sample, popmap_dict[sample]])

        output_df = pd.DataFrame(output_data, columns=["Ind", "Site"])
        output_df.to_csv(output_popmap, sep="\t", header=False, index=False)

        print(f"  Successfully matched {len(output_df)} samples")
        print(f"  Popmap generated with samples in VCF order")

    except FileNotFoundError:
        print(f"ERROR: Popmap file not found: {input_popmap}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read popmap file: {input_popmap}")
        print(f"  {e}")
        sys.exit(1)

else:
    # No popmap provided - generate from VCF sample names
    print("No popmap provided, generating from VCF sample names...")
    print(f"  Using separator: '{separator}'")

    output_data = []
    for sample in vcf_samples:
        # Split by separator and use first part as population
        parts = sample.split(separator)
        population = parts[0] if len(parts) > 1 else "UNKNOWN"
        output_data.append([sample, population])

    output_df = pd.DataFrame(output_data, columns=["Ind", "Site"])
    output_df.to_csv(output_popmap, sep="\t", header=False, index=False)

    print(f"  Generated popmap for {len(output_df)} samples")

print(f"\nPopmap written to: {output_popmap}")
