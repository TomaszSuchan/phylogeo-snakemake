#!/usr/bin/env python3
"""
Calculate VCF statistics:
- Number of chromosomes (unique values in column 1)
- Number of RAD fragments (unique first part of column 3, before "_pos")
- Number of variants (total rows)
- Number of samples (from VCF header)
"""

import sys
import gzip
import re
import os

def extract_rad_fragment(loc_id):
    """Extract RAD fragment ID from location ID (first part before '_pos')."""
    # Split by '_pos' and take the first part
    parts = loc_id.split('_pos')
    return parts[0] if parts else loc_id

def calculate_stats_vcf(vcf_file, output_file):
    """Calculate and write VCF statistics to output file."""
    chromosomes = set()
    rad_fragments = set()
    variant_count = 0
    num_samples = 0
    
    # Open VCF file (handle both gzipped and uncompressed)
    if vcf_file.endswith('.gz'):
        f = gzip.open(vcf_file, 'rt')
    else:
        f = open(vcf_file, 'r')
    
    try:
        with f:
            for line in f:
                # Check for #CHROM line to get sample count
                if line.startswith('#CHROM'):
                    # VCF format: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
                    # Samples start at column 9 (index 9)
                    fields = line.strip().split('\t')
                    if len(fields) > 9:
                        num_samples = len(fields) - 9  # Subtract 9 fixed columns
                    continue
                
                # Skip other header lines
                if line.startswith('#'):
                    continue
                
                # Parse VCF line
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                # Column 1: chromosome
                chrom = fields[0]
                chromosomes.add(chrom)
                
                # Column 3: location ID (extract RAD fragment)
                loc_id = fields[2]
                rad_fragment = extract_rad_fragment(loc_id)
                rad_fragments.add(rad_fragment)
                
                # Count variants
                variant_count += 1
    except Exception as e:
        print(f"Error reading VCF file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Get VCF file name and relative path
    vcf_filename = os.path.basename(vcf_file)
    vcf_abspath = os.path.abspath(vcf_file)
    
    # Calculate relative path from project root (assuming we're in results/{project}/stats_vcf/)
    # Try to find the project root by looking for common markers
    project_root = None
    path_parts = vcf_abspath.split(os.sep)
    if 'results' in path_parts:
        # Get everything up to and including 'results'
        results_idx = path_parts.index('results')
        project_root = os.sep.join(path_parts[:results_idx])
        vcf_relpath = os.path.relpath(vcf_abspath, project_root)
    else:
        # Fallback: use absolute path if we can't determine project root
        vcf_relpath = vcf_abspath
    
    # Write statistics to output file
    try:
        with open(output_file, 'w') as out:
            out.write(f"VCF file name: {vcf_filename}\n")
            out.write(f"VCF file location: {vcf_relpath}\n")
            out.write(f"Number of chromosomes: {len(chromosomes)}\n")
            out.write(f"Number of RAD fragments: {len(rad_fragments)}\n")
            out.write(f"Number of variants: {variant_count}\n")
            out.write(f"Number of samples: {num_samples}\n")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Statistics calculated:")
    print(f"  VCF file: {vcf_filename}")
    print(f"  VCF location: {vcf_relpath}")
    print(f"  Chromosomes: {len(chromosomes)}")
    print(f"  RAD fragments: {len(rad_fragments)}")
    print(f"  Variants: {variant_count}")
    print(f"  Samples: {num_samples}")

if __name__ == "__main__":
    # Check if running under Snakemake
    try:
        vcf_file = snakemake.input.vcf
        output_file = snakemake.output.stats
    except NameError:
        # Fallback to command-line arguments if not running under Snakemake
        if len(sys.argv) != 3:
            print("Usage: calculate_vcf_stats.py <input.vcf.gz> <output.txt>", file=sys.stderr)
            sys.exit(1)
        vcf_file = sys.argv[1]
        output_file = sys.argv[2]
    
    calculate_stats_vcf(vcf_file, output_file)

