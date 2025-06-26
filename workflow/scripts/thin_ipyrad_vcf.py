#!/usr/bin/env python3

import argparse
import random
import re
import sys
import gzip
from collections import defaultdict

try:
    import vcfpy
except ImportError:
    sys.exit("Error: This script requires vcfpy. Please install with:\npip install vcfpy")

def parse_arguments():
    parser = argparse.ArgumentParser(description='Select one SNP per RAD fragment from a VCF file.')
    parser.add_argument('--vcf', required=True, help='Input VCF file (.vcf or .vcf.gz)')
    parser.add_argument('--out', required=True, help='Output VCF file (.vcf for uncompressed, .vcf.gz for compressed)')
    parser.add_argument('--min-coverage', type=int, default=0, help='Minimum number of samples with data (default: 0)')
    parser.add_argument('--method', choices=['random', 'max_coverage', 'weighted'], default='max_coverage',
                        help='Selection method: random, max_coverage, or weighted by coverage (default: max_coverage)')
    parser.add_argument('--ties', choices=['first', 'random'], default='random',
                        help='How to handle ties in max_coverage method: take first SNP by position or random SNP (default)')
    parser.add_argument('--ns-tag', default='NS', help='INFO field tag for number of samples with data (default: NS)')
    parser.add_argument('--id-pattern', default=r'loc(\d+)_', 
                        help='Regex pattern to extract RAD fragment ID from variant ID (default: "loc(\\d+)_")')
    return parser.parse_args()

def get_coverage(record, ns_tag='NS'):
    """Extract sample coverage from variant record"""
    # Try the specified INFO tag first
    if ns_tag in record.INFO:
        value = record.INFO[ns_tag]
        if isinstance(value, list):
            return int(value[0]) if len(value) > 0 else 0
        return int(value) if value is not None else 0
    
    # Fallback: count non-missing genotypes
    coverage = 0
    for call in record.calls:
        if call.data.get('GT') and not call.is_het_missing() and not call.is_hom_missing():
            coverage += 1
    
    return coverage

def extract_rad_id(variant_id, pattern=r'loc(\d+)_'):
    """Extract RAD fragment ID from variant ID using regex pattern"""
    match = re.search(pattern, variant_id)
    if match:
        return match.group(1)  # Return the captured group
    return None

def process_vcf(vcf_file, min_coverage, ns_tag, id_pattern=r'loc(\d+)_'):
    """Process VCF grouping by RAD fragments"""
    # Group variants by RAD fragment
    rad_groups = defaultdict(list)
    variant_info = []
    unknown_ids = 0
    
    print(f"Reading {'compressed' if vcf_file.endswith('.gz') else 'uncompressed'} VCF file...")
    
    # Check if NS tag exists
    has_ns_tag = False
    variant_count = 0
    
    # Open VCF file
    if vcf_file.endswith('.gz'):
        file_handle = gzip.open(vcf_file, 'rt')
    else:
        file_handle = open(vcf_file, 'r')
    
    try:
        reader = vcfpy.Reader(file_handle)
        
        for record in reader:
            variant_count += 1
            chrom = record.CHROM
            pos = record.POS
            variant_id = record.ID[0] if record.ID else f"{chrom}:{pos}"
            
            # Check if NS tag exists (only need to check first few variants)
            if variant_count <= 10 and ns_tag in record.INFO:
                has_ns_tag = True
            
            coverage = get_coverage(record, ns_tag)
            
            # Store variant info
            variant_info.append({
                'chrom': chrom,
                'pos': pos,
                'id': variant_id,
                'coverage': coverage,
                'index': len(variant_info)
            })
            
            # Skip variants not meeting coverage requirement
            if coverage < min_coverage:
                continue
            
            # Extract RAD fragment ID
            rad_id = extract_rad_id(variant_id, id_pattern)
            if rad_id:
                group_id = f"rad_{rad_id}"
                rad_groups[group_id].append((chrom, pos, variant_id, coverage, len(variant_info) - 1))
            else:
                unknown_ids += 1
        
        reader.close()
        
    except Exception as e:
        sys.exit(f"Error reading VCF file: {e}")
    finally:
        file_handle.close()
    
    # Report coverage detection method
    if has_ns_tag:
        print(f"Using '{ns_tag}' tag from INFO fields for coverage calculation")
    else:
        print(f"Note: '{ns_tag}' tag not found in INFO fields, using genotype-based coverage calculation")
    
    # Sort SNPs within each RAD fragment group by chromosome and position
    print("Sorting SNPs within RAD fragments by position...")
    for group_id in rad_groups:
        # Sort by chromosome first, then by position
        rad_groups[group_id].sort(key=lambda x: (x[0], x[1]))
    
    return vcf_file, variant_info, rad_groups, unknown_ids

def select_snps(groups, method, ties):
    """
    Select one SNP per RAD fragment using the specified method.
    
    Parameters:
    -----------
    groups : dict
        SNP groups indexed by RAD fragment ID (sorted by position within each group)
    method : str
        Selection method: 'random', 'max_coverage', or 'weighted'
    ties : str
        How to handle ties in max_coverage: 'first' (earliest position) or 'random'
    """
    selected_indices = set()
    selected_data = []
    
    for group_id, snps in groups.items():
        if not snps:
            continue
            
        if method == "random":
            # Pure random selection
            chosen = random.choice(snps)
        elif method == "max_coverage":
            # Find maximum coverage value
            max_cov = max(snp[3] for snp in snps)
            
            # Get all SNPs with this maximum coverage
            max_snps = [snp for snp in snps if snp[3] == max_cov]
            
            if len(max_snps) == 1:
                # Only one SNP has the max coverage
                chosen = max_snps[0]
            else:
                # Handle ties based on ties parameter
                if ties == "random":
                    # Randomly select among SNPs with maximum coverage
                    chosen = random.choice(max_snps)
                else:  # "first" is the default
                    # Take the first SNP with maximum coverage (earliest position due to sorting)
                    chosen = max_snps[0]
        elif method == "weighted":
            # Weight by coverage (higher coverage = higher chance)
            weights = [snp[3] for snp in snps]
            total = sum(weights)
            if total > 0:
                probs = [w/total for w in weights]
                chosen = random.choices(snps, weights=probs, k=1)[0]
            else:
                chosen = random.choice(snps)
                
        # Add the index of the chosen SNP to selected indices
        selected_indices.add(chosen[4])
        selected_data.append(chosen)
    
    return selected_indices, selected_data

def write_output(vcf_file, variant_info, selected_indices, output_file):
    """Write output VCF using vcfpy"""
    compress = output_file.endswith('.gz')
    
    # Sort selected indices for efficient processing
    selected_set = set(selected_indices)
    
    # Open input VCF file
    if vcf_file.endswith('.gz'):
        input_handle = gzip.open(vcf_file, 'rt')
    else:
        input_handle = open(vcf_file, 'r')
    
    # Open output VCF file
    if compress:
        output_handle = gzip.open(output_file, 'wt')
    else:
        output_handle = open(output_file, 'w')
    
    try:
        reader = vcfpy.Reader(input_handle)
        writer = vcfpy.Writer(output_handle, reader.header)
        
        current_index = 0
        for record in reader:
            if current_index in selected_set:
                writer.write_record(record)
            current_index += 1
        
        reader.close()
        writer.close()
        
    except Exception as e:
        sys.exit(f"Error writing output VCF: {e}")
    finally:
        input_handle.close()
        output_handle.close()
    
    return output_file, compress

def main():
    args = parse_arguments()
    
    print(f"Parsing VCF file: {args.vcf} (grouping by RAD fragment)...")
    
    try:
        vcf_file, variant_info, rad_groups, unknown_ids = process_vcf(
            args.vcf, args.min_coverage, args.ns_tag, args.id_pattern)
        
        print(f"Found {len(variant_info)} SNPs in the VCF file.")
        print(f"Found {len(rad_groups)} RAD fragments with SNPs meeting the criteria.")
        
        if unknown_ids > 0:
            print(f"Warning: {unknown_ids} SNPs could not be assigned to a RAD fragment using pattern '{args.id_pattern}'")
        
        ties_method = ""
        if args.method == "max_coverage" and args.ties == "first":
            ties_method = ", resolving ties by taking SNP with earliest position"
        elif args.method == "max_coverage" and args.ties == "random":
            ties_method = ", resolving ties randomly"
        
        print(f"Selecting SNPs using the '{args.method}' method{ties_method}...")
        
        selected_indices, selected_data = select_snps(rad_groups, args.method, args.ties)
        print(f"Selected {len(selected_indices)} SNPs (one per RAD fragment).")
        
        # Write output VCF
        output_file, is_compressed = write_output(vcf_file, variant_info, selected_indices, args.out)
        
        compression_status = "compressed" if is_compressed else "uncompressed"
        print(f"Created {compression_status} filtered VCF: {output_file}")
        
        # Print summary statistics
        coverage_values = [snp[3] for snp in selected_data]
        if coverage_values:
            avg_coverage = sum(coverage_values) / len(coverage_values)
            min_coverage = min(coverage_values)
            max_coverage = max(coverage_values)
            print(f"\nSummary of selected SNPs:")
            print(f"  Average sample coverage: {avg_coverage:.2f}")
            print(f"  Minimum sample coverage: {min_coverage}")
            print(f"  Maximum sample coverage: {max_coverage}")
            
            # Count how many chromosomes are represented
            chroms = {snp[0] for snp in selected_data}
            print(f"  SNPs distributed across {len(chroms)} chromosomes")
        
        print("\nDone!")
        
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()