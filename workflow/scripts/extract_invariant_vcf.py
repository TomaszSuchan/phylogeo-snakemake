#!/usr/bin/env python3
"""
Extract invariant sites from RADseq loci file and generate VCF
A site is considered invariant if all samples with valid nucleotides (A,T,G,C)
have the same nucleotide, regardless of samples with N or ambiguous codes
"""

import sys
from collections import defaultdict
import argparse
from datetime import datetime

def parse_loci_file(filename):
    """Parse the loci file and return sequences by locus"""
    loci_data = []
    current_locus_data = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if '//' in line:  # Changed to handle lines that contain //
                # End of current locus
                if current_locus_data:
                    loci_data.append(current_locus_data)
                    current_locus_data = {}
            else:
                # Parse sequence line - only if it contains sequence data
                parts = line.split()
                if len(parts) >= 2 and not line.startswith('|'):
                    sample_id = parts[0]
                    sequence = parts[1]
                    # Only add if it looks like a sequence (contains ATGCN etc)
                    if any(c in sequence.upper() for c in 'ATGCN'):
                        current_locus_data[sample_id] = sequence.upper()
        
        # Don't forget the last locus if file doesn't end with //
        if current_locus_data:
            loci_data.append(current_locus_data)
    
    return loci_data

def find_invariant_sites(sequences):
    """Find positions where all samples with valid nucleotides have the same nucleotide"""
    if not sequences:
        return []
    
    # Get sequence length (assuming all are aligned and same length)
    seq_length = len(next(iter(sequences.values())))
    invariant_sites = []
    
    # Valid nucleotides (excluding ALL ambiguity codes)
    valid_nucs = {'A', 'T', 'G', 'C'}
    
    # IUPAC ambiguity codes: any of these makes a site variable
    # R=purine, Y=pyrimidine, S=strong, W=weak, K=keto, M=amino,
    # B=not A, D=not C, H=not G, V=not T, N=any
    ambiguity_codes = {'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'}
    
    for pos in range(seq_length):
        # Get all valid nucleotides at this position
        nucs_at_pos = []
        ambiguous_chars = []  # Track ambiguity codes at this position
        
        # Check ALL samples at this position
        for sample_id, seq in sequences.items():
            if pos < len(seq):
                char = seq[pos]
                if char in valid_nucs:
                    nucs_at_pos.append(char)
                elif char in ambiguity_codes:
                    ambiguous_chars.append(char)
        
        # Position is invariant if:
        # 1. At least one sample has a valid nucleotide
        # 2. All valid nucleotides are the same
        # 3. No ambiguity codes present (any ambiguity code indicates uncertainty and makes site variable)
        if nucs_at_pos and len(set(nucs_at_pos)) == 1 and len(ambiguous_chars) == 0:
            invariant_sites.append((pos, nucs_at_pos[0]))
    
    return invariant_sites

def read_samples_file(filename):
    """Read sample names from a file (one per line)"""
    samples = set()
    with open(filename, 'r') as f:
        for line in f:
            sample = line.strip()
            if sample:  # Skip empty lines
                samples.add(sample)
    return samples

def generate_vcf(loci_data, output_file, samples_filter=None):
    """Generate VCF file with invariant sites only
    
    Args:
        loci_data: List of dictionaries with sequences by locus
        output_file: Output VCF file path
        samples_filter: Optional set of sample names to include (only samples in both loci and filter)
    """
    
    # Get all unique sample names across all loci
    all_samples = set()
    for locus_seqs in loci_data:
        all_samples.update(locus_seqs.keys())
    
    # Apply filter if provided
    if samples_filter is not None:
        # Only include samples that are in both the loci file and the filter list
        all_samples = all_samples.intersection(samples_filter)
    
    samples_list = sorted(list(all_samples))
    
    with open(output_file, 'w') as vcf:
        # VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        vcf.write("##source=extract_invariant_sites.py\n")
        vcf.write("##INFO=<ID=INVARIANT,Number=0,Type=Flag,Description=\"Invariant site\">\n")
        vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        
        # Write column headers
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for sample in samples_list:
            vcf.write(f"\t{sample}")
        vcf.write("\n")
        
        # Process each locus
        for locus_idx, sequences in enumerate(loci_data):
            locus_num = locus_idx  # 0-based numbering
            
            # Filter sequences per locus to only process samples in the filter
            if samples_filter is not None:
                filtered_sequences = {sample: seq for sample, seq in sequences.items() 
                                    if sample in samples_filter}
            else:
                filtered_sequences = sequences
            
            invariant_sites = find_invariant_sites(filtered_sequences)
            
            for pos, ref_allele in invariant_sites:
                # Chromosome name
                chrom = f"RAD_{locus_num}"
                
                # Position (1-based in VCF)
                vcf_pos = pos + 1
                
                # ID formatted as requested
                var_id = f"loc{locus_num}_pos{pos}"
                
                # Write variant line
                vcf.write(f"{chrom}\t{vcf_pos}\t{var_id}\t{ref_allele}\t.\t.\tPASS\tINVARIANT\tGT")
                
                # Write genotypes for all samples
                # Use filtered_sequences if filter was applied, otherwise use original sequences
                sequences_to_use = filtered_sequences if samples_filter is not None else sequences
                
                for sample in samples_list:
                    if sample in sequences_to_use and pos < len(sequences_to_use[sample]):
                        nuc = sequences_to_use[sample][pos]
                        if nuc == ref_allele:
                            vcf.write("\t0/0")
                        elif nuc in {'A', 'T', 'G', 'C'}:
                            # This shouldn't happen for invariant sites
                            vcf.write("\t./.")
                        else:
                            # N or ambiguous code
                            vcf.write("\t./.")
                    else:
                        # Sample not in this locus
                        vcf.write("\t./.")
                
                vcf.write("\n")

def main():
    parser = argparse.ArgumentParser(
        description='Extract invariant sites from RADseq loci file and generate VCF.\n'
                    'A site is considered invariant if all samples with valid nucleotides (A,T,G,C)\n'
                    'have the same nucleotide, regardless of samples with N or ambiguous codes.'
    )
    parser.add_argument('input_file', help='Input loci file')
    parser.add_argument('-o', '--output', default='invariant_sites.vcf', 
                       help='Output VCF file (default: invariant_sites.vcf)')
    parser.add_argument('--samples-file', help='File containing sample names to include (one per line)')
    
    args = parser.parse_args()
    
    # Read samples file if provided
    samples_filter = None
    if args.samples_file:
        samples_filter = read_samples_file(args.samples_file)
        if not samples_filter:
            print(f"ERROR: No samples found in filter file: {args.samples_file}", file=sys.stderr)
            sys.exit(1)
    
    # Parse the loci file
    loci_data = parse_loci_file(args.input_file)
    
    # Check for critical errors
    if not loci_data:
        print(f"ERROR: No loci found in input file: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Collect all samples in loci
    all_samples_in_loci = set()
    for locus_seqs in loci_data:
        all_samples_in_loci.update(locus_seqs.keys())
    
    if not all_samples_in_loci:
        print(f"ERROR: No samples found in loci file: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Check if samples in filter file are missing from loci
    if samples_filter:
        samples_in_filter_not_in_loci = samples_filter - all_samples_in_loci
        if samples_in_filter_not_in_loci:
            print(f"ERROR: {len(samples_in_filter_not_in_loci)} sample(s) in filter file not found in loci file:", file=sys.stderr)
            for sample in sorted(samples_in_filter_not_in_loci):
                print(f"  {sample}", file=sys.stderr)
            sys.exit(1)
        
        # Check if any samples will be included after filtering
        filtered_samples = all_samples_in_loci.intersection(samples_filter)
        if not filtered_samples:
            print(f"ERROR: No samples match between filter file and loci file", file=sys.stderr)
            sys.exit(1)
    
    # Generate VCF
    generate_vcf(loci_data, args.output, samples_filter=samples_filter)
    
    print(f"VCF file written to: {args.output}")

if __name__ == "__main__":
    main()
