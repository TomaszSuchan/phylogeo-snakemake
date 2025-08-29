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

def find_invariant_sites(sequences, verbose=False):
    """Find positions where all samples with valid nucleotides have the same nucleotide"""
    if not sequences:
        return []
    
    # Get sequence length (assuming all are aligned and same length)
    seq_length = len(next(iter(sequences.values())))
    invariant_sites = []
    
    # Valid nucleotides (excluding ALL ambiguity codes)
    valid_nucs = {'A', 'T', 'G', 'C'}
    
    # Track statistics for debugging
    position_stats = []
    
    for pos in range(seq_length):
        # Get all valid nucleotides at this position
        nucs_at_pos = []
        invalid_count = 0
        
        # Check ALL samples at this position
        for sample_id, seq in sequences.items():
            if pos < len(seq):
                if seq[pos] in valid_nucs:
                    nucs_at_pos.append(seq[pos])
                else:
                    invalid_count += 1
            else:
                invalid_count += 1
        
        # Position is invariant if:
        # 1. At least one sample has a valid nucleotide
        # 2. All valid nucleotides are the same
        if nucs_at_pos and len(set(nucs_at_pos)) == 1:
            invariant_sites.append((pos, nucs_at_pos[0]))
        
        # Store stats for debugging
        if verbose and pos < 10:  # Show first 10 positions
            unique_nucs = set(nucs_at_pos) if nucs_at_pos else set()
            position_stats.append({
                'pos': pos,
                'valid_count': len(nucs_at_pos),
                'invalid_count': invalid_count,
                'unique_nucs': unique_nucs
            })
    
    if verbose and position_stats:
        print(f"    First 10 positions analysis:")
        for stat in position_stats:
            print(f"      Pos {stat['pos']}: {stat['valid_count']} valid samples, "
                  f"{stat['invalid_count']} with N/ambiguous, nucleotides: {stat['unique_nucs']}")
    
    return invariant_sites

def generate_vcf(loci_data, output_file):
    """Generate VCF file with invariant sites only"""
    
    # Get all unique sample names across all loci
    all_samples = set()
    for locus_seqs in loci_data:
        all_samples.update(locus_seqs.keys())
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
            invariant_sites = find_invariant_sites(sequences)
            
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
                for sample in samples_list:
                    if sample in sequences and pos < len(sequences[sample]):
                        nuc = sequences[sample][pos]
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
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Reading loci file: {args.input_file}")
    
    # Parse the loci file
    loci_data = parse_loci_file(args.input_file)
    
    if args.verbose:
        print(f"Found {len(loci_data)} loci")
        total_invariant = 0
        
        for locus_idx, sequences in enumerate(loci_data):
            locus_num = locus_idx  # 0-based numbering
            
            # Get sequence length and sample info
            if sequences:
                seq_length = len(next(iter(sequences.values())))
                # Show first few samples for debugging
                sample_list = list(sequences.keys())[:3]
                print(f"\n  Locus {locus_num} (RAD_{locus_num}):")
                print(f"    Samples: {len(sequences)} (first few: {', '.join(sample_list)}...)")
                print(f"    Sequence length: {seq_length} bp")
                
                # Find invariant sites with verbose output
                invariant_sites = find_invariant_sites(sequences, verbose=True)
                total_invariant += len(invariant_sites)
                print(f"    Invariant sites found: {len(invariant_sites)}")
            else:
                print(f"  Locus {locus_num} (RAD_{locus_num}): Empty")
        
        print(f"\nTotal invariant sites across all loci: {total_invariant}")
    
    # Generate VCF
    generate_vcf(loci_data, args.output)
    
    print(f"VCF file written to: {args.output}")

if __name__ == "__main__":
    main()
