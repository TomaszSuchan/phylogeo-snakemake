#!/usr/bin/env python3

import argparse
import sys
import gzip
import warnings

try:
    import vcfpy
except ImportError:
    sys.exit("Error: This script requires vcfpy. Please install with:\npip install vcfpy")

# Suppress vcfpy warnings about missing contig length fields (common in RAD-seq data)
warnings.filterwarnings("ignore", message="Field \"length\" not found in header line contig.*", module="vcfpy")

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert VCF file to Structure format.')
    parser.add_argument('--vcf', required=True, help='Input VCF file (.vcf or .vcf.gz)')
    parser.add_argument('--out', required=True, help='Output Structure file')
    return parser.parse_args()

def convert_to_structure(vcf_file, output_file):
    """Convert VCF to Structure format"""
    
    # Structure nucleotide encoding
    STRUCTURE = {'A': '1',
                 'C': '2',
                 'T': '3',
                 'G': '4',
                 'N': '-9',
                 '.': '-9'
                }
    
    print(f"Reading {'compressed' if vcf_file.endswith('.gz') else 'uncompressed'} VCF file...")
    
    # Open VCF file
    if vcf_file.endswith('.gz'):
        file_handle = gzip.open(vcf_file, 'rt')
    else:
        file_handle = open(vcf_file, 'r')
    
    # Store all variant data
    variants = []
    sample_names = []
    
    try:
        reader = vcfpy.Reader(file_handle)
        sample_names = reader.header.samples.names
        sample_n = len(sample_names)
        
        print(f"Found {sample_n} samples in VCF file.")
        
        # Process each variant
        variant_count = 0
        for record in reader:
            variant_count += 1
            
            # Store genotype data for each sample
            variant_genotypes = []
            
            for call in record.calls:
                # Use gt_alleles directly instead of GT field
                alleles = call.gt_alleles
                
                if alleles and len(alleles) == 2:
                    # Convert allele indices to actual bases
                    allele1_idx = alleles[0]
                    allele2_idx = alleles[1]
                    
                    # Convert to actual bases
                    if allele1_idx is None:
                        allele1 = '.'
                    elif allele1_idx == 0:
                        allele1 = record.REF
                    else:
                        allele1 = record.ALT[0].value  # Only one ALT for biallelic
                    
                    if allele2_idx is None:
                        allele2 = '.'
                    elif allele2_idx == 0:
                        allele2 = record.REF
                    else:
                        allele2 = record.ALT[0].value  # Only one ALT for biallelic
                else:
                    # Missing genotype
                    allele1 = '.'
                    allele2 = '.'
                
                variant_genotypes.append([allele1, allele2])
            
            variants.append(variant_genotypes)
        
        reader.close()
        
        print(f"Processed {variant_count} variants.")
        
    except Exception as e:
        sys.exit(f"Error reading VCF file: {e}")
    finally:
        file_handle.close()
    
    # Write Structure format output
    print(f"Writing Structure format to: {output_file}")
    
    try:
        with open(output_file, 'w') as output:
            # Write each sample (two lines per sample for diploid data)
            for sample_idx, sample_name in enumerate(sample_names):
                # First allele line
                output.write(f"{sample_name}\t\t\t\t")
                for variant in variants:
                    allele1 = variant[sample_idx][0]
                    structure_code = STRUCTURE.get(allele1, '-9')
                    output.write(f"\t{structure_code}")
                output.write('\n')
                
                # Second allele line
                output.write(f"{sample_name}\t\t\t\t")
                for variant in variants:
                    allele2 = variant[sample_idx][1]
                    structure_code = STRUCTURE.get(allele2, '-9')
                    output.write(f"\t{structure_code}")
                output.write('\n')
        
        print(f"Successfully created Structure file with {len(sample_names)} samples and {len(variants)} variants.")
        
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")

def main():
    args = parse_arguments()
    
    print(f"Converting VCF to Structure format...")
    print(f"Input VCF: {args.vcf}")
    print(f"Output Structure file: {args.out}")
    
    try:
        convert_to_structure(args.vcf, args.out)
        print("\nConversion completed successfully!")
        
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()