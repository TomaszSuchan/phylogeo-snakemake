#!/usr/bin/env python3
"""
Script to summarize Hardy-Weinberg equilibrium test results from vcftools --hardy output.
Provides summary statistics and identifies loci with significant deviations from HWE.
"""

import argparse
import sys
import numpy as np
import gzip
import re
from collections import defaultdict, Counter


def extract_rad_locus(variant_id):
    """
    Extract RAD locus ID from ipyrad variant ID.
    Examples: loc0_pos149 -> loc0, loc1_pos105 -> loc1
    """
    match = re.match(r'(loc\d+)_pos\d+', variant_id)
    if match:
        return match.group(1)
    return None


def build_vcf_position_map(vcf_file):
    """
    Build a mapping from (CHR, POS) to variant ID and RAD locus.
    """
    position_map = {}
    opener = gzip.open if vcf_file.endswith('.gz') else open

    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = fields[1]
            variant_id = fields[2]

            # Extract RAD locus from variant ID
            rad_locus = extract_rad_locus(variant_id)

            position_map[(chrom, pos)] = {
                'id': variant_id,
                'rad_locus': rad_locus
            }

    return position_map


def parse_hardy_file(hardy_file, position_map=None):
    """
    Parse vcftools .hwe file.

    Format:
    CHR  POS     OBS(HOM1/HET/HOM2)    E(HOM1/HET/HOM2)    ChiSq_HWE   P_HWE   P_HET_DEFICIT   P_HET_EXCESS
    """
    loci = []

    with open(hardy_file, 'r') as f:
        header = f.readline().strip()

        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 6:
                    chrom = parts[0]
                    pos = parts[1]
                    obs_counts = parts[2]
                    exp_counts = parts[3]
                    chisq = float(parts[4]) if parts[4] != 'nan' else np.nan
                    p_hwe = float(parts[5]) if parts[5] != 'nan' else np.nan
                    p_deficit = float(parts[6]) if len(parts) > 6 and parts[6] != 'nan' else np.nan
                    p_excess = float(parts[7]) if len(parts) > 7 and parts[7] != 'nan' else np.nan

                    # Get variant ID and RAD locus from position map if available
                    variant_id = None
                    rad_locus = None
                    if position_map and (chrom, pos) in position_map:
                        variant_id = position_map[(chrom, pos)]['id']
                        rad_locus = position_map[(chrom, pos)]['rad_locus']

                    loci.append({
                        'chrom': chrom,
                        'pos': pos,
                        'variant_id': variant_id,
                        'rad_locus': rad_locus,
                        'obs_counts': obs_counts,
                        'exp_counts': exp_counts,
                        'chisq': chisq,
                        'p_hwe': p_hwe,
                        'p_deficit': p_deficit,
                        'p_excess': p_excess
                    })

    return loci


def summarize_hwe(loci, alpha=0.05):
    """
    Generate summary statistics for HWE test results.
    """
    # Filter out loci with nan p-values
    valid_loci = [l for l in loci if not np.isnan(l['p_hwe'])]

    total_loci = len(loci)
    valid_loci_count = len(valid_loci)

    if valid_loci_count == 0:
        return {
            'total_loci': total_loci,
            'valid_loci': 0,
            'sig_deviations': 0,
            'sig_deficit': 0,
            'sig_excess': 0,
            'mean_p': np.nan,
            'median_p': np.nan,
            'top_deviations': [],
            'rad_locus_stats': {}
        }

    # Count significant deviations
    sig_deviations = sum(1 for l in valid_loci if l['p_hwe'] < alpha)
    sig_deficit = sum(1 for l in valid_loci if not np.isnan(l['p_deficit']) and l['p_deficit'] < alpha)
    sig_excess = sum(1 for l in valid_loci if not np.isnan(l['p_excess']) and l['p_excess'] < alpha)

    # Calculate statistics
    p_values = [l['p_hwe'] for l in valid_loci]
    mean_p = np.mean(p_values)
    median_p = np.median(p_values)

    # Find most significant deviations (top 100)
    sorted_loci = sorted(valid_loci, key=lambda x: x['p_hwe'])
    top_deviations = sorted_loci[:100]  # Top 100 most significant

    # Analyze RAD locus statistics
    rad_locus_stats = {}
    loci_with_rad = [l for l in valid_loci if l.get('rad_locus') is not None]

    if loci_with_rad:
        # Count total SNPs per RAD locus
        rad_locus_total = Counter(l['rad_locus'] for l in loci_with_rad)

        # Count significant deviations per RAD locus
        rad_locus_sig = Counter(l['rad_locus'] for l in loci_with_rad if l['p_hwe'] < alpha)

        # Calculate proportion of significant SNPs per locus
        for rad_locus, total_snps in rad_locus_total.items():
            sig_snps = rad_locus_sig.get(rad_locus, 0)
            rad_locus_stats[rad_locus] = {
                'total_snps': total_snps,
                'sig_snps': sig_snps,
                'prop_sig': sig_snps / total_snps if total_snps > 0 else 0
            }

    return {
        'total_loci': total_loci,
        'valid_loci': valid_loci_count,
        'sig_deviations': sig_deviations,
        'sig_deficit': sig_deficit,
        'sig_excess': sig_excess,
        'mean_p': mean_p,
        'median_p': median_p,
        'top_deviations': top_deviations,
        'rad_locus_stats': rad_locus_stats
    }


def write_report(summary, output_file, alpha=0.05, vcf_file=None):
    """
    Write a formatted summary report.
    """
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Hardy-Weinberg Equilibrium Test Summary\n")
        f.write("=" * 80 + "\n\n")

        if vcf_file:
            f.write(f"Input VCF: {vcf_file}\n")
        f.write(f"Significance threshold (alpha): {alpha}\n\n")

        f.write("=" * 80 + "\n")
        f.write("Overall Statistics\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Total loci tested: {summary['total_loci']}\n")
        f.write(f"Valid loci (with p-values): {summary['valid_loci']}\n")
        f.write(f"Loci with significant HWE deviation (p < {alpha}): {summary['sig_deviations']} ")
        f.write(f"({100.0 * summary['sig_deviations'] / summary['valid_loci']:.2f}%)\n")
        f.write(f"Loci with heterozygote deficit (p < {alpha}): {summary['sig_deficit']} ")
        f.write(f"({100.0 * summary['sig_deficit'] / summary['valid_loci']:.2f}%)\n")
        f.write(f"Loci with heterozygote excess (p < {alpha}): {summary['sig_excess']} ")
        f.write(f"({100.0 * summary['sig_excess'] / summary['valid_loci']:.2f}%)\n\n")

        f.write(f"Mean p-value: {summary['mean_p']:.6f}\n")
        f.write(f"Median p-value: {summary['median_p']:.6f}\n\n")

        f.write("=" * 80 + "\n")
        f.write("Interpretation\n")
        f.write("=" * 80 + "\n\n")

        deviation_pct = 100.0 * summary['sig_deviations'] / summary['valid_loci']

        if deviation_pct < 5:
            f.write("GOOD: Few loci deviate from HWE expectations. This suggests:\n")
            f.write("  - High genotyping accuracy\n")
            f.write("  - Random mating within populations\n")
            f.write("  - No strong selection at most loci\n\n")
        elif deviation_pct < 10:
            f.write("MODERATE: Some loci deviate from HWE. This may indicate:\n")
            f.write("  - Minor genotyping errors\n")
            f.write("  - Population substructure\n")
            f.write("  - Selection at some loci\n")
            f.write("  - Consider filtering loci with extreme deviations\n\n")
        else:
            f.write("WARNING: Many loci deviate from HWE. This suggests:\n")
            f.write("  - Potential genotyping issues (e.g., allelic dropout)\n")
            f.write("  - Strong population substructure\n")
            f.write("  - Selection or null alleles\n")
            f.write("  - Review genotyping quality and population structure\n\n")

        if summary['sig_deficit'] > summary['sig_excess'] * 2:
            f.write("NOTE: Heterozygote deficit is more common than excess.\n")
            f.write("This pattern often indicates:\n")
            f.write("  - Null alleles\n")
            f.write("  - Allelic dropout\n")
            f.write("  - Inbreeding or Wahlund effect\n\n")
        elif summary['sig_excess'] > summary['sig_deficit'] * 2:
            f.write("NOTE: Heterozygote excess is more common than deficit.\n")
            f.write("This pattern may indicate:\n")
            f.write("  - Selection favoring heterozygotes\n")
            f.write("  - Negative assortative mating\n")
            f.write("  - Mixing of divergent populations\n\n")

        f.write("=" * 80 + "\n")
        f.write("Top 100 SNPs with Most Significant HWE Deviations\n")
        f.write("=" * 80 + "\n\n")

        if summary['top_deviations']:
            # Check if we have variant IDs
            has_variant_ids = summary['top_deviations'][0].get('variant_id') is not None

            if has_variant_ids:
                f.write(f"{'Variant ID':<20} {'RAD Locus':<15} {'CHR':<10} {'POS':<10} {'P_HWE':<12} {'Observed':<20} {'Expected':<20}\n")
                f.write("-" * 120 + "\n")
                for locus in summary['top_deviations']:
                    variant_id = locus.get('variant_id', 'N/A')
                    rad_locus = locus.get('rad_locus', 'N/A')
                    f.write(f"{variant_id:<20} {rad_locus:<15} {locus['chrom']:<10} {locus['pos']:<10} {locus['p_hwe']:<12.6e} ")
                    f.write(f"{locus['obs_counts']:<20} {locus['exp_counts']:<20}\n")
            else:
                f.write(f"{'CHR':<15} {'POS':<10} {'P_HWE':<12} {'Observed':<20} {'Expected':<20}\n")
                f.write("-" * 80 + "\n")
                for locus in summary['top_deviations']:
                    f.write(f"{locus['chrom']:<15} {locus['pos']:<10} {locus['p_hwe']:<12.6e} ")
                    f.write(f"{locus['obs_counts']:<20} {locus['exp_counts']:<20}\n")
        else:
            f.write("No significant deviations found.\n")

        f.write("\n")

        # Add RAD locus analysis section if available
        if summary['rad_locus_stats']:
            f.write("=" * 80 + "\n")
            f.write("RAD Locus Analysis\n")
            f.write("=" * 80 + "\n\n")

            # Summary statistics for RAD loci
            total_rad_loci = len(summary['rad_locus_stats'])
            rad_loci_with_sig = sum(1 for stats in summary['rad_locus_stats'].values() if stats['sig_snps'] > 0)
            rad_loci_all_sig = sum(1 for stats in summary['rad_locus_stats'].values()
                                   if stats['sig_snps'] == stats['total_snps'] and stats['total_snps'] > 0)

            f.write(f"Total RAD loci analyzed: {total_rad_loci}\n")
            f.write(f"RAD loci with at least one significant SNP: {rad_loci_with_sig} ({100.0 * rad_loci_with_sig / total_rad_loci:.1f}%)\n")
            f.write(f"RAD loci where ALL SNPs are significant: {rad_loci_all_sig} ({100.0 * rad_loci_all_sig / total_rad_loci:.1f}%)\n\n")

            # List ALL RAD loci where 100% of SNPs deviate from HWE
            loci_100_percent = sorted(
                [(locus, stats) for locus, stats in summary['rad_locus_stats'].items()
                 if stats['sig_snps'] == stats['total_snps'] and stats['total_snps'] > 0],
                key=lambda x: x[0]  # Sort by locus name
            )

            if loci_100_percent:
                f.write("-" * 80 + "\n")
                f.write(f"RAD loci where ALL SNPs deviate from HWE ({len(loci_100_percent)} loci):\n")
                f.write("-" * 80 + "\n\n")
                f.write(f"{'RAD Locus':<15} {'Total SNPs':<12} {'Sig SNPs':<12} {'Proportion':<12}\n")
                f.write("-" * 80 + "\n")

                for rad_locus, stats in loci_100_percent:
                    f.write(f"{rad_locus:<15} {stats['total_snps']:<12} {stats['sig_snps']:<12} {stats['prop_sig']:<12.3f}\n")

                f.write("\n")
                f.write("NOTE: RAD loci where all SNPs deviate from HWE may indicate:\n")
                f.write("  - Paralogs (duplicated loci)\n")
                f.write("  - Loci under strong selection\n")
                f.write("  - Null alleles at these loci\n")
                f.write("  - Consider filtering these entire loci from downstream analyses\n\n")

            # Distribution of deviation proportions
            f.write("-" * 80 + "\n")
            f.write("Distribution of HWE deviation rates across RAD loci:\n")
            f.write("-" * 80 + "\n\n")

            # Create bins for proportions: 0%, 1-25%, 26-50%, 51-75%, 76-99%, 100%
            bins = {
                '0% (no deviations)': 0,
                '1-25%': 0,
                '26-50%': 0,
                '51-75%': 0,
                '76-99%': 0,
                '100% (all SNPs deviate)': 0
            }

            for locus, stats in summary['rad_locus_stats'].items():
                prop = stats['prop_sig']
                if prop == 0:
                    bins['0% (no deviations)'] += 1
                elif prop <= 0.25:
                    bins['1-25%'] += 1
                elif prop <= 0.50:
                    bins['26-50%'] += 1
                elif prop <= 0.75:
                    bins['51-75%'] += 1
                elif prop < 1.0:
                    bins['76-99%'] += 1
                else:  # prop == 1.0
                    bins['100% (all SNPs deviate)'] += 1

            f.write(f"{'Deviation Rate':<30} {'Number of Loci':<15} {'Percentage':<15}\n")
            f.write("-" * 80 + "\n")
            for bin_name in ['0% (no deviations)', '1-25%', '26-50%', '51-75%', '76-99%', '100% (all SNPs deviate)']:
                count = bins[bin_name]
                pct = 100.0 * count / total_rad_loci if total_rad_loci > 0 else 0
                f.write(f"{bin_name:<30} {count:<15} {pct:.1f}%\n")

            f.write("\n")

        f.write("=" * 80 + "\n")
        f.write("Recommendations\n")
        f.write("=" * 80 + "\n\n")

        f.write("1. Review loci with extreme HWE deviations (p < 0.001)\n")
        f.write("2. Check for population substructure using PCA or STRUCTURE\n")
        f.write("3. Consider filtering loci that deviate significantly from HWE\n")
        f.write("4. If many loci show deficit, investigate allelic dropout\n")
        f.write("5. Cross-reference with whoa results for genotyping accuracy\n")
        f.write("6. Perform HWE tests within populations if pooling multiple groups\n\n")


def main():
    parser = argparse.ArgumentParser(
        description='Summarize Hardy-Weinberg equilibrium test results'
    )
    parser.add_argument('--input', required=True,
                        help='Input .hwe file from vcftools --hardy')
    parser.add_argument('--vcf', required=False,
                        help='Optional VCF file to extract variant IDs and RAD locus information')
    parser.add_argument('--output', required=True,
                        help='Output summary report file')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Significance threshold (default: 0.05)')

    args = parser.parse_args()

    # Build VCF position map if VCF file provided
    position_map = None
    if args.vcf:
        print(f"Reading VCF file to extract variant IDs: {args.vcf}")
        position_map = build_vcf_position_map(args.vcf)
        print(f"  Extracted {len(position_map)} variant positions")

    # Parse HWE results
    print(f"Reading Hardy-Weinberg results from: {args.input}")
    loci = parse_hardy_file(args.input, position_map=position_map)

    if not loci:
        print("ERROR: No loci found in input file", file=sys.stderr)
        sys.exit(1)

    # Generate summary
    print(f"Summarizing {len(loci)} loci...")
    summary = summarize_hwe(loci, alpha=args.alpha)

    # Write report
    print(f"Writing report to: {args.output}")
    write_report(summary, args.output, alpha=args.alpha, vcf_file=args.vcf)

    print("Hardy-Weinberg summary completed successfully!")
    print(f"  Total loci: {summary['total_loci']}")
    print(f"  Significant deviations: {summary['sig_deviations']} ({100.0 * summary['sig_deviations'] / summary['valid_loci']:.1f}%)")

    if summary['rad_locus_stats']:
        print(f"  RAD loci analyzed: {len(summary['rad_locus_stats'])}")


if __name__ == '__main__':
    main()
