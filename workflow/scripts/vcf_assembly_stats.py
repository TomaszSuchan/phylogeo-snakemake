#!/usr/bin/env python3
"""
Comprehensive VCF assembly statistics for RADseq data.
Analyzes variants per RAD fragment, missing data, and other quality metrics.
"""

import argparse
import gzip
import re
from collections import defaultdict, Counter
import sys
import numpy as np


def parse_vcf_line(line):
    """Parse a VCF data line and extract relevant information."""
    fields = line.strip().split('\t')
    if len(fields) < 10:
        return None

    chrom = fields[0]
    pos = fields[1]
    variant_id = fields[2]
    ref = fields[3]
    alt = fields[4]
    qual = fields[5]
    filter_field = fields[6]
    info = fields[7]
    format_field = fields[8]
    samples = fields[9:]

    return {
        'chrom': chrom,
        'pos': pos,
        'id': variant_id,
        'ref': ref,
        'alt': alt,
        'qual': qual,
        'filter': filter_field,
        'info': info,
        'format': format_field,
        'samples': samples
    }


def extract_rad_fragment(variant_id, pattern):
    """Extract RAD fragment ID from variant ID using regex pattern."""
    match = re.match(pattern, variant_id)
    if match:
        return match.group(1)
    return None


def calculate_missing_data(samples):
    """Calculate missing data statistics for samples."""
    total = len(samples)
    missing = sum(1 for s in samples if s.startswith('./.') or s.startswith('.|.'))
    return missing, total, missing / total if total > 0 else 0


def parse_info_field(info):
    """Parse INFO field and extract key statistics."""
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def analyze_vcf(vcf_file, id_pattern):
    """Analyze VCF file and generate comprehensive statistics."""

    # Initialize counters
    total_variants = 0
    variants_per_fragment = defaultdict(list)
    chromosomes = set()
    quality_scores = []
    filter_counts = Counter()
    missing_data_per_variant = []

    # Sample-level stats
    sample_names = []
    sample_missing = defaultdict(int)
    sample_total = defaultdict(int)

    # Info field stats
    ns_values = []  # Number of samples with data
    dp_values = []  # INFO DP field values (for comparison)
    info_dp_description = None  # Store the DP field description from header

    # Per-variant depth statistics (calculated from FORMAT DP)
    variant_mean_depths = []  # Mean depth per variant across samples
    variant_depth_std = []  # Standard deviation of depth per variant

    # Per-sample depth tracking
    sample_depths = defaultdict(list)  # Track depth per sample

    # Heterozygosity tracking (for over-splitting/over-lumping detection)
    sample_het_sites = defaultdict(int)  # Heterozygous sites per sample
    sample_hom_sites = defaultdict(int)  # Homozygous sites per sample

    # Variant types
    snp_types = Counter()  # transitions vs transversions

    # Allele counts
    biallelic_count = 0
    multiallelic_count = 0
    allele_counts = Counter()  # distribution of number of alleles

    # Open file (handles both gzipped and plain text)
    opener = gzip.open if vcf_file.endswith('.gz') else open

    with opener(vcf_file, 'rt') as f:
        for line in f:
            line = line.strip()

            # Parse header to get sample names and INFO field descriptions
            if line.startswith('#CHROM'):
                fields = line.split('\t')
                if len(fields) > 9:
                    sample_names = fields[9:]
                continue

            # Capture INFO DP field description
            if line.startswith('##INFO=<ID=DP,'):
                # Extract description field
                desc_match = re.search(r'Description="([^"]+)"', line)
                if desc_match:
                    info_dp_description = desc_match.group(1)
                continue

            # Skip other header lines
            if line.startswith('#'):
                continue

            # Parse variant line
            variant = parse_vcf_line(line)
            if not variant:
                continue

            total_variants += 1

            # Extract RAD fragment
            fragment_id = extract_rad_fragment(variant['id'], id_pattern)
            if fragment_id:
                variants_per_fragment[fragment_id].append(variant['id'])

            # Track chromosomes
            chromosomes.add(variant['chrom'])

            # Quality scores
            if variant['qual'] != '.' and variant['qual'] != 'PASS':
                try:
                    quality_scores.append(float(variant['qual']))
                except ValueError:
                    pass

            # Filter field
            filter_counts[variant['filter']] += 1

            # Missing data
            missing, total, missing_prop = calculate_missing_data(variant['samples'])
            missing_data_per_variant.append(missing_prop)

            # Collect depths for this variant to calculate mean and SD
            variant_depths_this_site = []

            # Per-sample missing data and genotype analysis
            for i, sample_gt in enumerate(variant['samples']):
                sample_name = sample_names[i] if i < len(sample_names) else f"Sample_{i}"
                sample_total[sample_name] += 1

                # Parse genotype field (FORMAT:SAMPLE)
                gt_fields = sample_gt.split(':')
                genotype = gt_fields[0]

                if genotype.startswith('./.') or genotype.startswith('.|.'):
                    sample_missing[sample_name] += 1
                else:
                    # Track heterozygosity
                    alleles = genotype.replace('|', '/').split('/')
                    if len(alleles) == 2:
                        if alleles[0] != alleles[1]:
                            sample_het_sites[sample_name] += 1
                        else:
                            sample_hom_sites[sample_name] += 1

                    # Extract depth if available (DP field in genotype)
                    if len(gt_fields) > 1 and variant['format']:
                        format_keys = variant['format'].split(':')
                        if 'DP' in format_keys:
                            dp_idx = format_keys.index('DP')
                            if dp_idx < len(gt_fields):
                                try:
                                    depth = int(gt_fields[dp_idx])
                                    sample_depths[sample_name].append(depth)
                                    variant_depths_this_site.append(depth)
                                except (ValueError, IndexError):
                                    pass

            # Calculate mean and standard deviation of depth for this variant
            if variant_depths_this_site:
                variant_mean_depths.append(np.mean(variant_depths_this_site))
                if len(variant_depths_this_site) > 1:
                    variant_depth_std.append(np.std(variant_depths_this_site, ddof=1))
                else:
                    variant_depth_std.append(0.0)

            # Parse INFO field
            info_dict = parse_info_field(variant['info'])
            if 'NS' in info_dict:
                try:
                    ns_values.append(int(info_dict['NS']))
                except ValueError:
                    pass
            # NOTE: DP in INFO field is typically the sum of depths across all samples,
            # not the mean. This is why total DP values are higher than per-sample depths.
            if 'DP' in info_dict:
                try:
                    dp_values.append(int(info_dict['DP']))
                except ValueError:
                    pass

            # SNP type (transition vs transversion)
            ref = variant['ref']
            alt = variant['alt']
            if len(ref) == 1 and len(alt) == 1:
                snp_pair = tuple(sorted([ref, alt]))
                if snp_pair in [('A', 'G'), ('C', 'T')]:
                    snp_types['transition'] += 1
                elif all(b in 'ACGT' for b in snp_pair):
                    snp_types['transversion'] += 1

            # Count biallelic vs multiallelic variants
            # ALT field can contain comma-separated alternate alleles
            alt_alleles = variant['alt'].split(',')
            num_alleles = 1 + len(alt_alleles)  # 1 for ref + number of alt alleles
            allele_counts[num_alleles] += 1

            if num_alleles == 2:
                biallelic_count += 1
            elif num_alleles > 2:
                multiallelic_count += 1

    return {
        'total_variants': total_variants,
        'variants_per_fragment': dict(variants_per_fragment),
        'chromosomes': chromosomes,
        'quality_scores': quality_scores,
        'filter_counts': filter_counts,
        'missing_data_per_variant': missing_data_per_variant,
        'sample_names': sample_names,
        'sample_missing': sample_missing,
        'sample_total': sample_total,
        'sample_depths': dict(sample_depths),
        'sample_het_sites': dict(sample_het_sites),
        'sample_hom_sites': dict(sample_hom_sites),
        'ns_values': ns_values,
        'dp_values': dp_values,
        'info_dp_description': info_dp_description,
        'variant_mean_depths': variant_mean_depths,
        'variant_depth_std': variant_depth_std,
        'snp_types': snp_types,
        'biallelic_count': biallelic_count,
        'multiallelic_count': multiallelic_count,
        'allele_counts': allele_counts
    }


def write_report(stats, output_file, vcf_file=None):
    """Write comprehensive statistics report."""

    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("VCF ASSEMBLY STATISTICS REPORT\n")
        f.write("=" * 80 + "\n\n")

        if vcf_file:
            f.write(f"Input VCF: {vcf_file}\n\n")

        # Overall variant counts
        f.write("VARIANT COUNTS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total variants: {stats['total_variants']:,}\n")
        f.write(f"Number of RAD fragments: {len(stats['variants_per_fragment']):,}\n")
        f.write(f"Number of chromosomes/contigs: {len(stats['chromosomes']):,}\n\n")

        # Biallelic vs multiallelic
        f.write("BIALLELIC VS MULTIALLELIC VARIANTS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Biallelic variants (2 alleles): {stats['biallelic_count']:,} ({100 * stats['biallelic_count'] / stats['total_variants']:.2f}%)\n")
        f.write(f"Multiallelic variants (>2 alleles): {stats['multiallelic_count']:,} ({100 * stats['multiallelic_count'] / stats['total_variants']:.2f}%)\n")
        f.write("\nAllele count distribution:\n")
        for num_alleles in sorted(stats['allele_counts'].keys()):
            count = stats['allele_counts'][num_alleles]
            f.write(f"  {num_alleles} alleles: {count:,} variants ({100 * count / stats['total_variants']:.2f}%)\n")
        f.write("\n")

        # Variants per fragment
        f.write("VARIANTS PER RAD FRAGMENT\n")
        f.write("-" * 80 + "\n")
        variants_per_frag_counts = [len(v) for v in stats['variants_per_fragment'].values()]
        if variants_per_frag_counts:
            f.write(f"Mean variants per fragment: {sum(variants_per_frag_counts) / len(variants_per_frag_counts):.2f}\n")
            f.write(f"Min variants per fragment: {min(variants_per_frag_counts)}\n")
            f.write(f"Max variants per fragment: {max(variants_per_frag_counts)}\n")

            # Distribution
            frag_dist = Counter(variants_per_frag_counts)
            f.write("\nDistribution of variants per fragment:\n")
            for count in sorted(frag_dist.keys()):
                f.write(f"  {count} variant(s): {frag_dist[count]:,} fragments\n")
        f.write("\n")

        # Quality scores
        if stats['quality_scores']:
            f.write("QUALITY SCORES\n")
            f.write("-" * 80 + "\n")
            f.write(f"Mean quality: {sum(stats['quality_scores']) / len(stats['quality_scores']):.2f}\n")
            f.write(f"Min quality: {min(stats['quality_scores']):.2f}\n")
            f.write(f"Max quality: {max(stats['quality_scores']):.2f}\n\n")

        # Filter field
        f.write("FILTER FIELD STATISTICS\n")
        f.write("-" * 80 + "\n")
        for filter_val, count in stats['filter_counts'].most_common():
            f.write(f"{filter_val}: {count:,} variants ({100 * count / stats['total_variants']:.2f}%)\n")
        f.write("\n")

        # Missing data
        f.write("MISSING DATA STATISTICS\n")
        f.write("-" * 80 + "\n")
        if stats['missing_data_per_variant']:
            mean_missing = sum(stats['missing_data_per_variant']) / len(stats['missing_data_per_variant'])
            f.write(f"Mean missing data per variant: {100 * mean_missing:.2f}%\n")
            f.write(f"Min missing data: {100 * min(stats['missing_data_per_variant']):.2f}%\n")
            f.write(f"Max missing data: {100 * max(stats['missing_data_per_variant']):.2f}%\n\n")

            # Missing data distribution
            missing_bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            missing_dist = [0] * (len(missing_bins) - 1)
            for missing in stats['missing_data_per_variant']:
                for i in range(len(missing_bins) - 1):
                    if missing_bins[i] <= missing < missing_bins[i + 1]:
                        missing_dist[i] += 1
                        break
                if missing == 1.0:
                    missing_dist[-1] += 1

            f.write("Missing data distribution:\n")
            for i in range(len(missing_bins) - 1):
                pct_variants = 100 * missing_dist[i] / len(stats['missing_data_per_variant'])
                f.write(f"  {100 * missing_bins[i]:.0f}%-{100 * missing_bins[i+1]:.0f}%: {missing_dist[i]:,} variants ({pct_variants:.2f}%)\n")
        f.write("\n")

        # Per-sample missing data
        f.write("PER-SAMPLE MISSING DATA\n")
        f.write("-" * 80 + "\n")
        for sample in sorted(stats['sample_names']):
            if sample in stats['sample_total'] and stats['sample_total'][sample] > 0:
                missing_prop = stats['sample_missing'][sample] / stats['sample_total'][sample]
                f.write(f"{sample}: {100 * missing_prop:.2f}% ({stats['sample_missing'][sample]:,}/{stats['sample_total'][sample]:,})\n")
        f.write("\n")

        # NS (number of samples) stats
        if stats['ns_values']:
            f.write("NUMBER OF SAMPLES WITH DATA (NS)\n")
            f.write("-" * 80 + "\n")
            f.write(f"Mean NS: {sum(stats['ns_values']) / len(stats['ns_values']):.2f}\n")
            f.write(f"Min NS: {min(stats['ns_values'])}\n")
            f.write(f"Max NS: {max(stats['ns_values'])}\n\n")

        # Mean depth per variant (calculated from FORMAT DP)
        if stats['variant_mean_depths']:
            f.write("MEAN DEPTH PER VARIANT (calculated from FORMAT DP)\n")
            f.write("-" * 80 + "\n")
            f.write(f"Overall mean depth: {np.mean(stats['variant_mean_depths']):.2f}x\n")
            f.write(f"Overall SD: {np.std(stats['variant_mean_depths']):.2f}x\n")
            f.write(f"Min mean depth: {min(stats['variant_mean_depths']):.2f}x\n")
            f.write(f"Max mean depth: {max(stats['variant_mean_depths']):.2f}x\n")
            if stats['variant_depth_std']:
                f.write(f"\nMean within-variant SD: {np.mean(stats['variant_depth_std']):.2f}x\n")
            f.write("\n")

        # INFO DP field (for comparison if available)
        if stats['dp_values']:
            f.write("INFO DP FIELD (for comparison)\n")
            f.write("-" * 80 + "\n")
            if stats['info_dp_description']:
                f.write(f"VCF header description: {stats['info_dp_description']}\n\n")
            f.write(f"Mean INFO DP: {np.mean(stats['dp_values']):.2f}\n")
            f.write(f"Min INFO DP: {min(stats['dp_values'])}\n")
            f.write(f"Max INFO DP: {max(stats['dp_values'])}\n")
            # Calculate ratio to help understand if it's sum or mean
            if stats['variant_mean_depths'] and stats['ns_values']:
                mean_ratio = np.mean(stats['dp_values']) / np.mean(stats['variant_mean_depths'])
                mean_ns = np.mean(stats['ns_values'])
                f.write(f"\nRatio INFO_DP / FORMAT_DP_mean: {mean_ratio:.2f}\n")
                f.write(f"Mean number of samples with data: {mean_ns:.2f}\n")
                if abs(mean_ratio - mean_ns) < 1.0:
                    f.write("INFO DP appears to be the SUM of depths across samples\n")
                elif abs(mean_ratio - 1.0) < 0.2:
                    f.write("INFO DP appears to be the MEAN depth across samples\n")
            f.write("\n")

        # SNP types
        if stats['snp_types']:
            f.write("SNP TYPES\n")
            f.write("-" * 80 + "\n")
            total_typed = sum(stats['snp_types'].values())
            for snp_type, count in stats['snp_types'].most_common():
                f.write(f"{snp_type}: {count:,} ({100 * count / total_typed:.2f}%)\n")
            if stats['snp_types']['transition'] > 0 and stats['snp_types']['transversion'] > 0:
                ti_tv_ratio = stats['snp_types']['transition'] / stats['snp_types']['transversion']
                f.write(f"\nTransition/Transversion ratio: {ti_tv_ratio:.2f}\n")
        f.write("\n")

        # Per-sample depth statistics
        if stats['sample_depths']:
            f.write("PER-SAMPLE DEPTH STATISTICS\n")
            f.write("-" * 80 + "\n")
            for sample in sorted(stats['sample_names']):
                if sample in stats['sample_depths'] and stats['sample_depths'][sample]:
                    depths = stats['sample_depths'][sample]
                    mean_depth = sum(depths) / len(depths)
                    min_depth = min(depths)
                    max_depth = max(depths)
                    f.write(f"{sample}: mean={mean_depth:.2f}x, min={min_depth}x, max={max_depth}x\n")
            f.write("\n")

        # Observed heterozygosity per sample
        if stats['sample_het_sites'] or stats['sample_hom_sites']:
            f.write("OBSERVED HETEROZYGOSITY PER SAMPLE\n")
            f.write("-" * 80 + "\n")
            for sample in sorted(stats['sample_names']):
                het_count = stats['sample_het_sites'].get(sample, 0)
                hom_count = stats['sample_hom_sites'].get(sample, 0)
                total_called = het_count + hom_count
                if total_called > 0:
                    obs_het = het_count / total_called
                    f.write(f"{sample}: {obs_het:.4f} ({het_count:,}/{total_called:,} sites)\n")
            f.write("\n")

        # RAD-seq Assembly Sanity Checks (Literature-based)
        f.write("=" * 80 + "\n")
        f.write("RAD-SEQ ASSEMBLY SANITY CHECKS\n")
        f.write("=" * 80 + "\n")

        recommendations = []

        # Check 1: Coverage depth (Paris et al. 2017, Rochette et al. 2019)
        if stats['dp_values']:
            mean_dp = sum(stats['dp_values']) / len(stats['dp_values'])
            f.write(f"[1] COVERAGE DEPTH CHECK\n")
            f.write(f"    Mean depth: {mean_dp:.2f}x\n")
            if mean_dp < 10:
                f.write(f"    ⚠️  WARNING: Coverage <10x may lead to high genotyping error\n")
            elif mean_dp < 25:
                f.write(f"    ⚠️  CAUTION: Coverage <25x; >25x recommended for high accuracy\n")
            else:
                f.write(f"    ✓ PASS: Coverage ≥25x meets recommended threshold\n")
            f.write("\n")

        # Check 2: Missing data (Eaton 2014, Rivera-Colón & Catchen 2022)
        if stats['missing_data_per_variant']:
            mean_missing = sum(stats['missing_data_per_variant']) / len(stats['missing_data_per_variant'])
            f.write(f"[2] MISSING DATA CHECK\n")
            f.write(f"    Mean missing data per variant: {100 * mean_missing:.2f}%\n")
            high_missing_variants = sum(1 for m in stats['missing_data_per_variant'] if m > 0.5)
            pct_high_missing = 100 * high_missing_variants / len(stats['missing_data_per_variant'])
            f.write(f"    Variants with >50% missing: {pct_high_missing:.2f}%\n")
            if mean_missing > 0.5:
                f.write(f"    ⚠️  WARNING: High missing data may indicate poor assembly parameters\n")
            elif pct_high_missing > 30:
                f.write(f"    ⚠️  CAUTION: Many high-missingness variants; consider filtering\n")
            else:
                f.write(f"    ✓ PASS: Missing data levels acceptable\n")
            f.write("    NOTE: Avoid over-filtering on missingness (min_samples_locus)\n")
            f.write("\n")

        # Check 3: Ti/Tv ratio (DePristo et al. 2011, natural genomes)
        if stats['snp_types']['transition'] > 0 and stats['snp_types']['transversion'] > 0:
            ti_tv_ratio = stats['snp_types']['transition'] / stats['snp_types']['transversion']
            f.write(f"[3] TRANSITION/TRANSVERSION RATIO CHECK\n")
            f.write(f"    Ti/Tv ratio: {ti_tv_ratio:.2f}\n")
            if ti_tv_ratio < 1.5:
                f.write(f"    ⚠️  CAUTION: Ti/Tv <1.5 may indicate quality issues\n")
                f.write(f"    Expected range: 2.0-2.5 for most organisms\n")
            elif ti_tv_ratio > 3.5:
                f.write(f"    ⚠️  CAUTION: Ti/Tv >3.5 is unusually high\n")
            else:
                f.write(f"    ✓ PASS: Ti/Tv ratio within expected range (1.5-3.5)\n")
            f.write("\n")

        # Check 4: Observed heterozygosity (Eaton 2014, ipyrad guidelines)
        if stats['sample_het_sites'] and stats['sample_hom_sites']:
            all_het = sum(stats['sample_het_sites'].values())
            all_hom = sum(stats['sample_hom_sites'].values())
            if (all_het + all_hom) > 0:
                global_obs_het = all_het / (all_het + all_hom)
                f.write(f"[4] HETEROZYGOSITY CHECK (clustering threshold assessment)\n")
                f.write(f"    Global observed heterozygosity: {global_obs_het:.4f}\n")
                if global_obs_het < 0.001:
                    f.write(f"    ⚠️  WARNING: Het <0.001 suggests over-splitting (threshold too high)\n")
                    recommendations.append("Consider decreasing clustering threshold")
                elif global_obs_het > 0.05:
                    f.write(f"    ⚠️  WARNING: Het >0.05 suggests over-lumping (threshold too low)\n")
                    recommendations.append("Consider increasing clustering threshold")
                else:
                    f.write(f"    ✓ PASS: Heterozygosity within expected range (~0.001-0.05)\n")
                f.write("\n")

        # Check 5: Variants per RAD fragment
        if stats['variants_per_fragment']:
            variants_per_frag_counts = [len(v) for v in stats['variants_per_fragment'].values()]
            mean_vars = sum(variants_per_frag_counts) / len(variants_per_frag_counts)
            max_vars = max(variants_per_frag_counts)
            high_var_frags = sum(1 for v in variants_per_frag_counts if v > 10)

            f.write(f"[5] VARIANTS PER FRAGMENT CHECK\n")
            f.write(f"    Mean variants per fragment: {mean_vars:.2f}\n")
            f.write(f"    Max variants per fragment: {max_vars}\n")
            if max_vars > 20:
                f.write(f"    ⚠️  CAUTION: Fragments with many variants may indicate paralogs\n")
                recommendations.append("Consider paralog detection/filtering")
            if mean_vars > 5:
                f.write(f"    ⚠️  CAUTION: Mean >5 variants/fragment suggests over-lumping\n")
            else:
                f.write(f"    ✓ PASS: Variants per fragment distribution looks reasonable\n")
            f.write("\n")

        # Check 6: Sample-level quality uniformity
        if stats['sample_missing'] and stats['sample_total']:
            sample_missing_props = []
            for sample in stats['sample_names']:
                if stats['sample_total'][sample] > 0:
                    missing_prop = stats['sample_missing'][sample] / stats['sample_total'][sample]
                    sample_missing_props.append(missing_prop)

            if sample_missing_props:
                max_sample_missing = max(sample_missing_props)
                min_sample_missing = min(sample_missing_props)
                range_missing = max_sample_missing - min_sample_missing

                f.write(f"[6] SAMPLE QUALITY UNIFORMITY CHECK\n")
                f.write(f"    Range of missing data across samples: {100*range_missing:.2f}%\n")
                if max_sample_missing > 0.75:
                    poor_samples = [s for s in stats['sample_names']
                                   if stats['sample_total'][s] > 0 and
                                   stats['sample_missing'][s] / stats['sample_total'][s] > 0.75]
                    f.write(f"    ⚠️  WARNING: {len(poor_samples)} samples with >75% missing data\n")
                    f.write(f"    Poor samples (for ipyrad removal, paste as-is):\n")
                    f.write(f"    {' '.join(poor_samples)}\n")
                    recommendations.append("Consider removing low-quality samples")
                elif range_missing > 0.3:
                    f.write(f"    ⚠️  CAUTION: High variation in sample quality (>30% range)\n")
                else:
                    f.write(f"    ✓ PASS: Sample quality is relatively uniform\n")
                f.write("\n")



def main():
    parser = argparse.ArgumentParser(
        description='Generate comprehensive VCF assembly statistics for RADseq data'
    )
    parser.add_argument('--vcf', required=True, help='Input VCF file (can be gzipped)')
    parser.add_argument('--out', required=True, help='Output statistics report file')
    parser.add_argument(
        '--id-pattern',
        default=r'loc(\d+)_',
        help='Regex pattern to extract RAD fragment ID from variant ID (default: "loc(\\d+)_")'
    )

    args = parser.parse_args()

    print(f"Analyzing VCF file: {args.vcf}", file=sys.stderr)
    stats = analyze_vcf(args.vcf, args.id_pattern)

    print(f"Writing report to: {args.out}", file=sys.stderr)
    write_report(stats, args.out, vcf_file=args.vcf)

    print("Analysis complete!", file=sys.stderr)


if __name__ == '__main__':
    main()
