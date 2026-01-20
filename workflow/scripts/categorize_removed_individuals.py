#!/usr/bin/env python3
"""
Categorize removed individuals by relationship type based on KING kinship values.

Reads KING table and samples-to-keep list, identifies removed individuals,
finds their maximum KING value, and categorizes them into:
- Clones: KING > 0.354 (duplicates/monozygotic twins)
- 1st-degree: 0.177 < KING <= 0.354 (parents/siblings)
- 2nd-degree: 0.0884 < KING <= 0.177 (half-siblings/grandparents)
- Other: KING <= 0.0884 (but above filtering threshold)
"""

import argparse
import sys
from collections import defaultdict

def parse_king_table(king_file):
    """Parse KING table and return dictionary of max KING values per individual."""
    max_king = defaultdict(float)
    
    with open(king_file, 'r') as f:
        # Skip header
        header = f.readline().strip().split('\t')
        
        # Find column indices
        # plink2 KING table uses IID1, IID2, KINSHIP columns
        try:
            id1_idx = header.index('IID1')
            id2_idx = header.index('IID2')
            kinship_idx = header.index('KINSHIP')
        except ValueError:
            # Try alternative column names
            try:
                id1_idx = header.index('ID1')
                id2_idx = header.index('ID2')
                kinship_idx = header.index('KINSHIP')
            except ValueError:
                try:
                    # Try FID columns (family ID)
                    id1_idx = header.index('#FID1')
                    id2_idx = header.index('#FID2')
                    kinship_idx = header.index('KINSHIP')
                except ValueError:
                    sys.stderr.write(f"Error: Could not find required columns in KING table. Found columns: {header}\n")
                    sys.exit(1)
        
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) <= max(id1_idx, id2_idx, kinship_idx):
                continue
            
            id1 = fields[id1_idx]
            id2 = fields[id2_idx]
            try:
                kinship = float(fields[kinship_idx])
            except (ValueError, IndexError):
                continue
            
            # Update max KING for both individuals
            if kinship > max_king[id1]:
                max_king[id1] = kinship
            if kinship > max_king[id2]:
                max_king[id2] = kinship
    
    return max_king

def read_samples_to_keep(samples_file):
    """Read list of samples to keep."""
    samples_to_keep = set()
    with open(samples_file, 'r') as f:
        for line in f:
            sample = line.strip()
            # Skip header line if present
            if sample and sample.upper() not in ['IID', 'ID', 'SAMPLE', 'INDIVIDUAL']:
                samples_to_keep.add(sample)
    return samples_to_keep

def categorize_relationship(king_value, thresholds):
    """Categorize relationship based on KING value."""
    if king_value > thresholds['clone']:
        return 'clone'
    elif king_value > thresholds['first_degree']:
        return '1st-degree'
    elif king_value > thresholds['second_degree']:
        return '2nd-degree'
    else:
        return 'other'

def main():
    parser = argparse.ArgumentParser(
        description='Categorize removed individuals by relationship type'
    )
    parser.add_argument(
        '--king-table',
        required=True,
        help='KING table file (TSV format)'
    )
    parser.add_argument(
        '--samples-to-keep',
        required=True,
        help='File with list of samples to keep (one per line)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output TSV file with removed individuals and categories'
    )
    parser.add_argument(
        '--clone-threshold',
        type=float,
        default=0.354,
        help='KING threshold for clones (default: 0.354)'
    )
    parser.add_argument(
        '--first-degree-threshold',
        type=float,
        default=0.177,
        help='KING threshold for 1st-degree relatives (default: 0.177)'
    )
    parser.add_argument(
        '--second-degree-threshold',
        type=float,
        default=0.0884,
        help='KING threshold for 2nd-degree relatives (default: 0.0884)'
    )
    
    args = parser.parse_args()
    
    # Read samples to keep
    samples_to_keep = read_samples_to_keep(args.samples_to_keep)
    
    # Parse KING table
    print(f"Reading KING table from {args.king_table}...", file=sys.stderr)
    max_king = parse_king_table(args.king_table)
    
    # Get all samples from KING table
    all_samples = set(max_king.keys())
    
    # Find removed individuals
    removed_samples = all_samples - samples_to_keep
    
    if not removed_samples:
        print("No removed individuals found.", file=sys.stderr)
        # Create empty output file
        with open(args.output, 'w') as f:
            f.write("individual\tmax_KING\tcategory\n")
        return
    
    # Categorize removed individuals
    thresholds = {
        'clone': args.clone_threshold,
        'first_degree': args.first_degree_threshold,
        'second_degree': args.second_degree_threshold
    }
    
    results = []
    for sample in removed_samples:
        king_value = max_king.get(sample, 0.0)
        category = categorize_relationship(king_value, thresholds)
        results.append((sample, king_value, category))
    
    # Sort by category, then by KING value (descending)
    category_order = {'clone': 0, '1st-degree': 1, '2nd-degree': 2, 'other': 3}
    results.sort(key=lambda x: (category_order.get(x[2], 99), -x[1]))
    
    # Write output
    print(f"Writing results to {args.output}...", file=sys.stderr)
    with open(args.output, 'w') as f:
        f.write("individual\tmax_KING\tcategory\n")
        for sample, king_value, category in results:
            f.write(f"{sample}\t{king_value:.6f}\t{category}\n")
    
    print(f"Found {len(removed_samples)} removed individuals.", file=sys.stderr)
    category_counts = defaultdict(int)
    for _, _, category in results:
        category_counts[category] += 1
    for category, count in sorted(category_counts.items()):
        print(f"  {category}: {count}", file=sys.stderr)

if __name__ == '__main__':
    main()

