#!/usr/bin/env python3
"""
List individuals removed by PLINK 2 --king-cutoff, and why.

PLINK removes one member of each pair whose KING kinship exceeds the cutoff.
When several relatives remain, it removes the sample with more such relatives
(the partner of the first degree-1 sample, otherwise the first sample of
maximal degree, in VCF sample order). The two-column table records that
sample and the relative or relatives that were still above the cutoff when
it was removed.
"""

import argparse
import sys
from collections import defaultdict


def parse_king_table(king_file):
    """Return (sample_order, list of (id1, id2, kinship))."""
    pairs = []
    order = []
    seen = set()
    with open(king_file, "r") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if not header or header == [""]:
            return order, pairs
        try:
            id1_idx = header.index("IID1")
            id2_idx = header.index("IID2")
            kinship_idx = header.index("KINSHIP")
        except ValueError:
            try:
                id1_idx = header.index("ID1")
                id2_idx = header.index("ID2")
                kinship_idx = header.index("KINSHIP")
            except ValueError:
                sys.stderr.write(
                    "Error: Could not find IID1/IID2/KINSHIP in KING table. "
                    f"Found columns: {header}\n"
                )
                sys.exit(1)
        first = True
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(id1_idx, id2_idx, kinship_idx):
                continue
            id1 = fields[id1_idx]
            id2 = fields[id2_idx]
            try:
                kinship = float(fields[kinship_idx])
            except ValueError:
                continue
            pairs.append((id1, id2, kinship))
            # PLINK writes the higher sample index as IID1 and the lower as IID2,
            # so the first IID2 is sample 0 and each new IID1 is the next sample.
            if first:
                order.append(id2)
                seen.add(id2)
                first = False
            if id1 not in seen:
                order.append(id1)
                seen.add(id1)
    return order, pairs


def read_samples_to_keep(samples_file):
    samples_to_keep = set()
    with open(samples_file, "r") as handle:
        for line in handle:
            sample = line.strip()
            if sample and sample.upper() not in ["IID", "ID", "SAMPLE", "INDIVIDUAL"]:
                samples_to_keep.add(sample)
    return samples_to_keep


def relationship_category(king_value, thresholds):
    """Category labels expected by the removed-individuals barplot."""
    if king_value > thresholds["clone"]:
        return "clone"
    if king_value > thresholds["first_degree"]:
        return "1st-degree"
    if king_value > thresholds["second_degree"]:
        return "2nd-degree"
    return "other"


def relationship_phrase(king_value, thresholds):
    if king_value > thresholds["clone"]:
        return "duplicate"
    if king_value > thresholds["first_degree"]:
        return "1st-degree relative"
    if king_value > thresholds["second_degree"]:
        return "2nd-degree relative"
    return "relative above the KING cutoff"


def removal_reason(partners, thresholds):
    """partners: list of (sample, kinship), closest first."""
    parts = []
    for sample, kinship in partners:
        kind = relationship_phrase(kinship, thresholds)
        parts.append(f"{kind} of {sample} (KING kinship {kinship:.4f})")
    return "; ".join(parts)


def prune_related(sample_order, pairs, threshold):
    """Reproduce PLINK 2 KinshipPruneDestructive. Return removal records."""
    neighbors = defaultdict(set)
    kinship = {}
    for id1, id2, value in pairs:
        if value > threshold:
            neighbors[id1].add(id2)
            neighbors[id2].add(id1)
            kinship[(id1, id2)] = value
            kinship[(id2, id1)] = value
    active = {sample for sample in sample_order if neighbors[sample]}
    removed = []
    while active:
        degree_one = [
            sample for sample in sample_order
            if sample in active and len(neighbors[sample]) == 1
        ]
        if degree_one:
            prune = next(iter(neighbors[degree_one[0]]))
        else:
            max_degree = max(len(neighbors[sample]) for sample in active)
            prune = next(
                sample for sample in sample_order
                if sample in active and len(neighbors[sample]) == max_degree
            )
        partners = [
            (sample, kinship[(prune, sample)])
            for sample in sample_order
            if sample in neighbors[prune]
        ]
        partners.sort(key=lambda item: -item[1])
        removed.append((prune, partners))
        for partner, _value in partners:
            neighbors[partner].discard(prune)
            if not neighbors[partner]:
                active.discard(partner)
        neighbors[prune].clear()
        active.discard(prune)
    return removed


def main():
    parser = argparse.ArgumentParser(
        description="List individuals removed by KING cutoff and why"
    )
    parser.add_argument("--king-table", required=True)
    parser.add_argument("--samples-to-keep", required=True)
    parser.add_argument(
        "--output",
        required=True,
        help="Three-column table used by the removed-individuals barplot",
    )
    parser.add_argument(
        "--samples-removed",
        required=True,
        help="Two-column table: sample and the relative that caused the exclusion",
    )
    parser.add_argument("--king-threshold", type=float, default=0.0884)
    parser.add_argument("--clone-threshold", type=float, default=0.354)
    parser.add_argument("--first-degree-threshold", type=float, default=0.177)
    parser.add_argument("--second-degree-threshold", type=float, default=0.0884)
    args = parser.parse_args()

    thresholds = {
        "clone": args.clone_threshold,
        "first_degree": args.first_degree_threshold,
        "second_degree": args.second_degree_threshold,
    }
    samples_to_keep = read_samples_to_keep(args.samples_to_keep)
    sample_order, pairs = parse_king_table(args.king_table)
    removed = prune_related(sample_order, pairs, args.king_threshold)

    expected = set(sample_order) - samples_to_keep
    got = {sample for sample, _partners in removed}
    if got != expected:
        sys.stderr.write(
            "Error: KING cutoff reconstruction does not match samples_to_keep "
            f"({len(got)} reconstructed, {len(expected)} absent from the keep list).\n"
        )
        sys.exit(1)

    with open(args.samples_removed, "w") as handle:
        handle.write("sample\treason\n")
        for sample, partners in removed:
            handle.write(f"{sample}\t{removal_reason(partners, thresholds)}\n")

    with open(args.output, "w") as handle:
        handle.write("individual\tmax_KING\tcategory\n")
        rows = []
        for sample, partners in removed:
            max_kinship = partners[0][1]
            category = relationship_category(max_kinship, thresholds)
            rows.append((sample, max_kinship, category))
        category_order = {"clone": 0, "1st-degree": 1, "2nd-degree": 2, "other": 3}
        rows.sort(key=lambda row: (category_order.get(row[2], 99), -row[1], row[0]))
        for sample, max_kinship, category in rows:
            handle.write(f"{sample}\t{max_kinship:.6f}\t{category}\n")

    print(f"Removed {len(removed)} individuals.", file=sys.stderr)


if __name__ == "__main__":
    main()
