#!/usr/bin/env python
"""Convert a biallelic VCF and indpopdata table to TreeMix allele-count input."""

import argparse
import gzip
import os
import re
import sys
from collections import OrderedDict

import pandas as pd
import vcfpy


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build TreeMix allele-count input from a VCF and indpopdata column."
    )
    parser.add_argument("--vcf", required=True, help="Input biallelic VCF/VCF.GZ")
    parser.add_argument("--indpopdata", required=True, help="indpopdata.txt table")
    parser.add_argument("--population-column", required=True, help="indpopdata column to group samples")
    parser.add_argument("--output", required=True, help="Output TreeMix .frq.gz file")
    parser.add_argument("--clust", required=True, help="Output three-column sample/sample/group file")
    parser.add_argument("--populations", required=True, help="Output population label mapping table")
    parser.add_argument("--positions", required=True, help="Output retained variant positions")
    parser.add_argument("--summary", required=True, help="Output conversion summary")
    parser.add_argument(
        "--min-called-populations",
        type=int,
        default=2,
        help="Minimum populations with at least one called allele at a SNP",
    )
    return parser.parse_args()


def get_snakemake_args():
    return argparse.Namespace(
        vcf=snakemake.input.vcf,
        indpopdata=snakemake.input.indpopdata,
        population_column=snakemake.params.population_column,
        output=snakemake.output.frq,
        clust=snakemake.output.clust,
        populations=snakemake.output.populations,
        positions=snakemake.output.positions,
        summary=snakemake.output.summary,
        min_called_populations=int(snakemake.params.min_called_populations),
    )


def setup_logging(log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    log_handle = open(log_path, "w")
    sys.stdout = log_handle
    sys.stderr = log_handle
    return log_handle


def sanitize_population_label(value):
    label = str(value).strip()
    if not label:
        raise ValueError("Population labels cannot be empty")
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label)


def load_indpopdata(path, population_column):
    data = pd.read_csv(path, sep="\t", dtype=str)
    if "Ind" not in data.columns:
        raise ValueError(f"{path} must contain an 'Ind' column")
    if population_column not in data.columns:
        available = ", ".join(data.columns)
        raise ValueError(
            f"Configured TreeMix population column '{population_column}' is not in "
            f"{path}. Available columns: {available}"
        )
    if data["Ind"].duplicated().any():
        dupes = data.loc[data["Ind"].duplicated(), "Ind"].tolist()
        raise ValueError(f"Duplicate individual IDs in {path}: {', '.join(dupes[:10])}")

    missing = data[population_column].isna() | (data[population_column].str.strip() == "")
    if missing.any():
        bad = data.loc[missing, "Ind"].tolist()
        raise ValueError(
            f"Missing TreeMix population assignments in column '{population_column}' "
            f"for: {', '.join(bad[:10])}"
        )

    return data.set_index("Ind")[population_column].to_dict()


def parse_gt(gt):
    if gt is None:
        return []
    alleles = re.split(r"[/|]", str(gt))
    return [allele for allele in alleles if allele != "."]


def main():
    args = get_snakemake_args() if "snakemake" in globals() else parse_args()
    log_handle = None
    if "snakemake" in globals() and len(snakemake.log) > 0:
        log_handle = setup_logging(snakemake.log[0])

    try:
        for path in [args.output, args.clust, args.populations, args.positions, args.summary]:
            os.makedirs(os.path.dirname(path), exist_ok=True)

        print(f"Reading indpopdata: {args.indpopdata}")
        sample_to_pop_raw = load_indpopdata(args.indpopdata, args.population_column)
        sample_to_pop = {}
        pop_label_map = OrderedDict()

        reader = vcfpy.Reader.from_path(args.vcf)
        samples = reader.header.samples.names
        print(f"Reading VCF: {args.vcf}")
        print(f"Found {len(samples)} samples in VCF")

        missing_samples = [sample for sample in samples if sample not in sample_to_pop_raw]
        if missing_samples:
            raise ValueError(
                "Samples found in VCF but missing from indpopdata: "
                + ", ".join(missing_samples[:20])
            )

        populations_seen = set()
        for sample in samples:
            raw_label = sample_to_pop_raw[sample]
            safe_label = sanitize_population_label(raw_label)
            if safe_label in pop_label_map and pop_label_map[safe_label] != raw_label:
                raise ValueError(
                    "TreeMix population labels collide after sanitizing whitespace/special "
                    f"characters: '{pop_label_map[safe_label]}' and '{raw_label}'"
                )
            pop_label_map[safe_label] = raw_label
            sample_to_pop[sample] = safe_label
            populations_seen.add(safe_label)

        populations = sorted(populations_seen)

        with open(args.clust, "w") as out:
            for sample in samples:
                out.write(f"{sample}\t{sample}\t{sample_to_pop[sample]}\n")

        with open(args.populations, "w") as out:
            out.write("treemix_population\tindpopdata_population\tpopulation_column\tn_samples\n")
            for pop in populations:
                n_samples = sum(1 for sample in samples if sample_to_pop[sample] == pop)
                out.write(f"{pop}\t{pop_label_map[pop]}\t{args.population_column}\t{n_samples}\n")

        total_records = 0
        skipped_non_biallelic = 0
        skipped_not_polymorphic = 0
        skipped_low_called_pops = 0
        retained_records = 0

        with gzip.open(args.output, "wt") as frq, open(args.positions, "w") as pos_out:
            frq.write(" ".join(populations) + "\n")
            pos_out.write("scaffold_pos\tscaffold\tpos\tcumulative_pos\n")
            scaffold_offsets = {}
            current_offset = 0

            for record in reader:
                total_records += 1
                if len(record.ALT) != 1:
                    skipped_non_biallelic += 1
                    continue

                counts = {pop: [0, 0] for pop in populations}
                for call in record.calls:
                    pop = sample_to_pop.get(call.sample)
                    if pop is None:
                        continue
                    for allele in parse_gt(call.data.get("GT")):
                        if allele == "0":
                            counts[pop][0] += 1
                        elif allele == "1":
                            counts[pop][1] += 1

                global_ref = sum(value[0] for value in counts.values())
                global_alt = sum(value[1] for value in counts.values())
                if global_ref == 0 or global_alt == 0:
                    skipped_not_polymorphic += 1
                    continue

                called_pops = sum(1 for value in counts.values() if sum(value) > 0)
                if called_pops < args.min_called_populations:
                    skipped_low_called_pops += 1
                    continue

                frq.write(" ".join(f"{counts[pop][0]},{counts[pop][1]}" for pop in populations) + "\n")

                if record.CHROM not in scaffold_offsets:
                    scaffold_offsets[record.CHROM] = current_offset
                cumulative_pos = scaffold_offsets[record.CHROM] + int(record.POS)
                current_offset = max(current_offset, cumulative_pos)
                pos_out.write(
                    f"{record.CHROM}:{record.POS}\t{record.CHROM}\t{record.POS}\t{cumulative_pos}\n"
                )
                retained_records += 1

        reader.close()

        with open(args.summary, "w") as out:
            out.write("metric\tvalue\n")
            out.write(f"vcf\t{args.vcf}\n")
            out.write(f"indpopdata\t{args.indpopdata}\n")
            out.write(f"population_column\t{args.population_column}\n")
            out.write(f"samples\t{len(samples)}\n")
            out.write(f"populations\t{len(populations)}\n")
            out.write(f"total_records\t{total_records}\n")
            out.write(f"retained_records\t{retained_records}\n")
            out.write(f"skipped_non_biallelic\t{skipped_non_biallelic}\n")
            out.write(f"skipped_not_polymorphic\t{skipped_not_polymorphic}\n")
            out.write(f"skipped_low_called_populations\t{skipped_low_called_pops}\n")

        print(f"Wrote TreeMix input: {args.output}")
        print(f"Retained {retained_records} SNPs across {len(populations)} populations")
    finally:
        if log_handle is not None:
            log_handle.close()


if __name__ == "__main__":
    main()
