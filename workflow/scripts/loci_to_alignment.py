#!/usr/bin/env python3
"""Build a concatenated PHYLIP/FASTA alignment from ipyrad .loci for retained samples."""

import shutil
import sys
import tempfile
from pathlib import Path

# Snakemake wraps scripts and may copy them to a temp path; prefer scriptdir.
try:
    _SCRIPT_DIR = Path(snakemake.scriptdir)
except NameError:
    _SCRIPT_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(_SCRIPT_DIR / "extract-invariant-vcf"))
from extract_invariant_vcf import iter_loci_file  # noqa: E402


def _read_keep_samples(path):
    samples = []
    seen = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            sample = line.strip()
            if not sample or sample in seen:
                continue
            seen.add(sample)
            samples.append(sample)
    if not samples:
        raise SystemExit("No sample names found in {}".format(path))
    return samples


def _write_fasta(path, sample, seq_path):
    """Write single-line FASTA (required by rapidNJ -i fa)."""
    with open(path, "a", encoding="utf-8") as out, open(seq_path, encoding="utf-8") as seq_fh:
        out.write(">{}\n".format(sample))
        shutil.copyfileobj(seq_fh, out)
        out.write("\n")


def main(loci_path, samples_path, phy_path, fasta_path, summary_path):
    keep = _read_keep_samples(samples_path)
    keep_set = set(keep)
    tmpdir = Path(tempfile.mkdtemp(prefix="loci_to_alignment_"))
    try:
        seq_paths = {sample: tmpdir / "{}.seq".format(i) for i, sample in enumerate(keep)}
        handles = {sample: open(path, "w", encoding="utf-8") for sample, path in seq_paths.items()}
        loci_present = {sample: 0 for sample in keep}
        n_loci = 0
        n_sites = 0
        try:
            for locus in iter_loci_file(loci_path):
                sequences = locus["sequences"]
                present = [sample for sample in keep if sample in sequences]
                if not present:
                    continue
                lengths = {len(sequences[sample]) for sample in present}
                locus_len = max(lengths)
                fill = "N" * locus_len
                for sample in keep:
                    seq = sequences.get(sample)
                    if seq is None:
                        handles[sample].write(fill)
                        continue
                    loci_present[sample] += 1
                    if len(seq) == locus_len:
                        handles[sample].write(seq)
                    elif len(seq) < locus_len:
                        handles[sample].write(seq + ("N" * (locus_len - len(seq))))
                    else:
                        handles[sample].write(seq[:locus_len])
                n_loci += 1
                n_sites += locus_len
        finally:
            for handle in handles.values():
                handle.close()

        missing = [sample for sample, count in loci_present.items() if count == 0]
        if missing:
            preview = ", ".join(missing[:10])
            extra = "" if len(missing) <= 10 else " (+{} more)".format(len(missing) - 10)
            raise SystemExit(
                "{} sample(s) from the keep list never appear in {}: {}{}".format(
                    len(missing), loci_path, preview, extra
                )
            )
        if n_loci == 0 or n_sites == 0:
            raise SystemExit(
                "No loci retained after filtering {} to {} samples".format(loci_path, len(keep))
            )

        phy = Path(phy_path)
        fasta = Path(fasta_path)
        phy.parent.mkdir(parents=True, exist_ok=True)
        fasta.parent.mkdir(parents=True, exist_ok=True)
        if fasta.exists():
            fasta.unlink()

        with open(phy, "w", encoding="utf-8") as phy_fh:
            phy_fh.write("{} {}\n".format(len(keep), n_sites))
            for sample in keep:
                seq = seq_paths[sample].read_text(encoding="utf-8")
                phy_fh.write("{} {}\n".format(sample, seq))
                _write_fasta(fasta, sample, seq_paths[sample])

        with open(summary_path, "w", encoding="utf-8") as summary_fh:
            summary_fh.write("n_samples\tn_loci\tn_sites\tn_samples_in_keep_list\n")
            summary_fh.write("{}\t{}\t{}\t{}\n".format(len(keep), n_loci, n_sites, len(keep_set)))
        print(
            "Kept {} samples, {} loci, {} sites".format(len(keep), n_loci, n_sites),
            flush=True,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    log_handles = []
    try:
        log_path = snakemake.log[0] if getattr(snakemake, "log", None) else None
        if log_path:
            log_fh = open(log_path, "w", encoding="utf-8")
            log_handles.append(log_fh)
            sys.stdout = log_fh
            sys.stderr = log_fh
        print("Reading loci: {}".format(snakemake.input["loci"]), flush=True)
        print("Keep list: {}".format(snakemake.input["samples"]), flush=True)
        main(
            loci_path=snakemake.input["loci"],
            samples_path=snakemake.input["samples"],
            phy_path=snakemake.output["phy"],
            fasta_path=snakemake.output["fasta"],
            summary_path=snakemake.output["summary"],
        )
        print(
            "Wrote {} and {}".format(snakemake.output["phy"], snakemake.output["fasta"]),
            flush=True,
        )
    finally:
        for handle in log_handles:
            handle.close()
