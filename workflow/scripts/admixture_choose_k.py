"""Parse ADMIXTURE cross-validation errors from per-K log files."""

import re
import sys
from pathlib import Path

CV_PATTERNS = (
    re.compile(r"Cross-validation error:\s*([0-9.eE+-]+)"),
    re.compile(r"CV error[^:]*:\s*([0-9.eE+-]+)"),
)


def parse_cv_error(log_path: Path) -> float:
    text = log_path.read_text()
    for pattern in CV_PATTERNS:
        matches = pattern.findall(text)
        if matches:
            return float(matches[-1])
    raise ValueError(f"No ADMIXTURE CV error found in {log_path}")


def parse_k_from_log(log_path: Path) -> int:
    match = re.search(r"\.(\d+)\.log$", log_path.name)
    if not match:
        raise ValueError(f"Could not parse K from log filename: {log_path}")
    return int(match.group(1))


def main(log_files, choose_k_file, cv_summary_file) -> None:
    rows: list[tuple[int, float]] = []
    errors: list[str] = []

    for log_path in sorted(log_files, key=parse_k_from_log):
        k = parse_k_from_log(log_path)
        try:
            cv_error = parse_cv_error(log_path)
        except ValueError as exc:
            errors.append(str(exc))
            continue
        rows.append((k, cv_error))

    if not rows:
        message = "No ADMIXTURE CV errors could be parsed from log files."
        if errors:
            message += " " + " ".join(errors)
        raise RuntimeError(message)

    rows.sort(key=lambda item: item[0])
    best_k = min(rows, key=lambda item: item[1])[0]

    cv_summary_file.parent.mkdir(parents=True, exist_ok=True)
    with cv_summary_file.open("w", encoding="utf-8") as handle:
        handle.write("K\tmedian\tmin\tmax\n")
        for k, value in rows:
            handle.write(f"{k}\t{value}\t{value}\t{value}\n")

    choose_k_lines = [
        "Choose the lowest cross-validation (CV) error across K values.",
        "Smaller values indicate better predictive accuracy.",
        f"Suggested K (minimum CV error): {best_k}",
        "",
        "K\tCV_error",
    ]
    choose_k_lines.extend(f"{k}\t{value:g}" for k, value in rows)
    choose_k_file.write_text("\n".join(choose_k_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    try:
        log_files = [Path(p) for p in snakemake.input]
        choose_k_file = Path(snakemake.output["choose_k_results"])
        cv_summary_file = Path(snakemake.output["cv_summary"])
    except NameError:
        print("This script is intended to run under Snakemake.", file=sys.stderr)
        sys.exit(1)

    main(log_files, choose_k_file, cv_summary_file)
