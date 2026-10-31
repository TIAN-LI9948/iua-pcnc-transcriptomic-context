#!/usr/bin/env python3
"""Validate reproduced primary endpoints against the frozen reference outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from pathlib import Path


REQUIRED_NUMERIC_COLUMNS = {
    "difference_PC_minus_NC": "difference_PC_minus_NC",
    "bootstrap_ci95_low": "bootstrap_ci95_low",
    "bootstrap_ci95_high": "bootstrap_ci95_high",
    "exact_permutation_p_two_sided": "exact_permutation_p_two_sided",
    "BH_q_within_family": "exact_BH_q_within_family",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> int:
    repo = Path(__file__).resolve().parent.parent
    reproduced = repo / "reproduced_outputs" / "primary_endpoint_tests.tsv"
    expected_primary = repo / "expected_outputs" / "primary_endpoint_tests.tsv"
    expected_secondary = repo / "expected_outputs" / "secondary_endpoint_tests.tsv"

    actual_rows = read_tsv(reproduced)
    expected_rows = read_tsv(expected_primary) + read_tsv(expected_secondary)
    actual = {row["endpoint"]: row for row in actual_rows}
    expected = {row["endpoint"]: row for row in expected_rows}
    endpoint_order = [
        "PGR_logCPM",
        "frozen_candidate25_score",
        "ESR1_logCPM",
        "stable_epithelial_core_score",
    ]

    checks: list[dict[str, object]] = []
    all_pass = True
    for endpoint in endpoint_order:
        for actual_column, expected_column in REQUIRED_NUMERIC_COLUMNS.items():
            observed = float(actual[endpoint][actual_column])
            target = float(expected[endpoint][expected_column])
            difference = abs(observed - target)
            passed = difference <= 1e-12
            all_pass &= passed
            checks.append(
                {
                    "endpoint": endpoint,
                    "field": actual_column,
                    "reproduced": observed,
                    "expected": target,
                    "absolute_difference": difference,
                    "status": "PASS" if passed else "FAIL",
                }
            )

    manifest_rows: list[dict[str, object]] = []
    for folder in ("data", "scripts", "expected_outputs", "reproduced_outputs"):
        for path in sorted((repo / folder).rglob("*")):
            if path.is_file():
                manifest_rows.append(
                    {
                        "relative_path": path.relative_to(repo).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )

    metadata_dir = repo / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = metadata_dir / "file_manifest_sha256.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    report = {
        "status": "PASS" if all_pass else "FAIL",
        "tolerance": 1e-12,
        "endpoint_count": len(endpoint_order),
        "numeric_check_count": len(checks),
        "checks": checks,
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
        },
        "manifest": manifest_path.relative_to(repo).as_posix(),
    }
    report_path = metadata_dir / "reproduction_validation.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": len(checks), "manifest_files": len(manifest_rows)}))
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
