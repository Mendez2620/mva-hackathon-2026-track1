#!/usr/bin/env python3
"""Validate the public Track 1 CSV contract without accessing ground truth."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    return parser.parse_args()


def validate_variant(row: dict[str, str], suffix: str, required: bool) -> None:
    values = [row[f"chrom_{suffix}"], row[f"pos_{suffix}"], row[f"ref_{suffix}"], row[f"alt_{suffix}"]]
    if not any(values) and not required:
        return
    if not all(values):
        raise ValueError(f"Variant {suffix} is only partially populated")
    if not row[f"chrom_{suffix}"].startswith("chr"):
        raise ValueError(f"chrom_{suffix} must use the chr prefix")
    if int(row[f"pos_{suffix}"]) < 1:
        raise ValueError(f"pos_{suffix} must be positive")
    for allele_field in (f"ref_{suffix}", f"alt_{suffix}"):
        if row[allele_field] != row[allele_field].upper():
            raise ValueError(f"{allele_field} must be uppercase")


def main() -> None:
    args = parse_args()
    with args.csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValueError(f"Unexpected header: {reader.fieldnames}")
        rows = list(reader)
    if not 1 <= len(rows) <= 10:
        raise ValueError("Submission must contain 1 to 10 rows")

    seen: set[tuple[str, ...]] = set()
    previous_epcr = float("inf")
    for index, row in enumerate(rows, start=1):
        if row["proband_id"] != "PROBAND01":
            raise ValueError(f"Row {index}: only PROBAND01 is accepted")
        validate_variant(row, "1", required=True)
        validate_variant(row, "2", required=False)
        epcr = float(row["epcr"])
        if not 0 < epcr <= 1:
            raise ValueError(f"Row {index}: EPCR must be in (0, 1]")
        if epcr > previous_epcr:
            raise ValueError("Rows must be ranked by non-increasing EPCR")
        previous_epcr = epcr
        if row["finding_type"] not in {"primary", "secondary"}:
            raise ValueError(f"Row {index}: invalid finding_type")
        key = tuple(row[field] for field in FIELDS[1:9])
        if key in seen:
            raise ValueError(f"Row {index}: duplicate prediction")
        seen.add(key)
    print(f"OK: {len(rows)} ranked row(s), schema valid")


if __name__ == "__main__":
    main()
