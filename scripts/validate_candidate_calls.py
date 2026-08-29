#!/usr/bin/env python3
"""Validate call-level evidence for locally identified private candidates.

No network access is used.  Coordinates are taken from the private candidate
JSON rather than from command-line arguments so they do not appear in process
metadata or external requests.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gene", default="BUB1B")
    return parser.parse_args()


def serializable(value: Any) -> Any:
    if isinstance(value, tuple):
        return [serializable(item) for item in value]
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def main() -> None:
    args = parse_args()
    candidate_rows = json.loads(args.candidates.read_text(encoding="utf-8"))
    candidates = [
        row
        for row in candidate_rows
        if row.get("gene") == args.gene and row.get("inheritance_flag")
    ]
    vcf = pysam.VariantFile(str(args.vcf))
    sample_name = list(vcf.header.samples)[0]
    results: list[dict[str, Any]] = []

    for candidate in candidates:
        chrom = str(candidate["chrom"])
        position = int(candidate["pos"])
        ref = str(candidate["ref"])
        alt = str(candidate["alt"])
        matches = [
            record
            for record in vcf.fetch(chrom, position - 1, position)
            if record.pos == position and record.ref == ref and alt in (record.alts or ())
        ]
        if len(matches) != 1:
            raise RuntimeError(f"Expected exactly one local VCF match, observed {len(matches)}")
        record = matches[0]
        alt_index = list(record.alts or ()).index(alt) + 1
        sample = record.samples[sample_name]
        ad = sample.get("AD")
        total_ad = sum(value or 0 for value in ad or ())
        filters = list(record.filter.keys())

        nearby = list(vcf.fetch(chrom, max(0, position - 101), position + 100))
        results.append(
            {
                "gene": args.gene,
                "chrom": chrom,
                "pos": position,
                "ref": ref,
                "alt": alt,
                "record_id": record.id or ".",
                "filter": "PASS" if not filters or filters == ["PASS"] else filters,
                "qual": serializable(record.qual),
                "info": {key: serializable(value) for key, value in record.info.items()},
                "sample": {
                    key: serializable(value) for key, value in sample.items()
                },
                "selected_alt_depth": serializable(ad[alt_index] if ad and alt_index < len(ad) else None),
                "total_allele_depth": total_ad,
                "selected_alt_balance": round((ad[alt_index] or 0) / total_ad, 6)
                if ad and alt_index < len(ad) and total_ad
                else None,
                "nearby_201bp": {
                    "record_count": len(nearby),
                    "pass_count": sum(
                        not list(item.filter.keys()) or list(item.filter.keys()) == ["PASS"]
                        for item in nearby
                    ),
                },
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
