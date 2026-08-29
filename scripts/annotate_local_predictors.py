#!/usr/bin/env python3
"""Attach locally downloaded whole-dataset pathogenicity scores.

This utility never accesses the network.  It reads patient variants from the
private targeted-analysis JSON, looks up scores in complete public bigWig
files already present on disk, and writes another private JSON file.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import pyBigWig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--revel-g", type=Path)
    parser.add_argument("--alphamissense-g", type=Path)
    return parser.parse_args()


def score_at(handle: pyBigWig.bigWigFile, chrom: str, pos_1: int) -> float | None:
    values = handle.values(f"chr{chrom}", pos_1 - 1, pos_1)
    if not values or values[0] is None or not math.isfinite(values[0]):
        return None
    return round(float(values[0]), 6)


def main() -> None:
    args = parse_args()
    rows: list[dict[str, Any]] = json.loads(args.input.read_text(encoding="utf-8"))
    resources: list[tuple[str, str, Path]] = []
    if args.revel_g:
        resources.append(("revel", "G", args.revel_g))
    if args.alphamissense_g:
        resources.append(("alphamissense", "G", args.alphamissense_g))

    handles: list[tuple[str, str, pyBigWig.bigWigFile]] = []
    try:
        for label, alternate, path in resources:
            if not path.is_file():
                raise FileNotFoundError(path)
            handles.append((label, alternate, pyBigWig.open(str(path))))
        for row in rows:
            row.setdefault("local_predictors", {})
            if len(str(row["ref"])) != 1 or len(str(row["alt"])) != 1:
                continue
            for label, alternate, handle in handles:
                if row["alt"].upper() == alternate:
                    row["local_predictors"][label] = score_at(
                        handle, str(row["chrom"]), int(row["pos"])
                    )
    finally:
        for _, _, handle in handles:
            handle.close()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
