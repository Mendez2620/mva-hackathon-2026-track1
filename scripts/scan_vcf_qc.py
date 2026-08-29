#!/usr/bin/env python3
"""Stream a single-sample VCF and emit aggregate QC only.

No variant-level rows or coordinates are written. This keeps the intermediate
artifact suitable for sharing inside the analysis workspace.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path


PRIMARY_CONTIGS = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def add_hist(counter: Counter[str], value: int | None, bins: tuple[int, ...]) -> None:
    if value is None:
        counter["missing"] += 1
        return
    lower = 0
    for upper in bins:
        if value < upper:
            counter[f"{lower}-{upper - 1}"] += 1
            return
        lower = upper
    counter[f">={bins[-1]}"] += 1


def classify_allele(ref: str, alt: str) -> str:
    if alt.startswith("<") or "[" in alt or "]" in alt or alt == "*":
        return "symbolic_or_breakend"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == len(alt):
        return "mnv"
    if len(ref) < len(alt):
        return "insertion"
    return "deletion"


def main() -> None:
    args = parse_args()
    records = 0
    alleles = 0
    multiallelic = 0
    primary_records = 0
    nonprimary_records = 0
    dbsnp_records = 0
    transitions = 0
    transversions = 0
    high_quality_alt = 0
    called_alt = 0
    called_hom_ref = 0
    missing_gt = 0
    phased_gt = 0
    het_ab_sum = 0.0
    het_ab_count = 0
    out_of_order = 0
    previous_contig = None
    previous_pos = -1

    filters: Counter[str] = Counter()
    allele_types: Counter[str] = Counter()
    genotype_types: Counter[str] = Counter()
    chrom_records: Counter[str] = Counter()
    dp_hist: Counter[str] = Counter()
    gq_hist: Counter[str] = Counter()
    qual_hist: Counter[str] = Counter()
    header_counts: Counter[str] = Counter()
    references: list[str] = []
    samples: list[str] = []
    format_ids: list[str] = []
    info_ids: list[str] = []

    transition_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    with gzip.open(args.vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                header_counts[line[2:].split("=", 1)[0]] += 1
                if line.startswith("##reference="):
                    references.append(line.rstrip().split("=", 1)[1])
                elif line.startswith("##FORMAT=<ID="):
                    format_ids.append(line.split("##FORMAT=<ID=", 1)[1].split(",", 1)[0])
                elif line.startswith("##INFO=<ID="):
                    info_ids.append(line.split("##INFO=<ID=", 1)[1].split(",", 1)[0])
                continue
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                samples = columns[9:]
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos_text, record_id, ref, alt_text, qual_text, filt, info, fmt, sample = fields[:10]
            pos = int(pos_text)
            records += 1
            chrom_records[chrom] += 1
            if chrom in PRIMARY_CONTIGS:
                primary_records += 1
            else:
                nonprimary_records += 1

            if chrom == previous_contig and pos < previous_pos:
                out_of_order += 1
            previous_contig = chrom
            previous_pos = pos

            if record_id != "." or ";DB" in f";{info}" or info == "DB" or info.startswith("DB;"):
                dbsnp_records += 1

            filter_values = [filt] if filt in {"PASS", "."} else filt.split(";")
            filters.update(filter_values)

            alts = alt_text.split(",")
            alleles += len(alts)
            if len(alts) > 1:
                multiallelic += 1
            for alt in alts:
                kind = classify_allele(ref, alt)
                allele_types[kind] += 1
                if kind == "snv":
                    if (ref.upper(), alt.upper()) in transition_pairs:
                        transitions += 1
                    else:
                        transversions += 1

            try:
                qual_value = int(float(qual_text)) if qual_text != "." else None
            except ValueError:
                qual_value = None
            add_hist(qual_hist, qual_value, (10, 20, 30, 50, 100, 200, 500, 1000))

            format_keys = fmt.split(":")
            sample_values = sample.split(":")
            values = dict(zip(format_keys, sample_values))
            gt = values.get("GT", ".")
            if "|" in gt:
                phased_gt += 1
            normalized_gt = gt.replace("|", "/")
            genotype_types[normalized_gt] += 1
            allelic_gt = [x for x in normalized_gt.split("/") if x not in {".", ""}]
            if not allelic_gt:
                missing_gt += 1
            elif all(x == "0" for x in allelic_gt):
                called_hom_ref += 1
            elif any(x != "0" for x in allelic_gt):
                called_alt += 1

            def parse_int(name: str) -> int | None:
                raw = values.get(name, ".")
                try:
                    return int(raw) if raw not in {".", ""} else None
                except ValueError:
                    return None

            dp = parse_int("DP")
            gq = parse_int("GQ")
            add_hist(dp_hist, dp, (5, 10, 15, 20, 30, 40, 60, 100, 200))
            add_hist(gq_hist, gq, (10, 20, 30, 50, 70, 90, 99))

            if normalized_gt in {"0/1", "1/0"}:
                ad = values.get("AD", ".")
                if ad not in {".", ""}:
                    try:
                        depths = [int(x) for x in ad.split(",")]
                        total = sum(depths)
                        if total > 0 and len(depths) > 1:
                            het_ab_sum += depths[1] / total
                            het_ab_count += 1
                    except ValueError:
                        pass

            if filt in {"PASS", "."} and allelic_gt and any(x != "0" for x in allelic_gt):
                if (dp is None or dp >= 10) and (gq is None or gq >= 20):
                    high_quality_alt += 1

    ti_tv = transitions / transversions if transversions else math.nan
    result = {
        "input": args.vcf.name,
        "records": records,
        "alleles": alleles,
        "multiallelic_records": multiallelic,
        "primary_contig_records": primary_records,
        "nonprimary_contig_records": nonprimary_records,
        "dbsnp_or_named_records": dbsnp_records,
        "called_alt_records": called_alt,
        "called_hom_ref_records": called_hom_ref,
        "missing_gt_records": missing_gt,
        "phased_gt_records": phased_gt,
        "high_quality_alt_records": high_quality_alt,
        "transition_count": transitions,
        "transversion_count": transversions,
        "ti_tv_ratio": round(ti_tv, 4) if math.isfinite(ti_tv) else None,
        "mean_het_alt_balance": round(het_ab_sum / het_ab_count, 4) if het_ab_count else None,
        "het_alt_balance_records": het_ab_count,
        "out_of_order_records": out_of_order,
        "samples": samples,
        "references": references,
        "format_ids": format_ids,
        "info_ids": info_ids,
        "header_counts": dict(header_counts),
        "filters": dict(filters.most_common()),
        "allele_types": dict(allele_types.most_common()),
        "genotype_types": dict(genotype_types.most_common()),
        "chrom_records": dict(chrom_records.most_common()),
        "dp_histogram": dict(dp_hist),
        "gq_histogram": dict(gq_hist),
        "qual_histogram": dict(qual_hist),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
