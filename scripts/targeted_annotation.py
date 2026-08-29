#!/usr/bin/env python3
"""Privacy-preserving targeted annotation for the MVA hackathon.

The script downloads only public reference data for complete gene intervals.
Patient variants never leave the local process. Variant-level outputs are
written beneath a directory explicitly named ``private`` and must never be
committed or redistributed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pysam


UCSC_API = "https://api.genome.ucsc.edu"

# The first seven entries represent the currently named MVA1-MVA7 genes.
# BUB1/BUB3 are exploratory spindle-checkpoint genes with supportive evidence.
GENES: dict[str, dict[str, Any]] = {
    "BUB1B": {"inheritance": "AR", "evidence_tier": 5.0, "upstream": 60_000, "downstream": 10_000},
    "CEP57": {"inheritance": "AR", "evidence_tier": 4.5, "upstream": 10_000, "downstream": 10_000},
    "TRIP13": {"inheritance": "AR", "evidence_tier": 4.5, "upstream": 10_000, "downstream": 10_000},
    "CENATAC": {"inheritance": "AR", "evidence_tier": 3.5, "upstream": 10_000, "downstream": 10_000},
    "SLF2": {"inheritance": "AR", "evidence_tier": 3.0, "upstream": 10_000, "downstream": 10_000},
    "SMC5": {"inheritance": "AR", "evidence_tier": 3.0, "upstream": 10_000, "downstream": 10_000},
    "MAD1L1": {"inheritance": "AR", "evidence_tier": 5.0, "upstream": 10_000, "downstream": 10_000},
    "BUB1": {"inheritance": "AD/supportive", "evidence_tier": 2.0, "upstream": 10_000, "downstream": 10_000},
    "BUB3": {"inheritance": "AD/supportive", "evidence_tier": 2.0, "upstream": 10_000, "downstream": 10_000},
}

GNOMAD_TRACKS = {
    "gnomad_genomes_v4_1_1": "gnomadGenomesVariantsV4_1",
    "gnomad_exomes_v4_1_1": "gnomadExomesVariantsV4_1",
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=Path("work/public_cache"))
    parser.add_argument("--private-dir", type=Path, default=Path("work/private"))
    parser.add_argument("--summary", type=Path, default=Path("work/targeted_summary.json"))
    parser.add_argument("--refresh-public-cache", action="store_true")
    return parser.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def fetch_json(url: str, cache_path: Path, refresh: bool = False) -> dict[str, Any]:
    if cache_path.exists() and not refresh:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "mva-hackathon-local-analysis/0.1"})
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read().decode("utf-8")
            data = json.loads(payload)
            cache_path.write_text(json.dumps(data, separators=(",", ":")) + "\n", encoding="utf-8")
            return data
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Unable to fetch public reference URL after retries: {url}") from last_error


def ucsc_search(symbol: str, cache_dir: Path, refresh: bool) -> dict[str, Any]:
    url = f"{UCSC_API}/search?genome=hg38;search={urllib.parse.quote(symbol)}"
    return fetch_json(url, cache_dir / "search" / f"{symbol}.json", refresh)


def ucsc_track(
    track: str,
    chrom: str,
    start: int,
    end: int,
    cache_path: Path,
    refresh: bool,
) -> dict[str, Any]:
    url = (
        f"{UCSC_API}/getData/track?genome=hg38;track={track};"
        f"chrom=chr{chrom};start={start};end={end};maxItemsOutput=1000000"
    )
    return fetch_json(url, cache_path, refresh)


def ucsc_sequence(
    chrom: str,
    start: int,
    end: int,
    cache_path: Path,
    refresh: bool,
) -> str:
    url = f"{UCSC_API}/getData/sequence?genome=hg38;chrom=chr{chrom};start={start};end={end}"
    data = fetch_json(url, cache_path, refresh)
    return str(data["dna"]).upper()


def parse_display_position(position: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"chr([^:]+):(\d+)-(\d+)", position)
    if not match:
        raise ValueError(f"Unexpected UCSC position: {position}")
    chrom, start_1, end_1 = match.groups()
    return chrom, int(start_1) - 1, int(end_1)


def find_mane_transcript(symbol: str, cache_dir: Path, refresh: bool) -> dict[str, Any]:
    search = ucsc_search(symbol, cache_dir, refresh)
    candidates: list[tuple[str, int, int]] = []
    for group in search.get("positionMatches", []):
        if group.get("name") != "mane":
            continue
        for match in group.get("matches", []):
            if str(match.get("posName", "")).startswith(f"{symbol} ("):
                candidates.append(parse_display_position(match["position"]))
    if not candidates:
        raise RuntimeError(f"No MANE interval found for {symbol}")

    for chrom, start, end in candidates:
        cache_path = cache_dir / "mane" / f"{symbol}.{chrom}.{start}.{end}.json"
        data = ucsc_track("mane", chrom, max(0, start - 100), end + 100, cache_path, refresh)
        rows = [row for row in data.get("mane", []) if row.get("geneName2") == symbol]
        rows.sort(key=lambda row: (row.get("maneStat") != "MANE Select", row.get("chromStart", 0)))
        if rows:
            return rows[0]
    raise RuntimeError(f"MANE search returned no exact transcript row for {symbol}")


def split_ints(value: str) -> list[int]:
    return [int(item) for item in value.rstrip(",").split(",") if item]


def mane_to_transcript(row: dict[str, Any]) -> dict[str, Any]:
    tx_start = int(row["chromStart"])
    offsets = split_ints(str(row["chromStarts"]))
    sizes = split_ints(str(row["blockSizes"]))
    exon_starts = [tx_start + value for value in offsets]
    exon_ends = [start + size for start, size in zip(exon_starts, sizes)]
    return {
        "gene": row["geneName2"],
        "chrom": str(row["chrom"]).removeprefix("chr"),
        "strand": row["strand"],
        "transcript": row.get("ncbiId") or row.get("name"),
        "ensembl_transcript": row.get("name"),
        "mane_status": row.get("maneStat"),
        "tx_start": tx_start,
        "tx_end": int(row["chromEnd"]),
        "cds_start": int(row["thickStart"]),
        "cds_end": int(row["thickEnd"]),
        "exon_starts": exon_starts,
        "exon_ends": exon_ends,
    }


def scan_interval(transcript: dict[str, Any], gene_config: dict[str, Any]) -> tuple[int, int]:
    if transcript["strand"] == "+":
        return (
            max(0, transcript["tx_start"] - int(gene_config["upstream"])),
            transcript["tx_end"] + int(gene_config["downstream"]),
        )
    return (
        max(0, transcript["tx_start"] - int(gene_config["downstream"])),
        transcript["tx_end"] + int(gene_config["upstream"]),
    )


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(sequence: str) -> str:
    return "".join(CODON_TABLE.get(sequence[i : i + 3], "X") for i in range(0, len(sequence) - 2, 3))


def overlaps(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


@dataclass
class TranscriptModel:
    metadata: dict[str, Any]
    genomic_sequence: str

    def __post_init__(self) -> None:
        self.chrom = self.metadata["chrom"]
        self.strand = self.metadata["strand"]
        self.tx_start = int(self.metadata["tx_start"])
        self.tx_end = int(self.metadata["tx_end"])
        self.cds_start = int(self.metadata["cds_start"])
        self.cds_end = int(self.metadata["cds_end"])
        self.exons = list(zip(self.metadata["exon_starts"], self.metadata["exon_ends"]))

        coding_segments = [
            (max(start, self.cds_start), min(end, self.cds_end))
            for start, end in self.exons
            if max(start, self.cds_start) < min(end, self.cds_end)
        ]
        positions: list[int] = []
        if self.strand == "+":
            for start, end in coding_segments:
                positions.extend(range(start, end))
        else:
            for start, end in reversed(coding_segments):
                positions.extend(range(end - 1, start - 1, -1))
        self.cds_positions = positions
        self.cds_index = {position: index for index, position in enumerate(positions)}
        bases = [self.genomic_sequence[position - self.tx_start] for position in positions]
        self.cds_sequence = "".join(bases if self.strand == "+" else [revcomp(base) for base in bases])

    def _nearest_splice_boundary(self, var_start: int, var_end: int) -> tuple[int, str] | None:
        options: list[tuple[int, str]] = []
        for index, (exon_start, exon_end) in enumerate(self.exons):
            if index > 0:
                distance = exon_start - var_end + 1 if var_end <= exon_start else var_start - exon_start + 1
                label = "acceptor" if self.strand == "+" else "donor"
                options.append((abs(distance), label))
            if index < len(self.exons) - 1:
                distance = var_start - exon_end + 1 if var_start >= exon_end else exon_end - var_end + 1
                label = "donor" if self.strand == "+" else "acceptor"
                options.append((abs(distance), label))
        return min(options, default=None)

    def annotate(self, pos_1: int, ref: str, alt: str) -> dict[str, Any]:
        var_start = pos_1 - 1
        var_end = var_start + len(ref)
        result: dict[str, Any] = {
            "region_consequence": "intergenic",
            "coding_consequence": None,
            "hgvsc_approx": None,
            "hgvsp_approx": None,
            "reference_match": None,
        }
        if var_end <= self.tx_start:
            result["region_consequence"] = "upstream" if self.strand == "+" else "downstream"
            return result
        if var_start >= self.tx_end:
            result["region_consequence"] = "downstream" if self.strand == "+" else "upstream"
            return result

        exon_overlap = any(overlaps(var_start, var_end, start, end) for start, end in self.exons)
        genomic_positions = list(range(var_start, var_end))
        all_coding = bool(genomic_positions) and all(position in self.cds_index for position in genomic_positions)

        if all_coding and not alt.startswith("<") and "[" not in alt and "]" not in alt:
            return self._annotate_coding(pos_1, ref, alt, genomic_positions)
        if exon_overlap:
            result["region_consequence"] = "exonic_non_coding_or_utr"
            return result

        nearest = self._nearest_splice_boundary(var_start, var_end)
        if nearest and nearest[0] <= 8:
            distance, boundary = nearest
            result["region_consequence"] = (
                f"splice_{boundary}_site" if distance <= 2 else f"splice_{boundary}_region"
            )
            return result
        result["region_consequence"] = "intronic"
        return result

    def _annotate_coding(
        self,
        pos_1: int,
        ref: str,
        alt: str,
        genomic_positions: list[int],
    ) -> dict[str, Any]:
        indices = sorted(self.cds_index[position] for position in genomic_positions)
        first = indices[0]
        ref_tx = ref.upper() if self.strand == "+" else revcomp(ref.upper())
        alt_tx = alt.upper() if self.strand == "+" else revcomp(alt.upper())
        reference_match = self.cds_sequence[first : first + len(ref_tx)] == ref_tx
        result: dict[str, Any] = {
            "region_consequence": "coding",
            "coding_consequence": "coding_reference_mismatch" if not reference_match else None,
            "hgvsc_approx": None,
            "hgvsp_approx": None,
            "reference_match": reference_match,
        }
        c_start = first + 1
        c_end = first + len(ref_tx)
        transcript = self.metadata["transcript"]
        if len(ref_tx) == len(alt_tx) == 1:
            result["hgvsc_approx"] = f"{transcript}:c.{c_start}{ref_tx}>{alt_tx}"
        else:
            span = str(c_start) if c_start == c_end else f"{c_start}_{c_end}"
            result["hgvsc_approx"] = f"{transcript}:c.{span}delins{alt_tx}"
        if not reference_match:
            return result

        mutated = self.cds_sequence[:first] + alt_tx + self.cds_sequence[first + len(ref_tx) :]
        reference_protein = translate(self.cds_sequence)
        alternate_protein = translate(mutated)
        first_aa = min(first // 3, max(0, len(reference_protein) - 1))
        while (
            first_aa < len(reference_protein)
            and first_aa < len(alternate_protein)
            and reference_protein[first_aa] == alternate_protein[first_aa]
        ):
            first_aa += 1
        aa_position = first_aa + 1
        ref_aa = reference_protein[first_aa] if first_aa < len(reference_protein) else "X"
        alt_aa = alternate_protein[first_aa] if first_aa < len(alternate_protein) else "*"

        delta = len(alt_tx) - len(ref_tx)
        if delta % 3 != 0:
            result["coding_consequence"] = "frameshift_variant"
            stop_index = alternate_protein.find("*", first_aa)
            suffix = f"fsTer{stop_index - first_aa + 1}" if stop_index >= 0 else "fs"
            result["hgvsp_approx"] = (
                f"{self.metadata['transcript']}:p.{AA3.get(ref_aa, 'Xaa')}{aa_position}{suffix}"
            )
        elif len(ref_tx) != len(alt_tx):
            result["coding_consequence"] = "inframe_indel"
            result["hgvsp_approx"] = (
                f"{self.metadata['transcript']}:p.{AA3.get(ref_aa, 'Xaa')}{aa_position}delins"
                f"{AA3.get(alt_aa, 'Xaa')}"
            )
        elif ref_aa == alt_aa:
            result["coding_consequence"] = "synonymous_variant"
            result["hgvsp_approx"] = f"{self.metadata['transcript']}:p.{AA3.get(ref_aa, 'Xaa')}{aa_position}="
        elif alt_aa == "*":
            result["coding_consequence"] = "stop_gained"
            result["hgvsp_approx"] = (
                f"{self.metadata['transcript']}:p.{AA3.get(ref_aa, 'Xaa')}{aa_position}Ter"
            )
        elif ref_aa == "*":
            result["coding_consequence"] = "stop_lost"
            result["hgvsp_approx"] = (
                f"{self.metadata['transcript']}:p.Ter{aa_position}{AA3.get(alt_aa, 'Xaa')}ext"
            )
        elif aa_position == 1 and ref_aa == "M":
            result["coding_consequence"] = "start_lost"
            result["hgvsp_approx"] = (
                f"{self.metadata['transcript']}:p.Met1{AA3.get(alt_aa, 'Xaa')}"
            )
        else:
            result["coding_consequence"] = "missense_variant"
            result["hgvsp_approx"] = (
                f"{self.metadata['transcript']}:p.{AA3.get(ref_aa, 'Xaa')}{aa_position}"
                f"{AA3.get(alt_aa, 'Xaa')}"
            )
        return result


def variant_key(chrom: str, pos_1: int, ref: str, alt: str) -> tuple[str, int, str, str]:
    return chrom.removeprefix("chr"), int(pos_1), ref.upper(), alt.upper()


def gnomad_index(rows: Iterable[dict[str, Any]]) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    index: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    for row in rows:
        if not all(field in row for field in ("chrom", "ref", "alt")):
            continue
        pos_1 = row.get("_startPos", int(row.get("chromStart", 0)) + 1)
        index[variant_key(str(row["chrom"]), int(pos_1), str(row["ref"]), str(row["alt"]))] = row
    return index


def clinvar_index(rows: Iterable[dict[str, Any]]) -> tuple[
    dict[tuple[str, int, str, str], list[dict[str, Any]]], dict[str, list[dict[str, Any]]]
]:
    by_variant: dict[tuple[str, int, str, str], list[dict[str, Any]]] = {}
    by_rsid: dict[str, list[dict[str, Any]]] = {}
    pattern = re.compile(r"chr([^:]+):(\d+):([^>]+)>(.+)")
    for row in rows:
        description = str(row.get("vcfDesc", ""))
        match = pattern.fullmatch(description)
        if match:
            chrom, pos_1, ref, alt = match.groups()
            by_variant.setdefault(variant_key(chrom, int(pos_1), ref, alt), []).append(row)
        rsid = str(row.get("snpId", ""))
        if rsid:
            by_rsid.setdefault(rsid, []).append(row)
    return by_variant, by_rsid


def parse_float(value: Any) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def max_population_af(row: dict[str, Any]) -> float | None:
    values = [parse_float(row.get("AF")), parse_float(row.get("AF_grpmax"))]
    present = [value for value in values if value is not None]
    return max(present) if present else None


def severity(consequence: str | None, region: str) -> float:
    if consequence in {"stop_gained", "frameshift_variant", "start_lost", "stop_lost"}:
        return 5.0
    if region in {"splice_donor_site", "splice_acceptor_site"}:
        return 5.0
    if consequence in {"missense_variant", "inframe_indel"}:
        return 3.0
    if region in {"splice_donor_region", "splice_acceptor_region"}:
        return 2.5
    if consequence == "synonymous_variant":
        return 0.5
    return 0.0


def ranking_score(row: dict[str, Any], config: dict[str, Any]) -> float:
    score = float(config["evidence_tier"])
    score += severity(row.get("coding_consequence"), row["region_consequence"])
    gt = row.get("genotype", "")
    if config["inheritance"] == "AR":
        score += 3.0 if gt in {"1/1", "2/2"} else 1.0
    else:
        score += 1.5 if gt not in {"0/0", "./."} else 0.0

    af = row.get("max_gnomad_af")
    if af is None:
        score += 3.0
    elif af < 1e-4:
        score += 2.5
    elif af < 1e-3:
        score += 2.0
    elif af < 1e-2:
        score += 0.5
    else:
        score -= 5.0

    clinvar_text = " | ".join(row.get("clinvar_significance", [])).lower()
    if "pathogenic" in clinvar_text and "conflicting" not in clinvar_text:
        score += 5.0
    if "benign" in clinvar_text and "conflicting" not in clinvar_text:
        score -= 5.0
    if row["filter"] == "PASS":
        score += 1.0
    else:
        score -= 3.0
    if row["gene"] == "BUB1B" and row["region_consequence"] == "upstream":
        score += 1.0
    return round(score, 3)


def fmt_gt(gt: tuple[int | None, ...] | None) -> str:
    if not gt:
        return "./."
    return "/".join("." if allele is None else str(allele) for allele in gt)


def extract_patient_rows(
    vcf: pysam.VariantFile,
    transcript: dict[str, Any],
    model: TranscriptModel,
    interval: tuple[int, int],
) -> list[dict[str, Any]]:
    gene = transcript["gene"]
    chrom = transcript["chrom"]
    sample_name = list(vcf.header.samples)[0]
    rows: list[dict[str, Any]] = []
    for record in vcf.fetch(chrom, interval[0], interval[1]):
        sample = record.samples[sample_name]
        gt = sample.get("GT")
        present_alleles = {allele for allele in (gt or ()) if allele not in (None, 0)}
        for alt_index, alt in enumerate(record.alts or (), start=1):
            if alt_index not in present_alleles:
                continue
            ad = sample.get("AD")
            dp = sample.get("DP")
            gq = sample.get("GQ")
            alt_balance = None
            if ad and alt_index < len(ad) and sum(value or 0 for value in ad) > 0:
                alt_balance = (ad[alt_index] or 0) / sum(value or 0 for value in ad)
            filters = list(record.filter.keys())
            filter_text = "PASS" if not filters or filters == ["PASS"] else ";".join(filters)
            annotation = model.annotate(record.pos, record.ref, alt)
            rows.append(
                {
                    "gene": gene,
                    "chrom": chrom,
                    "pos": int(record.pos),
                    "id": record.id or ".",
                    "ref": record.ref,
                    "alt": alt,
                    "qual": float(record.qual) if record.qual is not None else None,
                    "filter": filter_text,
                    "genotype": fmt_gt(gt),
                    "dp": int(dp) if dp is not None else None,
                    "gq": int(gq) if gq is not None else None,
                    "allele_balance": round(alt_balance, 4) if alt_balance is not None else None,
                    "transcript": transcript["transcript"],
                    **annotation,
                }
            )
    return rows


def attach_public_annotations(
    rows: list[dict[str, Any]],
    public_data: dict[str, dict[str, Any]],
) -> None:
    gnomad_indexes = {
        label: gnomad_index(data.get(track, []))
        for label, (track, data) in public_data["gnomad"].items()
    }
    clinvar_track, clinvar_data = public_data["clinvar"]
    clin_by_variant, clin_by_rsid = clinvar_index(clinvar_data.get(clinvar_track, []))

    for row in rows:
        key = variant_key(row["chrom"], row["pos"], row["ref"], row["alt"])
        af_values: list[float] = []
        row["gnomad"] = {}
        for label, index in gnomad_indexes.items():
            match = index.get(key)
            if not match:
                row["gnomad"][label] = None
                continue
            af = max_population_af(match)
            if af is not None:
                af_values.append(af)
            row["gnomad"][label] = {
                "AF": parse_float(match.get("AF")),
                "AF_grpmax": parse_float(match.get("AF_grpmax")),
                "nhomalt": match.get("nhomalt"),
                "rsId": match.get("rsId"),
                "variation_type": match.get("variation_type"),
            }
        row["max_gnomad_af"] = max(af_values) if af_values else None

        matches = list(clin_by_variant.get(key, []))
        if not matches and str(row.get("id", "")).startswith("rs"):
            matches = list(clin_by_rsid.get(str(row["id"]), []))
        row["clinvar_significance"] = sorted({str(item.get("clinSign")) for item in matches if item.get("clinSign")})
        row["clinvar_vcv"] = sorted({str(item.get("vcvId")) for item in matches if item.get("vcvId")})
        row["clinvar_hgvs_protein"] = sorted({str(item.get("_hgvsProt")) for item in matches if item.get("_hgvsProt")})


def is_rare_damaging(row: dict[str, Any]) -> bool:
    af = row.get("max_gnomad_af")
    rare = af is None or af < 0.01
    damaging = severity(row.get("coding_consequence"), row["region_consequence"]) >= 2.5
    clinvar_pathogenic = any("pathogenic" in value.lower() for value in row.get("clinvar_significance", []))
    return row["filter"] == "PASS" and rare and (damaging or clinvar_pathogenic)


def mark_compound_heterozygotes(rows: list[dict[str, Any]]) -> None:
    by_gene: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_gene.setdefault(row["gene"], []).append(row)
    for gene, gene_rows in by_gene.items():
        if GENES[gene]["inheritance"] != "AR":
            continue
        heterozygous = [row for row in gene_rows if row["genotype"] in {"0/1", "0/2", "1/2"} and is_rare_damaging(row)]
        if len(heterozygous) >= 2:
            for row in heterozygous:
                row["inheritance_flag"] = "possible_compound_heterozygous_unphased"
        for row in gene_rows:
            if "inheritance_flag" not in row:
                row["inheritance_flag"] = "homozygous_recessive_candidate" if row["genotype"] == "1/1" and is_rare_damaging(row) else None


def write_private_outputs(rows: list[dict[str, Any]], private_dir: Path) -> tuple[Path, Path]:
    private_dir.mkdir(parents=True, exist_ok=True)
    json_path = private_dir / "targeted_variants.private.json"
    tsv_path = private_dir / "targeted_variants.private.tsv"
    json_path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    fields = [
        "gene", "chrom", "pos", "id", "ref", "alt", "filter", "genotype", "dp", "gq",
        "allele_balance", "transcript", "region_consequence", "coding_consequence",
        "hgvsc_approx", "hgvsp_approx", "max_gnomad_af", "clinvar_significance",
        "clinvar_vcv", "inheritance_flag", "ranking_score",
    ]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# PRIVATE CONTROLLED DATA - DO NOT COMMIT OR REDISTRIBUTE\n")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            flat = dict(row)
            for field in ("clinvar_significance", "clinvar_vcv"):
                flat[field] = "|".join(flat.get(field, []))
            writer.writerow(flat)
    return json_path, tsv_path


def main() -> None:
    args = parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    transcripts: dict[str, dict[str, Any]] = {}
    intervals: dict[str, tuple[int, int]] = {}
    models: dict[str, TranscriptModel] = {}
    provenance: dict[str, Any] = {"ucsc_api": UCSC_API, "genes": {}}

    for gene, config in GENES.items():
        mane = find_mane_transcript(gene, args.cache_dir, args.refresh_public_cache)
        transcript = mane_to_transcript(mane)
        interval = scan_interval(transcript, config)
        sequence = ucsc_sequence(
            transcript["chrom"],
            transcript["tx_start"],
            transcript["tx_end"],
            args.cache_dir / "sequence" / f"{gene}.{transcript['transcript']}.json",
            args.refresh_public_cache,
        )
        transcripts[gene] = transcript
        intervals[gene] = interval
        models[gene] = TranscriptModel(transcript, sequence)
        provenance["genes"][gene] = {
            "transcript": transcript,
            "scan_interval_0_based_half_open": [transcript["chrom"], interval[0], interval[1]],
        }

    vcf = pysam.VariantFile(str(args.vcf))
    all_rows: list[dict[str, Any]] = []
    for gene in GENES:
        transcript = transcripts[gene]
        interval = intervals[gene]
        gene_rows = extract_patient_rows(vcf, transcript, models[gene], interval)

        gnomad_data: dict[str, tuple[str, dict[str, Any]]] = {}
        for label, track in GNOMAD_TRACKS.items():
            cache_path = args.cache_dir / "tracks" / track / f"{gene}.{interval[0]}.{interval[1]}.json"
            data = ucsc_track(
                track,
                transcript["chrom"],
                interval[0],
                interval[1],
                cache_path,
                args.refresh_public_cache,
            )
            gnomad_data[label] = (track, data)

        clinvar_track = "clinvarMain"
        clinvar_cache = args.cache_dir / "tracks" / clinvar_track / f"{gene}.{interval[0]}.{interval[1]}.json"
        clinvar_data = ucsc_track(
            clinvar_track,
            transcript["chrom"],
            interval[0],
            interval[1],
            clinvar_cache,
            args.refresh_public_cache,
        )
        attach_public_annotations(
            gene_rows,
            {"gnomad": gnomad_data, "clinvar": (clinvar_track, clinvar_data)},
        )
        for row in gene_rows:
            row["ranking_score"] = ranking_score(row, GENES[gene])
        all_rows.extend(gene_rows)

    mark_compound_heterozygotes(all_rows)
    all_rows.sort(key=lambda row: (-float(row["ranking_score"]), row["gene"], row["pos"]))
    json_path, tsv_path = write_private_outputs(all_rows, args.private_dir)

    summary = {
        "vcf": args.vcf.name,
        "sample": list(vcf.header.samples)[0],
        "genes_analyzed": list(GENES),
        "transcripts": {gene: transcripts[gene]["transcript"] for gene in GENES},
        "variant_rows": len(all_rows),
        "rows_by_gene": {gene: sum(row["gene"] == gene for row in all_rows) for gene in GENES},
        "rare_damaging_rows_by_gene": {
            gene: sum(row["gene"] == gene and is_rare_damaging(row) for row in all_rows) for gene in GENES
        },
        "inheritance_flags": {
            gene: sorted({row["inheritance_flag"] for row in all_rows if row["gene"] == gene and row.get("inheritance_flag")})
            for gene in GENES
        },
        "public_cache_sha256": hashlib.sha256(
            json.dumps(provenance, sort_keys=True).encode("utf-8")
        ).hexdigest(),
        "private_outputs": [json_path.name, tsv_path.name],
        "privacy_warning": "Variant-level files are controlled data and must not be committed or redistributed.",
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.cache_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
