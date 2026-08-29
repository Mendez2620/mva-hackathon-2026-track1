# MVA Hackathon 2026 — Track 1

Privacy-preserving, phenotype-first variant prioritization for the Rare
Disease, Real Kid MVA Hackathon 2026.

## Result

The current Track 1 hypothesis is a single, rank-1 compound-heterozygous
**BUB1B** pair. The submission-ready CSV is
`results/Mendez26_mva-local-evidence-v1.csv`; the scientific rationale and
limitations are in `docs/Mendez26_track1_report_v1.md`.

## Data protection

This repository intentionally contains **no VCF, FASTQ, BAM/CRAM, clinical
DOCX, patient-level variant table, prompt containing bulk variants, or model
trained on the patient's genome**. Never add those files. The supplied
`.gitignore` excludes common controlled-data formats and all `data/` and
`work/` directories.

The code downloads public reference data by complete gene interval or as a
complete predictor file. Patient variants are matched locally and are never
sent to an annotation API.

## Reproduce

Python 3.12 was used for the submitted analysis.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Keep the gated inputs outside the repository or under the ignored `data/`
directory. Then run aggregate QC:

```bash
.venv/bin/python scripts/scan_vcf_qc.py \
  /secure/path/WGS_EX2312012_HGWCNDSX7.vcf.gz \
  --output work/vcf_qc.aggregate.json
```

Run the MANE-centered MVA gene analysis. Variant-level outputs are written to
the ignored `work/private/` directory:

```bash
.venv/bin/python scripts/targeted_annotation.py \
  --vcf /secure/path/WGS_EX2312012_HGWCNDSX7.vcf.gz
```

Download complete public alternate-G predictor tracks. These downloads total
about 0.9 GB and contain no patient data:

```bash
bash scripts/download_public_predictors.sh
```

Join predictor scores locally and validate call-level evidence:

```bash
.venv/bin/python scripts/annotate_local_predictors.py \
  --input work/private/targeted_variants.private.json \
  --output work/private/targeted_variants.scored.private.json \
  --revel-g work/public_cache/predictors/revel_G.bw \
  --alphamissense-g work/public_cache/predictors/alphaMissense_G.bw

.venv/bin/python scripts/validate_candidate_calls.py \
  --vcf /secure/path/WGS_EX2312012_HGWCNDSX7.vcf.gz \
  --candidates work/private/targeted_variants.scored.private.json \
  --output work/private/bub1b_call_validation.private.json
```

Validate the shareable Track 1 CSV schema:

```bash
.venv/bin/python scripts/validate_submission.py \
  results/Mendez26_mva-local-evidence-v1.csv
```

## Repository contents

- `scripts/scan_vcf_qc.py`: streaming aggregate QC without coordinates.
- `scripts/targeted_annotation.py`: local MANE consequence prediction and
  local joins to public ClinVar/gnomAD gene-interval downloads.
- `scripts/annotate_local_predictors.py`: offline REVEL and AlphaMissense
  bigWig lookup.
- `scripts/validate_candidate_calls.py`: local VCF evidence review.
- `scripts/validate_submission.py`: Track 1 format checks without ground truth.
- `results/phenotype_hpo.tsv`: faithful HPO transcription of the supplied
  phenotype.
- `results/Mendez26_mva-local-evidence-v1.csv`: submission-ready prediction.
- `docs/Mendez26_track1_report_v1.md`: methods, evidence, and limitations.

## Research-use disclaimer

This is a hackathon research result, not a clinical diagnosis or treatment
recommendation. Clinical use requires accredited confirmation, segregation
analysis, and qualified genetics review.
