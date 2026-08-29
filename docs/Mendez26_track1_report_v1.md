# Track 1 report — privacy-preserving phenotype-first prioritization

**Participant:** Mendez26  
**Challenge:** Rare Disease, Real Kid: The MVA Hackathon 2026  
**Proband:** PROBAND01  
**Assembly:** GRCh38  
**Analysis version:** 1.0 (2026-08-28)

## Executive finding

The leading hypothesis is a **compound-heterozygous BUB1B finding** consisting of one heterozygous nonsense variant and one heterozygous ultrarare missense variant in the BUBR1 kinase domain. The pair is ranked first with EPCR 0.97.

| Rank | Gene | GRCh38 variant | MANE consequence | Genotype | Main evidence |
|---:|---|---|---|---|---|
| 1a | BUB1B | chr15:40209701 T>G | NM_001211.6:c.2210T>G, p.Leu737Ter | heterozygous | High-quality call; ClinVar Pathogenic/Likely pathogenic with multiple submitters and no conflicts; predicted loss of function |
| 1b | BUB1B | chr15:40220612 T>G | NM_001211.6:c.3006T>G, p.Asn1002Lys | heterozygous | High-quality call; ultrarare with no observed homozygotes in the queried gnomAD v4.1.1 tracks; protein-kinase domain; AlphaMissense 0.9229; REVEL 0.472 |

The variants are unphased in the supplied VCF. Therefore, *in trans* configuration is a strong but unconfirmed inference and should be tested by parental segregation or an appropriate read-backed/long-read phasing method.

## Phenotype model

The controlled clinical document was transcribed into the following HPO terms without adding unreported features:

- HP:0002859 — rhabdomyosarcoma
- HP:0000121 — nephrocalcinosis
- HP:0004322 — short stature
- HP:0001508 — failure to thrive
- HP:0003202 — skeletal muscle atrophy
- HP:0001622 — premature birth
- HP:0001518 — small for gestational age / severe IUGR
- HP:0200067 — recurrent spontaneous abortion in the parental history

This constellation strongly prioritizes a chromosomal-instability disorder. Biallelic BUB1B disease is an established cause of mosaic variegated aneuploidy with growth restriction, congenital abnormalities, and childhood cancer. The original BUB1B report included individuals with embryonal rhabdomyosarcoma. A later allelic study notes that affected individuals commonly carry a nonsense allele together with a kinase-domain missense allele, which closely matches the observed pair.

## Methods

### 1. Privacy and data handling

All patient-level VCF parsing, genotype inspection, consequence prediction, and public-data matching were performed locally. No patient variant, coordinate list, VCF subset, or genomic sequence was submitted to an external annotation service. Public reference data were downloaded by complete gene interval or as complete predictor tracks, then matched locally. Private variant-level intermediates were isolated under a directory excluded from publication.

### 2. Input integrity and genome-wide call-set QC

The supplied single-sample VCF is VCFv4.2, aligned to the GRCh38 no-alt analysis set (no `chr` prefix), and indexed successfully. It contains 5,012,204 records and 5,079,064 alternate alleles. Of the records, 4,740,790 are PASS. Aggregate properties were internally consistent for a whole-genome germline call set: 4,087,464 SNVs, Ti/Tv 1.9413, mean heterozygous alternate balance 0.4782, and 4,409,613 genotypes with GQ at least 99. The file streamed to completion without a decompression or ordering error.

### 3. Targeted disease-gene screen

The first pass covered the seven currently named MVA genes — **BUB1B, CEP57, TRIP13, CENATAC, SLF2, SMC5, and MAD1L1** — plus exploratory spindle-checkpoint genes BUB1 and BUB3. For each gene, the MANE Select transcript and reference sequence were retrieved from UCSC. The local annotator evaluated coding, splice-region, UTR, intronic, and regulatory-window variants. BUB1B was scanned with an extended upstream window because pathogenic regulatory alleles have been reported in this disease.

Public ClinVar and gnomAD v4.1.1 exome/genome annotations were downloaded for complete gene intervals and joined locally on normalized GRCh38 REF/ALT keys. Whole public REVEL and AlphaMissense alternate-G bigWig tracks were downloaded before local score lookup.

### 4. Variant-level call validation

| Metric | p.Leu737Ter | p.Asn1002Lys |
|---|---:|---:|
| FILTER | PASS | PASS |
| QUAL | 708.77 | 344.77 |
| DP | 46 | 28 |
| AD (REF, ALT) | 21, 25 | 15, 13 |
| Alternate balance | 0.5435 | 0.4643 |
| GQ | 99 | 99 |
| MQ | 60 | 60 |
| QD | 15.41 | 12.31 |
| FS | 1.115 | 1.617 |
| SOR | 0.523 | 0.330 |

Both calls are isolated within their local 201-bp windows and show balanced heterozygous support with strong genotype-likelihood separation. These observations reduce, but do not eliminate, the risk of a short-read calling artifact.

## Evidence interpretation

### BUB1B p.Leu737Ter

ClinVar identifies the exact MANE variant as **Pathogenic/Likely pathogenic**, based on criteria from multiple submitters without conflicts, for BUB1B-related disease including mosaic variegated aneuploidy syndrome 1. It introduces a premature stop before the annotated protein-kinase domain. ClinVar also cautions that population frequency estimates at this site may be unreliable because of gnomAD quality metrics; consequently, pathogenicity is based on loss-of-function mechanism and clinical assertions rather than treating the population frequency as decisive.

### BUB1B p.Asn1002Lys

This variant was not present in the downloaded ClinVar interval. In the queried gnomAD v4.1.1 tracks it was observed only at approximately 6.84e-7 overall exome AF (maximum population AF approximately 8.99e-7), with no homozygotes. Residue 1002 lies inside the UniProt-annotated BUBR1 protein-kinase domain (residues 766–1050). AlphaMissense scored the substitution 0.9229, above its genome-wide likely-pathogenic threshold of 0.564; REVEL was intermediate at 0.472. These computational results are supportive rather than diagnostic.

An ACMG-style evidence inventory supports rarity, location in a disease-relevant functional domain, computational impact, and a highly concordant phenotype. However, formal application of evidence for occurrence *in trans* with a pathogenic allele is not justified without phasing. The missense allele should therefore remain explicitly described as a strong candidate/VUS pending segregation or functional confirmation.

### Pair-level coherence

The pair provides a parsimonious autosomal-recessive explanation for the phenotype and mirrors the characteristic BUB1B allelic architecture described in the literature: a truncating allele combined with a kinase-domain missense/hypomorphic allele. The rhabdomyosarcoma signal is particularly discriminating among MVA mechanisms. No comparably coherent rare damaging pair was detected in the other screened MVA genes using the MANE-centered analysis.

## Limitations and planned confirmation

1. **Phase is unknown.** Parental genotyping is the preferred confirmation that the two variants are in trans. Their genomic separation is too large for ordinary short paired-end read phasing.
2. **No orthogonal read review yet.** The available VCF supports both calls strongly, but BAM/CRAM or FASTQ reprocessing would permit visual review and independent recall.
3. **Current first pass is SNV/indel focused.** Structural variants, copy-number changes, repeat expansions, and mosaic calls require dedicated genome-wide methods.
4. **Transcript scope is MANE-centered.** Additional clinically relevant transcripts and deep-splice/regulatory models should be included in the full reproducibility pass.
5. **Computational scores are supporting evidence only.** The novel/uncurated missense allele needs segregation and ideally BUBR1 abundance, spindle-checkpoint, or chromosome-segregation functional evidence.

## Reproducibility and provenance

- Patient-level processing: Python 3.12, pysam 0.24.0, local consequence model.
- Reference: GRCh38; MANE Select transcript NM_001211.6.
- Public annotation snapshots queried on 2026-08-28 through UCSC Genome Browser public tracks: ClinVar, gnomAD genomes v4.1.1, gnomAD exomes v4.1.1, REVEL, and AlphaMissense.
- Public predictor files were downloaded in full for the required alternate allele before local lookup.
- The public challenge evaluator specifies a clinically validated compound-heterozygous answer and accepts `PROBAND01`; its ground-truth coordinates were not accessed or inferred from private challenge resources.

## References

1. Hanks S, et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. *Nature Genetics*. 2004. https://pubmed.ncbi.nlm.nih.gov/15475955/
2. Sieben CJ, et al. BubR1 allelic effects drive phenotypic heterogeneity in mosaic-variegated aneuploidy progeria syndrome. *Journal of Clinical Investigation*. 2020. https://pmc.ncbi.nlm.nih.gov/articles/PMC6934189/
3. ClinVar VCV000533901.9, BUB1B c.2210T>G (p.Leu737Ter). https://www.ncbi.nlm.nih.gov/clinvar/variation/533901/
4. UniProtKB O60566, BUB1B/BUBR1. https://www.uniprot.org/uniprotkb/O60566/entry
5. AlphaMissense score interpretation, EMBL-EBI. https://www.ebi.ac.uk/training/online/courses/alphafold/classifying-the-effects-of-missense-variants-using-alphamissense/understanding-pathogenicity-scores-from-alphamissense/
6. UCSC Genome Browser API documentation. https://genome.ucsc.edu/goldenPath/help/api.html
7. Hackathon Track 1 evaluation implementation. https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py

## Research-use disclaimer

This is a hackathon research hypothesis, not a clinical diagnostic report or treatment recommendation. Clinical use would require confirmation in an accredited laboratory, segregation analysis, and review by qualified genetics professionals.
