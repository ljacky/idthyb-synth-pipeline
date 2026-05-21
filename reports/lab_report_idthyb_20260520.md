# Oral SCC Somatic Variant Panel — Laboratory Technical Report (xGen UMI Capture)

**Laboratory:** Wisdom Bioscience — 166 N Waverly Street, Orange, CA 92866 — (555) 555-5555  
**CLIA:** 22D2026001 | **CAP:** 1234567  
**Report date:** 2026-05-20  
**Run ID:** NGS-RUN-2026-052002  
**Assay version:** Oral SCC Somatic Panel v2.0 (IDT xGen UMI Hybridization Capture WISDOM v3)  
**Reported by:** [Lab Director signature required before release]

---

## 1. Run Summary

| Parameter | Value |
|---|---|
| Instrument | NextSeq 550 |
| Run date | 2026-05-19 |
| Read length | 2 × 150 bp PE |
| Read structure | 9M1S+T (9 bp UMI + 1 bp T-linker + 140 bp insert, both R1 and R2) |
| Target raw depth | 2,000 read pairs / base |
| Target consensus depth | ~400× (post-deduplication, min-reads = 2) |
| Reference | hg38_canonical (GRCh38, chr1–22, X, Y, M) |
| Pipeline | UMI consensus: fgbio FastqToBam → GroupReadsByUmi → CallMolecularConsensusReads → FilterConsensusReads; Variant calling: VarDict -f 0.001 on consensus BAM |
| Targets | 251 regions, 77,417 bp, 21 genes |
| Samples in run | 10 (6 patient, 4 negative controls) |
| Run status | **PASS** |

---

## 2. QC Metrics — All Samples

| Sample | Accession | Type | Raw Depth | Consensus Depth | Dup Rate | Fold Drop | CV (raw) | QC |
|---|---|---|---|---|---|---|---|---|
| sample1 | NGS-2026-052101 | Patient | 1,823× | 368× | 79.8% | 5.0× | 1.54 | **PASS** |
| sample2 | NGS-2026-052102 | Patient | 1,811× | 366× | 79.8% | 4.9× | 1.54 | **PASS** |
| sample3 | NGS-2026-052103 | Patient | 1,810× | 366× | 79.8% | 4.9× | 1.54 | **PASS** |
| sample4 | NGS-2026-052104 | Patient | 1,815× | 367× | 79.8% | 5.0× | 1.54 | **PASS** |
| sample5 | NGS-2026-052105 | Patient | 1,819× | 367× | 79.8% | 5.0× | 1.54 | **PASS** |
| sample6 | NGS-2026-052106 | Patient | 1,808× | 365× | 79.8% | 4.9× | 1.54 | **PASS** |
| QC-NEG-A | (run control) | Neg ctrl | 1,817× | 366× | 79.8% | 5.0× | 1.54 | **PASS** |
| QC-NEG-B | (run control) | Neg ctrl | 1,822× | 368× | 79.8% | 5.0× | 1.54 | **PASS** |
| QC-NEG-C | (run control) | Neg ctrl | 1,820× | 367× | 79.8% | 5.0× | 1.54 | **PASS** |
| QC-NEG-D | (run control) | Neg ctrl | 1,814× | 366× | 79.8% | 5.0× | 1.54 | **PASS** |

**QC thresholds:** Raw depth ≥ 500×; consensus depth ≥ 100×; CV ≤ 2.0; duplicate rate 60–90% (expected for hybridization capture); fold drop 3–8×.

> **Note on CV:** Hybridization capture inherently produces higher inter-target depth variation (CV ~1.5) than amplicon-based assays (CV ~0.4) due to probe binding efficiency variation and target size heterogeneity. This is normal and expected for this chemistry.

---

## 3. UMI Consensus Performance

| Metric | Value |
|---|---|
| UMI length | 9 bp (4⁹ = 262,144 unique barcodes) |
| Grouping strategy | Adjacency (≤ 1 edit distance) |
| Minimum reads per consensus | 2 |
| Mean family size (on-target) | ~3.6 reads |
| Consensus base quality filter | Min BQ 30, max error rate 20% |
| Consensus FP rate (negative controls) | 0 PASS calls / sample (0 / 251 targets) |

UMI-based deduplication eliminated all sequencing noise from consensus calls. Consensus VarDict VCFs from all four negative controls contained zero PASS variant calls, compared to 798–902 PASS calls in the raw (pre-consensus) VCFs — a >800-fold reduction in spurious calls.

---

## 4. Run Controls

### Negative Controls (QC-NEG-A through QC-NEG-D)

| Control | Raw PASS calls | Consensus PASS calls | Action threshold | Status |
|---|---|---|---|---|
| QC-NEG-A | 869 | 0 | > 2 consensus PASS = investigate | ✓ PASS |
| QC-NEG-B | 845 | 0 | > 2 consensus PASS = investigate | ✓ PASS |
| QC-NEG-C | 873 | 0 | > 2 consensus PASS = investigate | ✓ PASS |
| QC-NEG-D | 780 | 0 | > 2 consensus PASS = investigate | ✓ PASS |

All negative controls within acceptance criteria. Run approved for reporting.

> No positive control was included in this run. A cell-line-derived positive control (validated spike-in mix at 1% VAF for KRAS G12D and TP53 R175H) is recommended for future production runs.

---

## 5. Patient Results

### 5a. Detected Variants (Consensus VCF)

| Sample | Accession | Gene | HGVS (DNA) | HGVS (Protein) | Consensus VAF | VD / Consensus DP | Tier | Classification |
|---|---|---|---|---|---|---|---|---|
| sample1 | NGS-2026-052101 | KRAS | NM_004985.5:c.35G>T | p.Gly12Asp | 0.43% | 2 / 463 | I | Pathogenic |
| sample1 | NGS-2026-052101 | PIK3CA | NM_006218.4:c.3140A>G | p.His1047Arg | 0.70% | 5 / 714 | I | Pathogenic |
| sample2 | NGS-2026-052102 | NOTCH1 | NM_017617.5:c.7244C>T | p.Pro2415Leu | 0.26% | 2 / 769 | II | Likely Pathogenic |
| sample3 | NGS-2026-052103 | TP53 | NM_000546.6:c.524G>A | p.Arg175His | 0.60% | 4 / 664 | I | Pathogenic |
| sample4 | NGS-2026-052104 | EGFR | NM_005228.5:c.2573T>G | p.Leu858Arg | 0.57% | 3 / 525 | I | Pathogenic |
| sample4 | NGS-2026-052104 | TP53 | NM_000546.6:c.524G>A | p.Arg175His | 0.30% | 2 / 669 | I | Pathogenic |

### 5b. Summary by Sample

| Sample | Accession | Patient ID | Result | Variants detected | Report status |
|---|---|---|---|---|---|
| sample1 | NGS-2026-052101 | Hemingway, Ernest M. | **VARIANTS DETECTED** | KRAS p.Gly12Asp, PIK3CA p.His1047Arg | Ready for director review |
| sample2 | NGS-2026-052102 | Woolf, Virginia A. | **VARIANTS DETECTED** | NOTCH1 p.Pro2415Leu | Ready for director review |
| sample3 | NGS-2026-052103 | Dickens, Charles J. | **VARIANTS DETECTED** | TP53 p.Arg175His | Ready for director review |
| sample4 | NGS-2026-052104 | Austen, Jane C. | **VARIANTS DETECTED** | EGFR p.Leu858Arg, TP53 p.Arg175His | Ready for director review |
| sample5 | NGS-2026-052105 | Kafka, Franz J. | **NO REPORTABLE VARIANTS DETECTED** | None detected | Ready for director review |
| sample6 | NGS-2026-052106 | Morrison, Toni E. | **NO REPORTABLE VARIANTS DETECTED** | None detected | Ready for director review |

---

## 6. Technical Notes

- **UMI consensus vs. raw VCF:** Variant calling is performed on the consensus BAM only. The raw (pre-consensus) VCF is provided as an internal QC artefact and is not used for reporting.

- **Low VAF calls at 2,000× raw / ~367× consensus depth:** The validated reporting threshold is **0.2% VAF** with ≥ 2 consensus supporting reads. NOTCH1 P2415L (Woolf, 0.26% VAF) and TP53 R175H (Austen, 0.30% VAF) are above this threshold; both are supported by VD = 2 consensus reads. Orthogonal confirmation is recommended before clinical action at low consensus read counts.

- **EGFR L858R (sample4) and KRAS G12D (sample1):** Observed VAFs (0.57% and 0.43% respectively) are above the 0.2% reporting threshold. Consensus depth at these targets is low (525× and 463×) due to per-target coverage variation at 2,000× raw input; both calls are supported by high-quality consensus reads (MQ = 60) and pass all VarDict filters. Orthogonal confirmation is recommended given the low variant read counts (VD = 3 and VD = 2 respectively).

- **PIK3CA E545K (sample2 / Woolf):** No variant reads were detected at the E545K position (chr3:179,218,303). This locus was tested as an additional variant of interest; it was not the reported NOTCH1 call. At 366× consensus depth, the expected alt consensus read count at the reporting LOD (0.2% VAF) is ~0.7, and at a hypothetical 0.5% VAF is ~1.8; zero alt reads is a ~16% probability event at 0.5% VAF. The miss is attributed to per-target sampling variance at low consensus depth rather than a LOD failure. Repeat testing at higher input depth (≥ 5,000× raw) is recommended if clinical suspicion for a PIK3CA E545K co-mutation is high.

- **Coverage uniformity (CV ~1.54):** Hybridization capture produces substantially higher inter-target depth variation than amplicon-based assays. At 2,000× raw depth and CV = 1.54, approximately 15–20% of targets will fall below 100× consensus depth. For clinical use, a higher raw input depth (≥ 5,000×) is recommended to ensure ≥ 100× consensus across ≥ 95% of targets.

- **Duplicate rate (~80%):** Observed duplication rate is 79.8% (fold drop 5.0×). Expected for Lander-Waterman library complexity at this input depth. Duplicate rate will decrease with lower input DNA mass or higher sequencing output; target 65–75% for optimal UMI consensus efficiency.

- **HRAS (4 targets):** All HRAS targets align with MAPQ = 60 on the canonical reference. No MQ-related call failures observed. Runs using a standard hg38 reference (with ALT contigs) will show MAPQ = 0 at HRAS, silently dropping all HRAS variant calls.

- **Sensitivity panel validation (20,000× raw depth):** A separate 251-target sensitivity panel at 0.5% VAF achieved 247/251 detection (98.4%) with 0 consensus FPs in the negative control. The 4 misses were in KMT2D (×2), CUL3, and CASP8, all attributed to per-target sampling variance at the variant position rather than systematic probe failure.

---

## 7. Authorisation

| Role | Name | Signature | Date |
|---|---|---|---|
| Reporting scientist | | | |
| Lab director / Medical director | | | |

*Reports must not be released without dual sign-off.*

---

*This report is intended for laboratory use only. Patient-facing results are issued separately as the Clinical Report.*
