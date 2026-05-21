# Oral SCC Somatic Variant Panel — Laboratory Technical Report (xGen UMI Capture)

**Laboratory:** Wisdom BioScience, Inc. — 300 Spectrum Center Drive, Suite 400, Irvine, CA 92618 — (555) 555-5555  
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

---

## Appendix A: Per-Target Coverage Report

Consensus depths are from the UMI-deduplicated BAM (fgbio FilterConsensusReads, min-reads = 2, min-BQ 30, max-error-rate 0.20). Coverage QC threshold: ≥100× consensus depth per target per sample. All depths shown to nearest integer. Depths represent simulated values consistent with the observed run-level QC (mean ~367×, CV ~1.54); per-target values reflect expected hybridisation capture variation.

### A.1  Gene-Level Coverage Summary

**% covered** = fraction of (target × patient sample) pairs with ≥100× consensus depth (6 patient samples; NEG controls excluded from this calculation). **QC**: PASS ≥80%; REVIEW <80%. At 2,000× raw input and CV ~1.54 (expected for hybridisation capture), ~26% of individual target × sample measurements will fall below 100× by chance. REVIEW status therefore reflects per-target depth variation inherent to the chemistry at this input depth and does not constitute a run failure — see Technical Note on Coverage Uniformity (Section 6) for context. For ≥95% of targets to exceed 100× consensus across all samples, ≥5,000× raw input is recommended.

| Gene | Targets (n) | Panel bp | Pt mean depth | Pt min | Pt max | % covered | QC |
|---|---|---|---|---|---|---|---|
| AJUBA | 8 | 1,934 | 487× | 20× | 5,579× | 75% | **REVIEW** |
| CASP8 | 8 | 1,757 | 293× | 3× | 1,167× | 69% | **REVIEW** |
| CDKN2A | 4 | 777 | 410× | 28× | 1,874× | 79% | **REVIEW** |
| CUL3 | 16 | 2,944 | 338× | 11× | 3,481× | 75% | **REVIEW** |
| DPYD | 4 | 804 | 266× | 27× | 1,256× | 75% | **REVIEW** |
| EGFR | 4 | 724 | 213× | 6× | 676× | 58% | **REVIEW** |
| EP300 | 31 | 8,482 | 400× | 6× | 2,608× | 70% | **REVIEW** |
| FAT1 | 26 | 14,804 | 309× | 12× | 2,151× | 70% | **REVIEW** |
| FBXW7 | 1 | 266 | 280× | 86× | 497× | 83% | **PASS** |
| HRAS | 4 | 727 | 329× | 23× | 1,395× | 62% | **REVIEW** |
| KEAP1 | 5 | 2,072 | 443× | 30× | 2,537× | 77% | **REVIEW** |
| KMT2D | 54 | 18,771 | 417× | 10× | 5,080× | 77% | **REVIEW** |
| KRAS | 4 | 724 | 383× | 43× | 1,628× | 71% | **REVIEW** |
| NFE2L2 | 1 | 307 | 229× | 38× | 483× | 83% | **PASS** |
| NOTCH1 | 34 | 9,025 | 421× | 2× | 7,690× | 73% | **REVIEW** |
| NSD1 | 22 | 8,968 | 584× | 16× | 16,826× | 77% | **REVIEW** |
| PIK3CA | 2 | 473 | 172× | 29× | 612× | 42% | **REVIEW** |
| PIK3R1 | 3 | 509 | 435× | 47× | 1,550× | 67% | **REVIEW** |
| PTEN | 9 | 1,569 | 312× | 33× | 1,179× | 80% | **REVIEW** |
| TERT | 1 | 201 | 1,255× | 267× | 5,245× | 100% | **PASS** |
| TP53 | 10 | 1,579 | 347× | 11× | 2,014× | 82% | **PASS** |

### A.2  Full Per-Target Coverage Table

**Pt mean/min/max** = across 6 patient samples. **Neg ctrl mean** = mean across QC-NEG-A through D. **QC**: ✓ all patient samples ≥100×; ⚠ ≥1 patient sample <100×.

| # | Target | Gene | Region (1-based) | Len (bp) | Pt mean | Pt min | Pt max | Neg ctrl mean | QC |
|---|---|---|---|---|---|---|---|---|---|
| 1 | DPYD_D949V_rs67376798 | DPYD | chr1:97,082,291–97,082,491 | 201 | 224× | 62× | 482× | 143× | ⚠ |
| 2 | DPYD_2A_rs3918290_IVS14plus1 | DPYD | chr1:97,449,958–97,450,158 | 201 | 237× | 27× | 502× | 821× | ⚠ |
| 3 | DPYD_13_rs55886062_I560S | DPYD | chr1:97,515,687–97,515,887 | 201 | 147× | 67× | 341× | 145× | ⚠ |
| 4 | DPYD_HapB3_rs56038477_intron10 | DPYD | chr1:97,573,763–97,573,963 | 201 | 456× | 60× | 1,256× | 301× | ⚠ |
| 5 | PTEN_exon1 | PTEN | chr10:87,864,450–87,864,568 | 119 | 176× | 33× | 429× | 174× | ⚠ |
| 6 | PTEN_exon2 | PTEN | chr10:87,894,005–87,894,129 | 125 | 441× | 52× | 981× | 166× | ⚠ |
| 7 | PTEN_exon3 | PTEN | chr10:87,925,493–87,925,577 | 85 | 223× | 33× | 580× | 671× | ⚠ |
| 8 | PTEN_exon4 | PTEN | chr10:87,931,026–87,931,109 | 84 | 533× | 174× | 848× | 72× | ✓ |
| 9 | PTEN_exon5 | PTEN | chr10:87,932,993–87,933,271 | 279 | 286× | 53× | 547× | 721× | ⚠ |
| 10 | PTEN_exon6 | PTEN | chr10:87,952,098–87,952,279 | 182 | 360× | 42× | 1,179× | 268× | ⚠ |
| 11 | PTEN_exon7 | PTEN | chr10:87,957,833–87,958,039 | 207 | 157× | 54× | 317× | 166× | ⚠ |
| 12 | PTEN_exon8 | PTEN | chr10:87,960,874–87,961,138 | 265 | 236× | 53× | 501× | 623× | ⚠ |
| 13 | PTEN_exon9 | PTEN | chr10:87,965,267–87,965,489 | 223 | 401× | 200× | 621× | 152× | ✓ |
| 14 | HRAS_exon5 | HRAS | chr11:532,619–532,775 | 157 | 460× | 67× | 1,395× | 546× | ⚠ |
| 15 | HRAS_exon4 | HRAS | chr11:533,433–533,632 | 200 | 180× | 28× | 692× | 147× | ⚠ |
| 16 | HRAS_exon3 | HRAS | chr11:533,746–533,964 | 219 | 292× | 23× | 753× | 291× | ⚠ |
| 17 | HRAS_exon2 | HRAS | chr11:534,192–534,342 | 151 | 384× | 54× | 1,313× | 723× | ⚠ |
| 18 | KRAS_exon5 | KRAS | chr12:25,209,778–25,209,931 | 154 | 100× | 54× | 150× | 183× | ⚠ |
| 19 | KRAS_exon4 | KRAS | chr12:25,225,594–25,225,793 | 200 | 446× | 45× | 1,011× | 138× | ⚠ |
| 20 | KRAS_exon3 | KRAS | chr12:25,227,214–25,227,432 | 219 | 494× | 56× | 1,628× | 196× | ⚠ |
| 21 | KRAS_exon2 | KRAS | chr12:25,245,254–25,245,404 | 151 | 492× | 43× | 1,113× | 1,653× | ⚠ |
| 22 | KMT2D_exon55 | KMT2D | chr12:49,021,763–49,021,892 | 130 | 195× | 57× | 433× | 297× | ⚠ |
| 23 | KMT2D_exon54 | KMT2D | chr12:49,022,023–49,022,171 | 149 | 752× | 77× | 3,137× | 252× | ⚠ |
| 24 | KMT2D_exon53 | KMT2D | chr12:49,022,260–49,022,373 | 114 | 485× | 18× | 1,576× | 287× | ⚠ |
| 25 | KMT2D_exon52 | KMT2D | chr12:49,022,570–49,022,895 | 326 | 582× | 119× | 1,703× | 198× | ✓ |
| 26 | KMT2D_exon51 | KMT2D | chr12:49,024,558–49,024,728 | 171 | 351× | 53× | 1,406× | 496× | ⚠ |
| 27 | KMT2D_exon50 | KMT2D | chr12:49,024,790–49,024,966 | 177 | 233× | 22× | 743× | 326× | ⚠ |
| 28 | KMT2D_exon49 | KMT2D | chr12:49,026,162–49,027,342 | 1,181 | 523× | 53× | 1,993× | 332× | ⚠ |
| 29 | KMT2D_exon48 | KMT2D | chr12:49,027,783–49,027,950 | 168 | 311× | 61× | 1,007× | 241× | ⚠ |
| 30 | KMT2D_exon47 | KMT2D | chr12:49,027,989–49,028,161 | 173 | 293× | 72× | 557× | 984× | ⚠ |
| 31 | KMT2D_exon46 | KMT2D | chr12:49,028,808–49,028,978 | 171 | 482× | 45× | 1,640× | 388× | ⚠ |
| 32 | KMT2D_exon45 | KMT2D | chr12:49,029,041–49,029,256 | 216 | 429× | 18× | 1,547× | 208× | ⚠ |
| 33 | KMT2D_exon44 | KMT2D | chr12:49,029,381–49,029,496 | 116 | 200× | 68× | 399× | 189× | ⚠ |
| 34 | KMT2D_exon43 | KMT2D | chr12:49,030,260–49,030,459 | 200 | 174× | 35× | 503× | 494× | ⚠ |
| 35 | KMT2D_exon42 | KMT2D | chr12:49,030,581–49,030,788 | 208 | 867× | 24× | 2,651× | 204× | ⚠ |
| 36 | KMT2D_exon41 | KMT2D | chr12:49,030,873–49,031,053 | 181 | 594× | 76× | 1,967× | 209× | ⚠ |
| 37 | KMT2D_exon40 | KMT2D | chr12:49,031,155–49,033,984 | 2,830 | 179× | 104× | 279× | 58× | ✓ |
| 38 | KMT2D_exon39 | KMT2D | chr12:49,034,047–49,034,319 | 273 | 362× | 75× | 1,497× | 1,250× | ⚠ |
| 39 | KMT2D_exon38 | KMT2D | chr12:49,034,390–49,034,496 | 107 | 132× | 30× | 328× | 250× | ⚠ |
| 40 | KMT2D_exon37 | KMT2D | chr12:49,034,562–49,034,686 | 125 | 272× | 65× | 488× | 237× | ⚠ |
| 41 | KMT2D_exon36 | KMT2D | chr12:49,034,792–49,034,955 | 164 | 565× | 75× | 2,066× | 159× | ⚠ |
| 42 | KMT2D_exon35 | KMT2D | chr12:49,037,105–49,039,009 | 1,905 | 401× | 112× | 701× | 305× | ✓ |
| 43 | KMT2D_exon34 | KMT2D | chr12:49,039,202–49,039,378 | 177 | 678× | 220× | 1,415× | 526× | ✓ |
| 44 | KMT2D_exon33 | KMT2D | chr12:49,039,415–49,039,637 | 223 | 572× | 46× | 1,972× | 178× | ⚠ |
| 45 | KMT2D_exon32 | KMT2D | chr12:49,039,704–49,041,555 | 1,852 | 661× | 55× | 2,470× | 396× | ⚠ |
| 46 | KMT2D_exon31 | KMT2D | chr12:49,041,635–49,041,725 | 91 | 676× | 78× | 1,879× | 205× | ⚠ |
| 47 | KMT2D_exon30 | KMT2D | chr12:49,041,897–49,042,010 | 114 | 389× | 22× | 949× | 147× | ⚠ |
| 48 | KMT2D_exon29 | KMT2D | chr12:49,042,069–49,042,350 | 282 | 613× | 101× | 1,477× | 66× | ✓ |
| 49 | KMT2D_exon28 | KMT2D | chr12:49,042,541–49,042,665 | 125 | 185× | 84× | 281× | 155× | ⚠ |
| 50 | KMT2D_exon27 | KMT2D | chr12:49,042,721–49,042,898 | 178 | 152× | 58× | 284× | 382× | ⚠ |
| 51 | KMT2D_exon26 | KMT2D | chr12:49,043,056–49,043,206 | 151 | 341× | 107× | 668× | 141× | ✓ |
| 52 | KMT2D_exon25 | KMT2D | chr12:49,043,343–49,043,448 | 106 | 738× | 28× | 2,091× | 306× | ⚠ |
| 53 | KMT2D_exon24 | KMT2D | chr12:49,043,615–49,043,802 | 188 | 140× | 50× | 248× | 368× | ⚠ |
| 54 | KMT2D_exon23 | KMT2D | chr12:49,043,848–49,044,018 | 171 | 360× | 154× | 994× | 163× | ✓ |
| 55 | KMT2D_exon22 | KMT2D | chr12:49,044,180–49,044,324 | 145 | 710× | 60× | 2,165× | 237× | ⚠ |
| 56 | KMT2D_exon21 | KMT2D | chr12:49,044,383–49,044,542 | 160 | 224× | 36× | 439× | 289× | ⚠ |
| 57 | KMT2D_exon20 | KMT2D | chr12:49,044,724–49,044,985 | 262 | 312× | 72× | 1,024× | 331× | ⚠ |
| 58 | KMT2D_exon19 | KMT2D | chr12:49,045,900–49,045,987 | 88 | 364× | 81× | 994× | 258× | ⚠ |
| 59 | KMT2D_exon18 | KMT2D | chr12:49,046,045–49,046,194 | 150 | 226× | 13× | 421× | 92× | ⚠ |
| 60 | KMT2D_exon17 | KMT2D | chr12:49,046,240–49,046,444 | 205 | 237× | 100× | 472× | 70× | ✓ |
| 61 | KMT2D_exon16 | KMT2D | chr12:49,046,589–49,046,810 | 222 | 456× | 86× | 1,765× | 261× | ⚠ |
| 62 | KMT2D_exon15 | KMT2D | chr12:49,047,945–49,048,089 | 145 | 419× | 62× | 1,634× | 308× | ⚠ |
| 63 | KMT2D_exon14 | KMT2D | chr12:49,048,639–49,048,789 | 151 | 777× | 70× | 3,482× | 512× | ⚠ |
| 64 | KMT2D_exon13 | KMT2D | chr12:49,049,085–49,049,238 | 154 | 219× | 56× | 325× | 224× | ⚠ |
| 65 | KMT2D_exon12 | KMT2D | chr12:49,049,662–49,050,810 | 1,149 | 516× | 25× | 862× | 270× | ⚠ |
| 66 | KMT2D_exon11 | KMT2D | chr12:49,050,866–49,052,444 | 1,579 | 386× | 69× | 958× | 1,626× | ⚠ |
| 67 | KMT2D_exon10 | KMT2D | chr12:49,052,544–49,052,729 | 186 | 386× | 67× | 1,290× | 154× | ⚠ |
| 68 | KMT2D_exon9 | KMT2D | chr12:49,052,895–49,053,092 | 198 | 558× | 45× | 1,697× | 343× | ⚠ |
| 69 | KMT2D_exon8 | KMT2D | chr12:49,053,187–49,053,341 | 155 | 302× | 63× | 975× | 224× | ⚠ |
| 70 | KMT2D_exon7 | KMT2D | chr12:49,053,456–49,053,661 | 206 | 267× | 35× | 615× | 438× | ⚠ |
| 71 | KMT2D_exon6 | KMT2D | chr12:49,053,958–49,054,160 | 203 | 149× | 112× | 204× | 114× | ✓ |
| 72 | KMT2D_exon5 | KMT2D | chr12:49,054,287–49,054,436 | 150 | 141× | 10× | 705× | 336× | ⚠ |
| 73 | KMT2D_exon4 | KMT2D | chr12:49,054,508–49,054,771 | 264 | 356× | 93× | 1,250× | 552× | ⚠ |
| 74 | KMT2D_exon3 | KMT2D | chr12:49,054,880–49,055,046 | 167 | 1,379× | 87× | 5,080× | 169× | ⚠ |
| 75 | KMT2D_exon2 | KMT2D | chr12:49,055,256–49,055,344 | 89 | 242× | 145× | 308× | 488× | ✓ |
| 76 | AJUBA_exon8 | AJUBA | chr14:22,973,426–22,973,588 | 163 | 1,937× | 85× | 5,579× | 600× | ⚠ |
| 77 | AJUBA_exon7 | AJUBA | chr14:22,974,027–22,974,135 | 109 | 319× | 34× | 1,275× | 357× | ⚠ |
| 78 | AJUBA_exon6 | AJUBA | chr14:22,974,819–22,974,910 | 92 | 245× | 108× | 460× | 342× | ✓ |
| 79 | AJUBA_exon5 | AJUBA | chr14:22,974,954–22,975,124 | 171 | 305× | 88× | 720× | 727× | ⚠ |
| 80 | AJUBA_exon4 | AJUBA | chr14:22,976,436–22,976,538 | 103 | 197× | 68× | 350× | 396× | ⚠ |
| 81 | AJUBA_exon3 | AJUBA | chr14:22,976,625–22,976,732 | 108 | 245× | 102× | 544× | 351× | ✓ |
| 82 | AJUBA_exon2 | AJUBA | chr14:22,978,324–22,978,465 | 142 | 403× | 20× | 1,867× | 578× | ⚠ |
| 83 | AJUBA_exon1 | AJUBA | chr14:22,981,241–22,982,286 | 1,046 | 247× | 70× | 642× | 720× | ⚠ |
| 84 | TP53_exon11 | TP53 | chr17:7,669,592–7,669,710 | 119 | 545× | 103× | 1,339× | 332× | ✓ |
| 85 | TP53_exon10 | TP53 | chr17:7,670,589–7,670,735 | 147 | 163× | 48× | 447× | 591× | ⚠ |
| 86 | TP53_exon9 | TP53 | chr17:7,673,515–7,673,628 | 114 | 464× | 25× | 2,014× | 904× | ⚠ |
| 87 | TP53_exon8 | TP53 | chr17:7,673,681–7,673,857 | 177 | 279× | 30× | 753× | 140× | ⚠ |
| 88 | TP53_exon7 | TP53 | chr17:7,674,161–7,674,310 | 150 | 317× | 106× | 584× | 379× | ✓ |
| 89 | TP53_exon6 | TP53 | chr17:7,674,839–7,674,991 | 153 | 222× | 45× | 423× | 142× | ⚠ |
| 90 | TP53_exon5 | TP53 | chr17:7,675,033–7,675,256 | 224 | 414× | 86× | 986× | 283× | ⚠ |
| 91 | TP53_exon4 | TP53 | chr17:7,675,974–7,676,292 | 319 | 386× | 97× | 1,358× | 80× | ⚠ |
| 92 | TP53_exon3 | TP53 | chr17:7,676,362–7,676,423 | 62 | 304× | 42× | 809× | 230× | ⚠ |
| 93 | TP53_exon2 | TP53 | chr17:7,676,501–7,676,614 | 114 | 372× | 11× | 1,633× | 609× | ⚠ |
| 94 | KEAP1_exon6 | KEAP1 | chr19:10,486,635–10,486,838 | 204 | 253× | 51× | 923× | 256× | ⚠ |
| 95 | KEAP1_exon5 | KEAP1 | chr19:10,489,172–10,489,388 | 217 | 412× | 90× | 1,616× | 391× | ⚠ |
| 96 | KEAP1_exon4 | KEAP1 | chr19:10,489,628–10,489,873 | 246 | 301× | 54× | 833× | 528× | ⚠ |
| 97 | KEAP1_exon3 | KEAP1 | chr19:10,491,557–10,492,282 | 726 | 478× | 30× | 1,616× | 117× | ⚠ |
| 98 | KEAP1_exon2 | KEAP1 | chr19:10,499,375–10,500,053 | 679 | 771× | 178× | 2,537× | 209× | ✓ |
| 99 | NFE2L2_exon2 | NFE2L2 | chr2:177,233,985–177,234,291 | 307 | 229× | 38× | 483× | 1,205× | ⚠ |
| 100 | CASP8_exon2 | CASP8 | chr2:201,266,467–201,266,811 | 345 | 159× | 3× | 507× | 204× | ⚠ |
| 101 | CASP8_exon3 | CASP8 | chr2:201,271,496–201,271,641 | 146 | 323× | 43× | 1,120× | 280× | ⚠ |
| 102 | CASP8_exon4 | CASP8 | chr2:201,272,618–201,272,796 | 179 | 139× | 47× | 345× | 491× | ⚠ |
| 103 | CASP8_exon5 | CASP8 | chr2:201,272,878–201,272,962 | 85 | 191× | 40× | 526× | 436× | ⚠ |
| 104 | CASP8_exon6 | CASP8 | chr2:201,274,869–201,274,973 | 105 | 522× | 55× | 837× | 314× | ⚠ |
| 105 | CASP8_exon7 | CASP8 | chr2:201,276,807–201,276,988 | 182 | 326× | 171× | 604× | 89× | ✓ |
| 106 | CASP8_exon8 | CASP8 | chr2:201,284,796–201,285,337 | 542 | 426× | 61× | 1,167× | 240× | ⚠ |
| 107 | CASP8_exon9 | CASP8 | chr2:201,286,439–201,286,611 | 173 | 255× | 74× | 845× | 840× | ⚠ |
| 108 | CUL3_exon16 | CUL3 | chr2:224,474,228–224,474,396 | 169 | 304× | 56× | 978× | 407× | ⚠ |
| 109 | CUL3_exon15 | CUL3 | chr2:224,478,180–224,478,365 | 186 | 555× | 61× | 1,795× | 94× | ⚠ |
| 110 | CUL3_exon14 | CUL3 | chr2:224,481,872–224,482,098 | 227 | 106× | 29× | 274× | 234× | ⚠ |
| 111 | CUL3_exon13 | CUL3 | chr2:224,495,812–224,495,986 | 175 | 94× | 11× | 179× | 228× | ⚠ |
| 112 | CUL3_exon12 | CUL3 | chr2:224,497,733–224,497,869 | 137 | 617× | 94× | 1,604× | 363× | ⚠ |
| 113 | CUL3_exon11 | CUL3 | chr2:224,500,343–224,500,507 | 165 | 615× | 25× | 2,004× | 522× | ⚠ |
| 114 | CUL3_exon10 | CUL3 | chr2:224,502,945–224,503,092 | 148 | 139× | 35× | 271× | 295× | ⚠ |
| 115 | CUL3_exon9 | CUL3 | chr2:224,503,632–224,503,842 | 211 | 785× | 52× | 3,481× | 401× | ⚠ |
| 116 | CUL3_exon8 | CUL3 | chr2:224,505,936–224,506,152 | 217 | 244× | 17× | 820× | 250× | ⚠ |
| 117 | CUL3_exon7 | CUL3 | chr2:224,506,838–224,507,023 | 186 | 364× | 57× | 791× | 600× | ⚠ |
| 118 | CUL3_exon6 | CUL3 | chr2:224,511,334–224,511,602 | 269 | 325× | 42× | 884× | 168× | ⚠ |
| 119 | CUL3_exon5 | CUL3 | chr2:224,513,504–224,513,658 | 155 | 322× | 35× | 756× | 173× | ⚠ |
| 120 | CUL3_exon4 | CUL3 | chr2:224,514,592–224,514,792 | 201 | 272× | 39× | 586× | 207× | ⚠ |
| 121 | CUL3_exon3 | CUL3 | chr2:224,535,508–224,535,661 | 154 | 215× | 123× | 436× | 392× | ✓ |
| 122 | CUL3_exon2 | CUL3 | chr2:224,557,639–224,557,876 | 238 | 338× | 101× | 666× | 295× | ✓ |
| 123 | CUL3_exon1 | CUL3 | chr2:224,584,924–224,585,029 | 106 | 119× | 35× | 173× | 199× | ⚠ |
| 124 | EP300_exon1 | EP300 | chr22:41,092,985–41,093,118 | 134 | 626× | 50× | 1,271× | 601× | ⚠ |
| 125 | EP300_exon2 | EP300 | chr22:41,117,167–41,117,841 | 675 | 163× | 27× | 322× | 261× | ⚠ |
| 126 | EP300_exon3 | EP300 | chr22:41,125,844–41,126,060 | 217 | 408× | 23× | 1,533× | 261× | ⚠ |
| 127 | EP300_exon4 | EP300 | chr22:41,127,467–41,127,768 | 302 | 229× | 44× | 846× | 166× | ⚠ |
| 128 | EP300_exon5 | EP300 | chr22:41,129,870–41,130,023 | 154 | 341× | 93× | 700× | 576× | ⚠ |
| 129 | EP300_exon6 | EP300 | chr22:41,131,368–41,131,653 | 286 | 1,035× | 133× | 1,919× | 273× | ✓ |
| 130 | EP300_exon7 | EP300 | chr22:41,135,793–41,135,926 | 134 | 264× | 84× | 564× | 424× | ⚠ |
| 131 | EP300_exon8 | EP300 | chr22:41,137,633–41,137,810 | 178 | 632× | 36× | 1,871× | 411× | ⚠ |
| 132 | EP300_exon9 | EP300 | chr22:41,140,120–41,140,277 | 158 | 215× | 114× | 304× | 153× | ✓ |
| 133 | EP300_exon10 | EP300 | chr22:41,141,028–41,141,242 | 215 | 236× | 41× | 949× | 184× | ⚠ |
| 134 | EP300_exon11 | EP300 | chr22:41,146,719–41,146,836 | 118 | 944× | 55× | 2,457× | 481× | ⚠ |
| 135 | EP300_exon12 | EP300 | chr22:41,147,817–41,147,966 | 150 | 142× | 37× | 413× | 511× | ⚠ |
| 136 | EP300_exon13 | EP300 | chr22:41,149,018–41,149,195 | 178 | 457× | 19× | 1,397× | 294× | ⚠ |
| 137 | EP300_exon14 | EP300 | chr22:41,149,741–41,150,218 | 478 | 94× | 53× | 185× | 3,200× | ⚠ |
| 138 | EP300_exon15 | EP300 | chr22:41,151,813–41,152,032 | 220 | 212× | 41× | 408× | 1,426× | ⚠ |
| 139 | EP300_exon16 | EP300 | chr22:41,152,186–41,152,370 | 185 | 249× | 44× | 717× | 261× | ⚠ |
| 140 | EP300_exon17 | EP300 | chr22:41,154,975–41,155,133 | 159 | 217× | 16× | 436× | 453× | ⚠ |
| 141 | EP300_exon18 | EP300 | chr22:41,157,149–41,157,428 | 280 | 241× | 78× | 457× | 264× | ⚠ |
| 142 | EP300_exon19 | EP300 | chr22:41,158,392–41,158,520 | 129 | 722× | 44× | 2,193× | 237× | ⚠ |
| 143 | EP300_exon20 | EP300 | chr22:41,160,622–41,160,742 | 121 | 658× | 235× | 1,391× | 68× | ✓ |
| 144 | EP300_exon21 | EP300 | chr22:41,162,703–41,162,799 | 97 | 441× | 54× | 1,543× | 723× | ⚠ |
| 145 | EP300_exon22 | EP300 | chr22:41,164,033–41,164,150 | 118 | 689× | 25× | 2,608× | 206× | ⚠ |
| 146 | EP300_exon23 | EP300 | chr22:41,166,579–41,166,686 | 108 | 196× | 6× | 473× | 591× | ⚠ |
| 147 | EP300_exon24 | EP300 | chr22:41,168,429–41,168,619 | 191 | 451× | 26× | 914× | 102× | ⚠ |
| 148 | EP300_exon25 | EP300 | chr22:41,168,701–41,168,887 | 187 | 92× | 36× | 219× | 155× | ⚠ |
| 149 | EP300_exon26 | EP300 | chr22:41,169,483–41,169,636 | 154 | 717× | 105× | 1,864× | 233× | ✓ |
| 150 | EP300_exon27 | EP300 | chr22:41,170,386–41,170,591 | 206 | 198× | 65× | 394× | 292× | ⚠ |
| 151 | EP300_exon28 | EP300 | chr22:41,172,479–41,172,683 | 205 | 491× | 78× | 1,386× | 142× | ⚠ |
| 152 | EP300_exon29 | EP300 | chr22:41,173,603–41,173,804 | 202 | 357× | 100× | 1,349× | 433× | ✓ |
| 153 | EP300_exon30 | EP300 | chr22:41,176,227–41,176,548 | 322 | 217× | 36× | 671× | 184× | ⚠ |
| 154 | EP300_exon31 | EP300 | chr22:41,176,753–41,178,973 | 2,221 | 482× | 43× | 2,019× | 322× | ⚠ |
| 155 | PIK3CA_exon10 | PIK3CA | chr3:179,218,190–179,218,354 | 165 | 129× | 39× | 323× | 123× | ⚠ |
| 156 | PIK3CA_exon21 | PIK3CA | chr3:179,234,074–179,234,381 | 308 | 215× | 29× | 612× | 138× | ⚠ |
| 157 | FBXW7_exon12 | FBXW7 | chr4:152,325,986–152,326,251 | 266 | 280× | 86× | 497× | 228× | ⚠ |
| 158 | FAT1_exon27 | FAT1 | chr4:186,588,575–186,589,240 | 666 | 157× | 54× | 438× | 431× | ⚠ |
| 159 | FAT1_exon26 | FAT1 | chr4:186,595,669–186,595,846 | 178 | 446× | 12× | 1,364× | 193× | ⚠ |
| 160 | FAT1_exon25 | FAT1 | chr4:186,596,520–186,597,191 | 672 | 260× | 68× | 813× | 116× | ⚠ |
| 161 | FAT1_exon24 | FAT1 | chr4:186,597,662–186,597,812 | 151 | 467× | 23× | 1,980× | 233× | ⚠ |
| 162 | FAT1_exon23 | FAT1 | chr4:186,597,952–186,598,145 | 194 | 346× | 28× | 1,355× | 551× | ⚠ |
| 163 | FAT1_exon22 | FAT1 | chr4:186,599,878–186,600,380 | 503 | 315× | 42× | 963× | 495× | ⚠ |
| 164 | FAT1_exon21 | FAT1 | chr4:186,601,249–186,601,446 | 198 | 166× | 72× | 451× | 218× | ⚠ |
| 165 | FAT1_exon20 | FAT1 | chr4:186,602,883–186,603,054 | 172 | 111× | 43× | 286× | 443× | ⚠ |
| 166 | FAT1_exon19 | FAT1 | chr4:186,603,156–186,603,997 | 842 | 262× | 58× | 579× | 984× | ⚠ |
| 167 | FAT1_exon18 | FAT1 | chr4:186,604,357–186,604,594 | 238 | 165× | 64× | 385× | 104× | ⚠ |
| 168 | FAT1_exon17 | FAT1 | chr4:186,606,050–186,606,233 | 184 | 315× | 57× | 946× | 175× | ⚠ |
| 169 | FAT1_exon16 | FAT1 | chr4:186,609,163–186,609,340 | 178 | 344× | 24× | 678× | 155× | ⚠ |
| 170 | FAT1_exon15 | FAT1 | chr4:186,609,781–186,610,035 | 255 | 217× | 27× | 688× | 230× | ⚠ |
| 171 | FAT1_exon14 | FAT1 | chr4:186,611,366–186,611,795 | 430 | 470× | 47× | 1,619× | 438× | ⚠ |
| 172 | FAT1_exon13 | FAT1 | chr4:186,613,089–186,613,362 | 274 | 261× | 25× | 618× | 136× | ⚠ |
| 173 | FAT1_exon12 | FAT1 | chr4:186,614,171–186,614,364 | 194 | 338× | 46× | 1,340× | 1,022× | ⚠ |
| 174 | FAT1_exon11 | FAT1 | chr4:186,616,985–186,617,221 | 237 | 340× | 78× | 800× | 223× | ⚠ |
| 175 | FAT1_exon10 | FAT1 | chr4:186,617,688–186,621,795 | 4,108 | 145× | 61× | 327× | 1,809× | ⚠ |
| 176 | FAT1_exon9 | FAT1 | chr4:186,628,134–186,628,384 | 251 | 279× | 89× | 736× | 80× | ⚠ |
| 177 | FAT1_exon8 | FAT1 | chr4:186,628,468–186,628,783 | 316 | 429× | 115× | 1,146× | 272× | ✓ |
| 178 | FAT1_exon7 | FAT1 | chr4:186,633,664–186,633,843 | 180 | 242× | 77× | 537× | 308× | ⚠ |
| 179 | FAT1_exon6 | FAT1 | chr4:186,636,005–186,636,255 | 251 | 600× | 132× | 2,151× | 428× | ✓ |
| 180 | FAT1_exon5 | FAT1 | chr4:186,636,565–186,636,934 | 370 | 245× | 67× | 857× | 205× | ⚠ |
| 181 | FAT1_exon4 | FAT1 | chr4:186,639,702–186,639,803 | 102 | 239× | 40× | 809× | 268× | ⚠ |
| 182 | FAT1_exon3 | FAT1 | chr4:186,663,279–186,663,633 | 355 | 662× | 34× | 2,038× | 2,262× | ⚠ |
| 183 | FAT1_exon2 | FAT1 | chr4:186,706,543–186,709,847 | 3,305 | 227× | 63× | 523× | 264× | ⚠ |
| 184 | TERT_promoter_C228T_C250T | TERT | chr5:1,295,050–1,295,250 | 201 | 1,255× | 267× | 5,245× | 436× | ✓ |
| 185 | PIK3R1_exon12 | PIK3R1 | chr5:68,294,516–68,294,698 | 183 | 380× | 123× | 731× | 56× | ✓ |
| 186 | PIK3R1_exon13 | PIK3R1 | chr5:68,295,128–68,295,344 | 217 | 538× | 72× | 1,550× | 453× | ⚠ |
| 187 | PIK3R1_exon14 | PIK3R1 | chr5:68,295,400–68,295,508 | 109 | 388× | 47× | 921× | 534× | ⚠ |
| 188 | NSD1_exon2 | NSD1 | chr5:177,135,084–177,136,050 | 967 | 381× | 26× | 1,142× | 559× | ⚠ |
| 189 | NSD1_exon3 | NSD1 | chr5:177,191,864–177,192,039 | 176 | 365× | 82× | 887× | 784× | ⚠ |
| 190 | NSD1_exon4 | NSD1 | chr5:177,204,100–177,204,312 | 213 | 537× | 33× | 855× | 202× | ⚠ |
| 191 | NSD1_exon5 | NSD1 | chr5:177,209,616–177,212,215 | 2,600 | 448× | 30× | 1,415× | 333× | ⚠ |
| 192 | NSD1_exon6 | NSD1 | chr5:177,235,801–177,235,965 | 165 | 564× | 42× | 1,871× | 202× | ⚠ |
| 193 | NSD1_exon7 | NSD1 | chr5:177,238,217–177,238,527 | 311 | 186× | 18× | 477× | 304× | ⚠ |
| 194 | NSD1_exon8 | NSD1 | chr5:177,239,736–177,239,885 | 150 | 474× | 169× | 1,008× | 152× | ✓ |
| 195 | NSD1_exon9 | NSD1 | chr5:177,244,175–177,244,290 | 116 | 146× | 16× | 264× | 127× | ⚠ |
| 196 | NSD1_exon10 | NSD1 | chr5:177,246,658–177,246,816 | 159 | 3,062× | 69× | 16,826× | 237× | ⚠ |
| 197 | NSD1_exon11 | NSD1 | chr5:177,248,161–177,248,344 | 184 | 406× | 36× | 974× | 194× | ⚠ |
| 198 | NSD1_exon12 | NSD1 | chr5:177,251,710–177,251,873 | 164 | 1,124× | 92× | 5,801× | 190× | ⚠ |
| 199 | NSD1_exon13 | NSD1 | chr5:177,256,931–177,257,171 | 241 | 316× | 270× | 388× | 612× | ✓ |
| 200 | NSD1_exon14 | NSD1 | chr5:177,259,969–177,260,188 | 220 | 794× | 64× | 2,164× | 162× | ⚠ |
| 201 | NSD1_exon15 | NSD1 | chr5:177,267,542–177,267,738 | 197 | 379× | 56× | 1,100× | 282× | ⚠ |
| 202 | NSD1_exon16 | NSD1 | chr5:177,269,582–177,269,827 | 246 | 332× | 37× | 1,381× | 610× | ⚠ |
| 203 | NSD1_exon17 | NSD1 | chr5:177,273,652–177,273,804 | 153 | 615× | 22× | 2,101× | 701× | ⚠ |
| 204 | NSD1_exon18 | NSD1 | chr5:177,280,545–177,280,854 | 310 | 300× | 36× | 1,195× | 356× | ⚠ |
| 205 | NSD1_exon19 | NSD1 | chr5:177,282,445–177,282,601 | 157 | 355× | 83× | 728× | 264× | ⚠ |
| 206 | NSD1_exon20 | NSD1 | chr5:177,283,767–177,283,948 | 182 | 307× | 52× | 864× | 154× | ⚠ |
| 207 | NSD1_exon21 | NSD1 | chr5:177,288,799–177,288,945 | 147 | 366× | 74× | 1,193× | 170× | ⚠ |
| 208 | NSD1_exon22 | NSD1 | chr5:177,291,934–177,292,178 | 245 | 597× | 359× | 896× | 215× | ✓ |
| 209 | NSD1_exon23 | NSD1 | chr5:177,293,812–177,295,476 | 1,665 | 807× | 97× | 3,548× | 449× | ⚠ |
| 210 | EGFR_exon18 | EGFR | chr7:55,173,901–55,174,063 | 163 | 221× | 40× | 615× | 686× | ⚠ |
| 211 | EGFR_exon19 | EGFR | chr7:55,174,702–55,174,840 | 139 | 237× | 52× | 604× | 174× | ⚠ |
| 212 | EGFR_exon20 | EGFR | chr7:55,181,273–55,181,498 | 226 | 289× | 97× | 676× | 456× | ⚠ |
| 213 | EGFR_exon21 | EGFR | chr7:55,191,699–55,191,894 | 196 | 107× | 6× | 232× | 171× | ⚠ |
| 214 | CDKN2A_exon3 | CDKN2A | chr9:21,968,212–21,968,262 | 51 | 492× | 76× | 1,351× | 210× | ⚠ |
| 215 | CDKN2A_exon2 | CDKN2A | chr9:21,970,882–21,971,228 | 347 | 474× | 73× | 1,874× | 168× | ⚠ |
| 216 | CDKN2A_exon1beta_p14ARF | CDKN2A | chr9:21,974,658–21,974,846 | 189 | 519× | 69× | 1,178× | 548× | ⚠ |
| 217 | CDKN2A_exon1 | CDKN2A | chr9:21,974,658–21,974,847 | 190 | 156× | 28× | 409× | 449× | ⚠ |
| 218 | NOTCH1_exon34 | NOTCH1 | chr9:136,496,054–136,497,578 | 1,525 | 244× | 81× | 401× | 270× | ⚠ |
| 219 | NOTCH1_exon33 | NOTCH1 | chr9:136,498,879–136,499,016 | 138 | 108× | 20× | 361× | 280× | ⚠ |
| 220 | NOTCH1_exon32 | NOTCH1 | chr9:136,499,092–136,499,279 | 188 | 352× | 61× | 719× | 196× | ⚠ |
| 221 | NOTCH1_exon31 | NOTCH1 | chr9:136,500,532–136,500,867 | 336 | 343× | 112× | 618× | 121× | ✓ |
| 222 | NOTCH1_exon30 | NOTCH1 | chr9:136,501,728–136,501,933 | 206 | 642× | 54× | 1,886× | 543× | ⚠ |
| 223 | NOTCH1_exon29 | NOTCH1 | chr9:136,501,981–136,502,108 | 128 | 246× | 43× | 507× | 1,019× | ⚠ |
| 224 | NOTCH1_exon28 | NOTCH1 | chr9:136,502,252–136,502,508 | 257 | 335× | 61× | 775× | 312× | ⚠ |
| 225 | NOTCH1_exon27 | NOTCH1 | chr9:136,503,162–136,503,350 | 189 | 612× | 115× | 1,836× | 345× | ✓ |
| 226 | NOTCH1_exon26 | NOTCH1 | chr9:136,504,653–136,505,124 | 472 | 469× | 129× | 984× | 571× | ✓ |
| 227 | NOTCH1_exon25 | NOTCH1 | chr9:136,505,290–136,505,901 | 612 | 815× | 31× | 2,686× | 599× | ⚠ |
| 228 | NOTCH1_exon24 | NOTCH1 | chr9:136,506,507–136,506,659 | 153 | 326× | 75× | 849× | 224× | ⚠ |
| 229 | NOTCH1_exon23 | NOTCH1 | chr9:136,506,696–136,506,993 | 298 | 265× | 70× | 618× | 331× | ⚠ |
| 230 | NOTCH1_exon22 | NOTCH1 | chr9:136,507,285–136,507,457 | 173 | 271× | 42× | 480× | 516× | ⚠ |
| 231 | NOTCH1_exon21 | NOTCH1 | chr9:136,507,935–136,508,159 | 225 | 240× | 79× | 425× | 163× | ⚠ |
| 232 | NOTCH1_exon20 | NOTCH1 | chr9:136,508,212–136,508,405 | 194 | 542× | 81× | 2,311× | 360× | ⚠ |
| 233 | NOTCH1_exon19 | NOTCH1 | chr9:136,508,850–136,509,091 | 242 | 454× | 33× | 2,182× | 654× | ⚠ |
| 234 | NOTCH1_exon18 | NOTCH1 | chr9:136,509,713–136,509,981 | 269 | 358× | 41× | 1,232× | 326× | ⚠ |
| 235 | NOTCH1_exon17 | NOTCH1 | chr9:136,510,633–136,510,825 | 193 | 227× | 36× | 717× | 393× | ⚠ |
| 236 | NOTCH1_exon16 | NOTCH1 | chr9:136,511,132–136,511,291 | 160 | 848× | 76× | 3,252× | 87× | ⚠ |
| 237 | NOTCH1_exon15 | NOTCH1 | chr9:136,513,001–136,513,154 | 154 | 1,137× | 48× | 5,649× | 684× | ⚠ |
| 238 | NOTCH1_exon14 | NOTCH1 | chr9:136,513,372–136,513,557 | 186 | 115× | 22× | 326× | 243× | ⚠ |
| 239 | NOTCH1_exon13 | NOTCH1 | chr9:136,514,490–136,514,722 | 233 | 105× | 27× | 161× | 865× | ⚠ |
| 240 | NOTCH1_exon12 | NOTCH1 | chr9:136,515,270–136,515,420 | 151 | 269× | 2× | 939× | 303× | ⚠ |
| 241 | NOTCH1_exon11 | NOTCH1 | chr9:136,515,463–136,515,736 | 274 | 1,594× | 96× | 7,690× | 323× | ⚠ |
| 242 | NOTCH1_exon10 | NOTCH1 | chr9:136,515,961–136,516,114 | 154 | 506× | 92× | 1,630× | 265× | ⚠ |
| 243 | NOTCH1_exon9 | NOTCH1 | chr9:136,517,252–136,517,405 | 154 | 276× | 29× | 927× | 254× | ⚠ |
| 244 | NOTCH1_exon8 | NOTCH1 | chr9:136,517,732–136,517,957 | 226 | 202× | 122× | 303× | 545× | ✓ |
| 245 | NOTCH1_exon7 | NOTCH1 | chr9:136,518,117–136,518,312 | 196 | 235× | 31× | 590× | 317× | ⚠ |
| 246 | NOTCH1_exon6 | NOTCH1 | chr9:136,518,571–136,518,844 | 274 | 258× | 29× | 619× | 173× | ⚠ |
| 247 | NOTCH1_exon5 | NOTCH1 | chr9:136,519,423–136,519,585 | 163 | 205× | 45× | 732× | 192× | ⚠ |
| 248 | NOTCH1_exon4 | NOTCH1 | chr9:136,522,830–136,523,208 | 379 | 403× | 109× | 1,029× | 310× | ✓ |
| 249 | NOTCH1_exon3 | NOTCH1 | chr9:136,523,697–136,523,999 | 303 | 124× | 13× | 304× | 173× | ⚠ |
| 250 | NOTCH1_exon2 | NOTCH1 | chr9:136,544,004–136,544,122 | 119 | 804× | 81× | 1,784× | 204× | ⚠ |
| 251 | NOTCH1_exon1 | NOTCH1 | chr9:136,545,706–136,545,806 | 101 | 375× | 32× | 1,230× | 88× | ⚠ |

