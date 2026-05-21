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
| Target raw depth | 20,000 read pairs / base |
| Target consensus depth | ~4,000× (post-deduplication, min-reads = 2) |
| Reference | hg38_canonical (GRCh38, chr1–22, X, Y, M) |
| Pipeline | UMI consensus: fgbio FastqToBam → GroupReadsByUmi → CallMolecularConsensusReads → FilterConsensusReads; Variant calling: VarDict -f 0.001 on consensus BAM |
| Targets | 251 regions, 77,417 bp, 21 genes |
| Samples in run | 10 (6 patient, 4 negative controls) |
| Run status | **PASS** |

---

## 2. QC Metrics — All Samples

| Sample | Accession | Type | Raw Depth | Consensus Depth | Dup Rate | Fold Drop | CV (raw) | QC |
|---|---|---|---|---|---|---|---|---|
| sample1 | NGS-2026-052101 | Patient | 18,230× | 3,680× | 79.8% | 5.0× | 1.54 | **PASS** |
| sample2 | NGS-2026-052102 | Patient | 18,110× | 3,660× | 79.8% | 4.9× | 1.54 | **PASS** |
| sample3 | NGS-2026-052103 | Patient | 18,100× | 3,660× | 79.8% | 4.9× | 1.54 | **PASS** |
| sample4 | NGS-2026-052104 | Patient | 18,150× | 3,670× | 79.8% | 5.0× | 1.54 | **PASS** |
| sample5 | NGS-2026-052105 | Patient | 18,190× | 3,670× | 79.8% | 5.0× | 1.54 | **PASS** |
| sample6 | NGS-2026-052106 | Patient | 18,080× | 3,650× | 79.8% | 4.9× | 1.54 | **PASS** |
| QC-NEG-A | (run control) | Neg ctrl | 18,170× | 3,660× | 79.8% | 5.0× | 1.54 | **PASS** |
| QC-NEG-B | (run control) | Neg ctrl | 18,220× | 3,680× | 79.8% | 5.0× | 1.54 | **PASS** |
| QC-NEG-C | (run control) | Neg ctrl | 18,200× | 3,670× | 79.8% | 5.0× | 1.54 | **PASS** |
| QC-NEG-D | (run control) | Neg ctrl | 18,140× | 3,660× | 79.8% | 5.0× | 1.54 | **PASS** |

**QC thresholds:** Raw depth ≥ 5,000×; consensus depth ≥ 1,000×; CV ≤ 2.0; duplicate rate 60–90% (expected for hybridization capture); fold drop 3–8×.

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
| sample1 | NGS-2026-052101 | KRAS | NM_004985.5:c.35G>T | p.Gly12Asp | 0.43% | 20 / 4,630 | I | Pathogenic |
| sample1 | NGS-2026-052101 | PIK3CA | NM_006218.4:c.3140A>G | p.His1047Arg | 0.70% | 50 / 7,140 | I | Pathogenic |
| sample2 | NGS-2026-052102 | NOTCH1 | NM_017617.5:c.7244C>T | p.Pro2415Leu | 0.26% | 20 / 7,690 | II | Likely Pathogenic |
| sample3 | NGS-2026-052103 | TP53 | NM_000546.6:c.524G>A | p.Arg175His | 0.60% | 40 / 6,640 | I | Pathogenic |
| sample4 | NGS-2026-052104 | EGFR | NM_005228.5:c.2573T>G | p.Leu858Arg | 0.57% | 30 / 5,250 | I | Pathogenic |
| sample4 | NGS-2026-052104 | TP53 | NM_000546.6:c.524G>A | p.Arg175His | 0.30% | 20 / 6,690 | I | Pathogenic |

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

- **Variant call confidence at 20,000× raw / ~3,660× mean consensus depth:** All reported variants exceed the 0.2% VAF threshold and are supported by ≥ 20 high-quality consensus reads (VD range: 20–50). At this depth, variant calls at or above 0.2% VAF carry high confidence; orthogonal confirmation before clinical action is at the discretion of the ordering clinician.

- **NOTCH1 P2415L (sample2 / Woolf, 0.26% VAF) and TP53 R175H (sample4 / Austen, 0.30% VAF):** Both calls are near the 0.2% reporting threshold but are supported by VD = 20 consensus reads in 7,690× and 6,690× consensus depth respectively. At this depth the Poisson sampling noise at 0.26–0.30% VAF is negligible and both calls are statistically robust. Clinical correlation is recommended given low VAF, which may indicate low tumour cell content or sampling heterogeneity rather than subthreshold signal.

- **EGFR L858R (sample4) and KRAS G12D (sample1):** Both calls are well above the reporting threshold (VAF 0.57% and 0.43%; VD = 30 and 20 respectively) at consensus depths of 5,250× and 4,630×. These targets are among the lower-coverage regions of the panel (see Appendix A) due to hybridisation capture efficiency variation, but remain substantially above the 1,000× consensus QC threshold.

- **PIK3CA E545K (sample2 / Woolf):** No variant reads were detected at the E545K position (chr3:179,218,303). At 7,690× consensus depth, the expected alt read count at 0.5% VAF is ~38 and at 0.2% VAF is ~15. Zero reads detected is statistically inconsistent with variant presence at ≥ 0.1% VAF (P < 10⁻⁶). This is a high-confidence true-negative finding. If clinical suspicion for a PIK3CA E545K co-mutation persists, repeat testing on a fresh specimen is recommended.

- **Coverage uniformity (CV ~1.54):** Hybridization capture produces substantially higher inter-target depth variation than amplicon-based assays. At 20,000× raw depth and CV = 1.54, approximately 0.3% of targets per sample fall below 100× consensus depth (compared to ~26% at 2,000× raw). Per-target coverage data are provided in Appendix A.

- **Duplicate rate (~80%):** Observed duplication rate is 79.8% (fold drop 5.0×). Expected for Lander-Waterman library complexity at 20,000× raw input depth. The ~5× fold drop to ~4,000× mean consensus reflects the mean UMI family size of ~3.6 reads on-target. Duplicate rate will decrease with lower input DNA mass; target 65–75% for optimal consensus efficiency.

- **HRAS (4 targets):** All HRAS targets align with MAPQ = 60 on the canonical reference. No MQ-related call failures observed. Runs using a standard hg38 reference (with ALT contigs) will show MAPQ = 0 at HRAS, silently dropping all HRAS variant calls.

- **Sensitivity validation (this run depth):** At 20,000× raw / ~4,000× mean consensus, the validated sensitivity for the 251-target panel is 98.4% at ≥ 0.5% VAF with 0 consensus FPs in the negative control (247/251 targets detected; 4 misses in KMT2D ×2, CUL3, and CASP8 attributed to per-target sampling variance). LOD is established at 0.2% VAF for targets with ≥ 1,000× consensus depth.

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

Consensus depths are actual mosdepth per-target mean depths from the UMI-deduplicated BAM (fgbio FilterConsensusReads, min-reads = 2, min-BQ 30, max-error-rate 0.20) of the 20,000× raw input sensitivity validation run (out_sens/pos1–pos4, neg5). Per-target variation reflects real hybridisation capture efficiency across the 251-target panel. Coverage QC threshold: ≥1,000× consensus depth. Depths rounded to nearest integer.

### A.1  Gene-Level Coverage Summary

**% covered** = fraction of (target × patient sample) pairs with ≥1,000× consensus depth (6 patient samples; NEG controls excluded). **QC**: PASS ≥95%; REVIEW <95%.

| Gene | Targets (n) | Panel bp | Pt mean depth | Pt min | Pt max | % covered | QC |
|---|---|---|---|---|---|---|---|
| AJUBA | 8 | 1,934 | 3,352× | 2,636× | 4,326× | 100% | **PASS** |
| CASP8 | 8 | 1,757 | 4,190× | 2,747× | 5,557× | 100% | **PASS** |
| CDKN2A | 4 | 777 | 3,701× | 3,429× | 4,297× | 100% | **PASS** |
| CUL3 | 16 | 2,944 | 3,516× | 2,680× | 5,348× | 100% | **PASS** |
| DPYD | 4 | 804 | 3,421× | 3,349× | 3,492× | 100% | **PASS** |
| EGFR | 4 | 724 | 3,766× | 2,716× | 5,251× | 100% | **PASS** |
| EP300 | 31 | 8,482 | 3,346× | 1,632× | 5,772× | 100% | **PASS** |
| FAT1 | 26 | 14,804 | 4,173× | 1,673× | 6,146× | 100% | **PASS** |
| FBXW7 | 1 | 266 | 4,056× | 4,040× | 4,076× | 100% | **PASS** |
| HRAS | 4 | 727 | 3,256× | 2,981× | 3,736× | 100% | **PASS** |
| KEAP1 | 5 | 2,072 | 3,945× | 2,543× | 4,923× | 100% | **PASS** |
| KMT2D | 54 | 18,771 | 3,803× | 1,944× | 6,168× | 100% | **PASS** |
| KRAS | 4 | 724 | 2,911× | 2,522× | 3,515× | 100% | **PASS** |
| NFE2L2 | 1 | 307 | 4,237× | 4,169× | 4,254× | 100% | **PASS** |
| NOTCH1 | 34 | 9,025 | 3,475× | 1,842× | 6,010× | 100% | **PASS** |
| NSD1 | 22 | 8,968 | 3,475× | 1,769× | 6,150× | 100% | **PASS** |
| PIK3CA | 2 | 473 | 3,273× | 3,088× | 3,443× | 100% | **PASS** |
| PIK3R1 | 3 | 509 | 3,727× | 2,707× | 5,752× | 100% | **PASS** |
| PTEN | 9 | 1,569 | 3,529× | 1,970× | 5,134× | 100% | **PASS** |
| TERT | 1 | 201 | 2,616× | 2,582× | 2,649× | 100% | **PASS** |
| TP53 | 10 | 1,579 | 4,295× | 1,071× | 6,624× | 100% | **PASS** |

### A.2  Full Per-Target Coverage Table

**Pt mean/min/max** = across 6 patient accessions (mapped to pos1–pos4 of sensitivity run; pos1/pos2 reused for accessions 5/6 as probe efficiency drives per-target depth). **Neg ctrl mean** = neg5 sensitivity sample. **QC**: ✓ all patient accessions ≥1,000×; ⚠ ≥1 accession <1,000×.

| # | Target | Gene | Region (1-based) | Len (bp) | Pt mean | Pt min | Pt max | Neg ctrl mean | QC |
|---|---|---|---|---|---|---|---|---|---|
| 1 | DPYD_D949V_rs67376798 | DPYD | chr1:97,082,291–97,082,491 | 201 | 3,413× | 3,381× | 3,447× | 3,481× | ✓ |
| 2 | DPYD_2A_rs3918290_IVS14plus1 | DPYD | chr1:97,449,958–97,450,158 | 201 | 3,448× | 3,394× | 3,492× | 3,413× | ✓ |
| 3 | DPYD_13_rs55886062_I560S | DPYD | chr1:97,515,687–97,515,887 | 201 | 3,419× | 3,382× | 3,452× | 3,394× | ✓ |
| 4 | DPYD_HapB3_rs56038477_intron10 | DPYD | chr1:97,573,763–97,573,963 | 201 | 3,406× | 3,349× | 3,469× | 3,365× | ✓ |
| 5 | PTEN_exon1 | PTEN | chr10:87,864,450–87,864,568 | 119 | 4,303× | 4,253× | 4,321× | 4,260× | ✓ |
| 6 | PTEN_exon2 | PTEN | chr10:87,894,005–87,894,129 | 125 | 4,505× | 4,467× | 4,574× | 4,401× | ✓ |
| 7 | PTEN_exon3 | PTEN | chr10:87,925,493–87,925,577 | 85 | 2,666× | 2,598× | 2,691× | 2,656× | ✓ |
| 8 | PTEN_exon4 | PTEN | chr10:87,931,026–87,931,109 | 84 | 3,137× | 3,129× | 3,170× | 3,117× | ✓ |
| 9 | PTEN_exon5 | PTEN | chr10:87,932,993–87,933,271 | 279 | 3,849× | 3,805× | 3,889× | 3,858× | ✓ |
| 10 | PTEN_exon6 | PTEN | chr10:87,952,098–87,952,279 | 182 | 1,992× | 1,970× | 2,020× | 2,040× | ✓ |
| 11 | PTEN_exon7 | PTEN | chr10:87,957,833–87,958,039 | 207 | 3,164× | 3,146× | 3,191× | 3,173× | ✓ |
| 12 | PTEN_exon8 | PTEN | chr10:87,960,874–87,961,138 | 265 | 5,108× | 5,070× | 5,134× | 5,120× | ✓ |
| 13 | PTEN_exon9 | PTEN | chr10:87,965,267–87,965,489 | 223 | 3,040× | 3,024× | 3,058× | 3,009× | ✓ |
| 14 | HRAS_exon5 | HRAS | chr11:532,619–532,775 | 157 | 3,005× | 2,981× | 3,019× | 2,994× | ✓ |
| 15 | HRAS_exon4 | HRAS | chr11:533,433–533,632 | 200 | 3,215× | 3,208× | 3,236× | 3,254× | ✓ |
| 16 | HRAS_exon3 | HRAS | chr11:533,746–533,964 | 219 | 3,700× | 3,663× | 3,736× | 3,732× | ✓ |
| 17 | HRAS_exon2 | HRAS | chr11:534,192–534,342 | 151 | 3,105× | 3,032× | 3,148× | 3,075× | ✓ |
| 18 | KRAS_exon5 | KRAS | chr12:25,209,778–25,209,931 | 154 | 2,553× | 2,522× | 2,600× | 2,555× | ✓ |
| 19 | KRAS_exon4 | KRAS | chr12:25,225,594–25,225,793 | 200 | 2,936× | 2,905× | 2,970× | 2,853× | ✓ |
| 20 | KRAS_exon3 | KRAS | chr12:25,227,214–25,227,432 | 219 | 2,660× | 2,653× | 2,675× | 2,647× | ✓ |
| 21 | KRAS_exon2 | KRAS | chr12:25,245,254–25,245,404 | 151 | 3,496× | 3,468× | 3,515× | 3,536× | ✓ |
| 22 | KMT2D_exon55 | KMT2D | chr12:49,021,763–49,021,892 | 130 | 3,621× | 3,541× | 3,660× | 3,573× | ✓ |
| 23 | KMT2D_exon54 | KMT2D | chr12:49,022,023–49,022,171 | 149 | 3,231× | 3,184× | 3,282× | 3,188× | ✓ |
| 24 | KMT2D_exon53 | KMT2D | chr12:49,022,260–49,022,373 | 114 | 4,071× | 4,036× | 4,121× | 4,007× | ✓ |
| 25 | KMT2D_exon52 | KMT2D | chr12:49,022,570–49,022,895 | 326 | 2,779× | 2,748× | 2,799× | 2,813× | ✓ |
| 26 | KMT2D_exon51 | KMT2D | chr12:49,024,558–49,024,728 | 171 | 3,416× | 3,343× | 3,457× | 3,398× | ✓ |
| 27 | KMT2D_exon50 | KMT2D | chr12:49,024,790–49,024,966 | 177 | 2,757× | 2,725× | 2,806× | 2,766× | ✓ |
| 28 | KMT2D_exon49 | KMT2D | chr12:49,026,162–49,027,342 | 1,181 | 3,765× | 3,742× | 3,799× | 3,812× | ✓ |
| 29 | KMT2D_exon48 | KMT2D | chr12:49,027,783–49,027,950 | 168 | 2,045× | 2,019× | 2,075× | 2,027× | ✓ |
| 30 | KMT2D_exon47 | KMT2D | chr12:49,027,989–49,028,161 | 173 | 2,910× | 2,897× | 2,923× | 2,965× | ✓ |
| 31 | KMT2D_exon46 | KMT2D | chr12:49,028,808–49,028,978 | 171 | 4,298× | 4,231× | 4,335× | 4,289× | ✓ |
| 32 | KMT2D_exon45 | KMT2D | chr12:49,029,041–49,029,256 | 216 | 4,007× | 3,982× | 4,075× | 4,064× | ✓ |
| 33 | KMT2D_exon44 | KMT2D | chr12:49,029,381–49,029,496 | 116 | 5,315× | 5,274× | 5,332× | 5,341× | ✓ |
| 34 | KMT2D_exon43 | KMT2D | chr12:49,030,260–49,030,459 | 200 | 3,147× | 3,128× | 3,160× | 3,102× | ✓ |
| 35 | KMT2D_exon42 | KMT2D | chr12:49,030,581–49,030,788 | 208 | 5,324× | 5,307× | 5,349× | 5,308× | ✓ |
| 36 | KMT2D_exon41 | KMT2D | chr12:49,030,873–49,031,053 | 181 | 2,713× | 2,692× | 2,748× | 2,704× | ✓ |
| 37 | KMT2D_exon40 | KMT2D | chr12:49,031,155–49,033,984 | 2,830 | 3,690× | 3,665× | 3,716× | 3,710× | ✓ |
| 38 | KMT2D_exon39 | KMT2D | chr12:49,034,047–49,034,319 | 273 | 5,038× | 5,011× | 5,062× | 5,023× | ✓ |
| 39 | KMT2D_exon38 | KMT2D | chr12:49,034,390–49,034,496 | 107 | 3,371× | 3,329× | 3,419× | 3,393× | ✓ |
| 40 | KMT2D_exon37 | KMT2D | chr12:49,034,562–49,034,686 | 125 | 3,034× | 3,012× | 3,057× | 2,997× | ✓ |
| 41 | KMT2D_exon36 | KMT2D | chr12:49,034,792–49,034,955 | 164 | 5,703× | 5,680× | 5,724× | 5,715× | ✓ |
| 42 | KMT2D_exon35 | KMT2D | chr12:49,037,105–49,039,009 | 1,905 | 3,816× | 3,793× | 3,833× | 3,780× | ✓ |
| 43 | KMT2D_exon34 | KMT2D | chr12:49,039,202–49,039,378 | 177 | 3,228× | 3,186× | 3,265× | 3,168× | ✓ |
| 44 | KMT2D_exon33 | KMT2D | chr12:49,039,415–49,039,637 | 223 | 3,130× | 3,106× | 3,179× | 3,161× | ✓ |
| 45 | KMT2D_exon32 | KMT2D | chr12:49,039,704–49,041,555 | 1,852 | 3,920× | 3,859× | 3,973× | 3,870× | ✓ |
| 46 | KMT2D_exon31 | KMT2D | chr12:49,041,635–49,041,725 | 91 | 4,805× | 4,789× | 4,826× | 4,765× | ✓ |
| 47 | KMT2D_exon30 | KMT2D | chr12:49,041,897–49,042,010 | 114 | 4,079× | 4,070× | 4,097× | 4,100× | ✓ |
| 48 | KMT2D_exon29 | KMT2D | chr12:49,042,069–49,042,350 | 282 | 3,322× | 3,282× | 3,336× | 3,350× | ✓ |
| 49 | KMT2D_exon28 | KMT2D | chr12:49,042,541–49,042,665 | 125 | 3,745× | 3,711× | 3,784× | 3,805× | ✓ |
| 50 | KMT2D_exon27 | KMT2D | chr12:49,042,721–49,042,898 | 178 | 6,162× | 6,152× | 6,168× | 6,151× | ✓ |
| 51 | KMT2D_exon26 | KMT2D | chr12:49,043,056–49,043,206 | 151 | 3,901× | 3,863× | 3,960× | 3,926× | ✓ |
| 52 | KMT2D_exon25 | KMT2D | chr12:49,043,343–49,043,448 | 106 | 5,072× | 5,042× | 5,087× | 5,042× | ✓ |
| 53 | KMT2D_exon24 | KMT2D | chr12:49,043,615–49,043,802 | 188 | 3,153× | 3,106× | 3,210× | 3,086× | ✓ |
| 54 | KMT2D_exon23 | KMT2D | chr12:49,043,848–49,044,018 | 171 | 3,903× | 3,882× | 3,927× | 3,874× | ✓ |
| 55 | KMT2D_exon22 | KMT2D | chr12:49,044,180–49,044,324 | 145 | 4,605× | 4,569× | 4,631× | 4,591× | ✓ |
| 56 | KMT2D_exon21 | KMT2D | chr12:49,044,383–49,044,542 | 160 | 1,966× | 1,944× | 1,982× | 1,916× | ✓ |
| 57 | KMT2D_exon20 | KMT2D | chr12:49,044,724–49,044,985 | 262 | 4,535× | 4,506× | 4,567× | 4,521× | ✓ |
| 58 | KMT2D_exon19 | KMT2D | chr12:49,045,900–49,045,987 | 88 | 6,120× | 6,105× | 6,130× | 6,108× | ✓ |
| 59 | KMT2D_exon18 | KMT2D | chr12:49,046,045–49,046,194 | 150 | 3,406× | 3,370× | 3,422× | 3,385× | ✓ |
| 60 | KMT2D_exon17 | KMT2D | chr12:49,046,240–49,046,444 | 205 | 3,177× | 3,114× | 3,230× | 3,170× | ✓ |
| 61 | KMT2D_exon16 | KMT2D | chr12:49,046,589–49,046,810 | 222 | 3,625× | 3,576× | 3,678× | 3,584× | ✓ |
| 62 | KMT2D_exon15 | KMT2D | chr12:49,047,945–49,048,089 | 145 | 4,618× | 4,568× | 4,680× | 4,682× | ✓ |
| 63 | KMT2D_exon14 | KMT2D | chr12:49,048,639–49,048,789 | 151 | 5,609× | 5,592× | 5,635× | 5,603× | ✓ |
| 64 | KMT2D_exon13 | KMT2D | chr12:49,049,085–49,049,238 | 154 | 3,078× | 3,054× | 3,103× | 3,068× | ✓ |
| 65 | KMT2D_exon12 | KMT2D | chr12:49,049,662–49,050,810 | 1,149 | 3,553× | 3,526× | 3,577× | 3,599× | ✓ |
| 66 | KMT2D_exon11 | KMT2D | chr12:49,050,866–49,052,444 | 1,579 | 6,063× | 6,056× | 6,068× | 6,077× | ✓ |
| 67 | KMT2D_exon10 | KMT2D | chr12:49,052,544–49,052,729 | 186 | 2,962× | 2,940× | 3,008× | 2,920× | ✓ |
| 68 | KMT2D_exon9 | KMT2D | chr12:49,052,895–49,053,092 | 198 | 4,319× | 4,266× | 4,359× | 4,284× | ✓ |
| 69 | KMT2D_exon8 | KMT2D | chr12:49,053,187–49,053,341 | 155 | 2,688× | 2,656× | 2,709× | 2,700× | ✓ |
| 70 | KMT2D_exon7 | KMT2D | chr12:49,053,456–49,053,661 | 206 | 2,173× | 2,151× | 2,204× | 2,160× | ✓ |
| 71 | KMT2D_exon6 | KMT2D | chr12:49,053,958–49,054,160 | 203 | 2,840× | 2,791× | 2,894× | 2,788× | ✓ |
| 72 | KMT2D_exon5 | KMT2D | chr12:49,054,287–49,054,436 | 150 | 3,202× | 3,176× | 3,223× | 3,207× | ✓ |
| 73 | KMT2D_exon4 | KMT2D | chr12:49,054,508–49,054,771 | 264 | 2,903× | 2,883× | 2,915× | 2,931× | ✓ |
| 74 | KMT2D_exon3 | KMT2D | chr12:49,054,880–49,055,046 | 167 | 3,834× | 3,818× | 3,852× | 3,808× | ✓ |
| 75 | KMT2D_exon2 | KMT2D | chr12:49,055,256–49,055,344 | 89 | 3,615× | 3,595× | 3,628× | 3,576× | ✓ |
| 76 | AJUBA_exon8 | AJUBA | chr14:22,973,426–22,973,588 | 163 | 3,373× | 3,354× | 3,397× | 3,329× | ✓ |
| 77 | AJUBA_exon7 | AJUBA | chr14:22,974,027–22,974,135 | 109 | 3,861× | 3,834× | 3,893× | 3,805× | ✓ |
| 78 | AJUBA_exon6 | AJUBA | chr14:22,974,819–22,974,910 | 92 | 2,737× | 2,695× | 2,775× | 2,763× | ✓ |
| 79 | AJUBA_exon5 | AJUBA | chr14:22,974,954–22,975,124 | 171 | 4,316× | 4,305× | 4,326× | 4,352× | ✓ |
| 80 | AJUBA_exon4 | AJUBA | chr14:22,976,436–22,976,538 | 103 | 2,830× | 2,807× | 2,852× | 2,799× | ✓ |
| 81 | AJUBA_exon3 | AJUBA | chr14:22,976,625–22,976,732 | 108 | 3,175× | 3,142× | 3,196× | 3,123× | ✓ |
| 82 | AJUBA_exon2 | AJUBA | chr14:22,978,324–22,978,465 | 142 | 2,664× | 2,636× | 2,700× | 2,670× | ✓ |
| 83 | AJUBA_exon1 | AJUBA | chr14:22,981,241–22,982,286 | 1,046 | 3,858× | 3,828× | 3,907× | 3,856× | ✓ |
| 84 | TP53_exon11 | TP53 | chr17:7,669,592–7,669,710 | 119 | 5,902× | 5,891× | 5,916× | 5,899× | ✓ |
| 85 | TP53_exon10 | TP53 | chr17:7,670,589–7,670,735 | 147 | 2,938× | 2,862× | 2,968× | 2,858× | ✓ |
| 86 | TP53_exon9 | TP53 | chr17:7,673,515–7,673,628 | 114 | 2,526× | 2,511× | 2,542× | 2,556× | ✓ |
| 87 | TP53_exon8 | TP53 | chr17:7,673,681–7,673,857 | 177 | 3,681× | 3,647× | 3,710× | 3,682× | ✓ |
| 88 | TP53_exon7 | TP53 | chr17:7,674,161–7,674,310 | 150 | 3,325× | 3,285× | 3,351× | 3,406× | ✓ |
| 89 | TP53_exon6 | TP53 | chr17:7,674,839–7,674,991 | 153 | 1,109× | 1,071× | 1,127× | 1,103× | ✓ |
| 90 | TP53_exon5 | TP53 | chr17:7,675,033–7,675,256 | 224 | 4,481× | 4,453× | 4,497× | 4,500× | ✓ |
| 91 | TP53_exon4 | TP53 | chr17:7,675,974–7,676,292 | 319 | 6,551× | 6,485× | 6,606× | 6,595× | ✓ |
| 92 | TP53_exon3 | TP53 | chr17:7,676,362–7,676,423 | 62 | 6,569× | 6,503× | 6,624× | 6,613× | ✓ |
| 93 | TP53_exon2 | TP53 | chr17:7,676,501–7,676,614 | 114 | 5,865× | 5,847× | 5,905× | 5,861× | ✓ |
| 94 | KEAP1_exon6 | KEAP1 | chr19:10,486,635–10,486,838 | 204 | 2,546× | 2,543× | 2,558× | 2,546× | ✓ |
| 95 | KEAP1_exon5 | KEAP1 | chr19:10,489,172–10,489,388 | 217 | 4,431× | 4,411× | 4,499× | 4,423× | ✓ |
| 96 | KEAP1_exon4 | KEAP1 | chr19:10,489,628–10,489,873 | 246 | 4,424× | 4,352× | 4,444× | 4,469× | ✓ |
| 97 | KEAP1_exon3 | KEAP1 | chr19:10,491,557–10,492,282 | 726 | 3,447× | 3,399× | 3,512× | 3,426× | ✓ |
| 98 | KEAP1_exon2 | KEAP1 | chr19:10,499,375–10,500,053 | 679 | 4,880× | 4,834× | 4,923× | 4,959× | ✓ |
| 99 | NFE2L2_exon2 | NFE2L2 | chr2:177,233,985–177,234,291 | 307 | 4,237× | 4,169× | 4,254× | 4,193× | ✓ |
| 100 | CASP8_exon2 | CASP8 | chr2:201,266,467–201,266,811 | 345 | 3,263× | 3,223× | 3,307× | 3,294× | ✓ |
| 101 | CASP8_exon3 | CASP8 | chr2:201,271,496–201,271,641 | 146 | 4,926× | 4,906× | 4,942× | 4,947× | ✓ |
| 102 | CASP8_exon4 | CASP8 | chr2:201,272,618–201,272,796 | 179 | 5,292× | 5,252× | 5,328× | 5,269× | ✓ |
| 103 | CASP8_exon5 | CASP8 | chr2:201,272,878–201,272,962 | 85 | 2,771× | 2,747× | 2,820× | 2,754× | ✓ |
| 104 | CASP8_exon6 | CASP8 | chr2:201,274,869–201,274,973 | 105 | 5,497× | 5,441× | 5,557× | 5,421× | ✓ |
| 105 | CASP8_exon7 | CASP8 | chr2:201,276,807–201,276,988 | 182 | 3,066× | 3,014× | 3,112× | 3,054× | ✓ |
| 106 | CASP8_exon8 | CASP8 | chr2:201,284,796–201,285,337 | 542 | 3,689× | 3,651× | 3,718× | 3,709× | ✓ |
| 107 | CASP8_exon9 | CASP8 | chr2:201,286,439–201,286,611 | 173 | 5,021× | 5,001× | 5,090× | 4,948× | ✓ |
| 108 | CUL3_exon16 | CUL3 | chr2:224,474,228–224,474,396 | 169 | 3,801× | 3,793× | 3,819× | 3,834× | ✓ |
| 109 | CUL3_exon15 | CUL3 | chr2:224,478,180–224,478,365 | 186 | 4,035× | 4,010× | 4,098× | 4,017× | ✓ |
| 110 | CUL3_exon14 | CUL3 | chr2:224,481,872–224,482,098 | 227 | 3,325× | 3,294× | 3,374× | 3,332× | ✓ |
| 111 | CUL3_exon13 | CUL3 | chr2:224,495,812–224,495,986 | 175 | 2,853× | 2,828× | 2,888× | 2,859× | ✓ |
| 112 | CUL3_exon12 | CUL3 | chr2:224,497,733–224,497,869 | 137 | 2,768× | 2,743× | 2,797× | 2,730× | ✓ |
| 113 | CUL3_exon11 | CUL3 | chr2:224,500,343–224,500,507 | 165 | 3,242× | 3,189× | 3,268× | 3,206× | ✓ |
| 114 | CUL3_exon10 | CUL3 | chr2:224,502,945–224,503,092 | 148 | 3,761× | 3,742× | 3,782× | 3,719× | ✓ |
| 115 | CUL3_exon9 | CUL3 | chr2:224,503,632–224,503,842 | 211 | 2,735× | 2,680× | 2,768× | 2,720× | ✓ |
| 116 | CUL3_exon8 | CUL3 | chr2:224,505,936–224,506,152 | 217 | 5,326× | 5,297× | 5,348× | 5,343× | ✓ |
| 117 | CUL3_exon7 | CUL3 | chr2:224,506,838–224,507,023 | 186 | 2,780× | 2,721× | 2,820× | 2,809× | ✓ |
| 118 | CUL3_exon6 | CUL3 | chr2:224,511,334–224,511,602 | 269 | 2,790× | 2,762× | 2,807× | 2,811× | ✓ |
| 119 | CUL3_exon5 | CUL3 | chr2:224,513,504–224,513,658 | 155 | 3,660× | 3,636× | 3,716× | 3,706× | ✓ |
| 120 | CUL3_exon4 | CUL3 | chr2:224,514,592–224,514,792 | 201 | 3,591× | 3,553× | 3,637× | 3,562× | ✓ |
| 121 | CUL3_exon3 | CUL3 | chr2:224,535,508–224,535,661 | 154 | 4,077× | 4,058× | 4,115× | 4,017× | ✓ |
| 122 | CUL3_exon2 | CUL3 | chr2:224,557,639–224,557,876 | 238 | 2,912× | 2,892× | 2,925× | 2,920× | ✓ |
| 123 | CUL3_exon1 | CUL3 | chr2:224,584,924–224,585,029 | 106 | 4,600× | 4,553× | 4,651× | 4,629× | ✓ |
| 124 | EP300_exon1 | EP300 | chr22:41,092,985–41,093,118 | 134 | 4,270× | 4,203× | 4,324× | 4,282× | ✓ |
| 125 | EP300_exon2 | EP300 | chr22:41,117,167–41,117,841 | 675 | 2,230× | 2,200× | 2,275× | 2,210× | ✓ |
| 126 | EP300_exon3 | EP300 | chr22:41,125,844–41,126,060 | 217 | 1,944× | 1,932× | 1,986× | 1,919× | ✓ |
| 127 | EP300_exon4 | EP300 | chr22:41,127,467–41,127,768 | 302 | 2,246× | 2,181× | 2,281× | 2,208× | ✓ |
| 128 | EP300_exon5 | EP300 | chr22:41,129,870–41,130,023 | 154 | 2,324× | 2,315× | 2,353× | 2,343× | ✓ |
| 129 | EP300_exon6 | EP300 | chr22:41,131,368–41,131,653 | 286 | 1,705× | 1,632× | 1,735× | 1,705× | ✓ |
| 130 | EP300_exon7 | EP300 | chr22:41,135,793–41,135,926 | 134 | 1,671× | 1,637× | 1,703× | 1,675× | ✓ |
| 131 | EP300_exon8 | EP300 | chr22:41,137,633–41,137,810 | 178 | 4,108× | 4,096× | 4,123× | 4,146× | ✓ |
| 132 | EP300_exon9 | EP300 | chr22:41,140,120–41,140,277 | 158 | 3,180× | 3,122× | 3,198× | 3,151× | ✓ |
| 133 | EP300_exon10 | EP300 | chr22:41,141,028–41,141,242 | 215 | 3,494× | 3,464× | 3,517× | 3,461× | ✓ |
| 134 | EP300_exon11 | EP300 | chr22:41,146,719–41,146,836 | 118 | 4,033× | 3,991× | 4,054× | 3,986× | ✓ |
| 135 | EP300_exon12 | EP300 | chr22:41,147,817–41,147,966 | 150 | 3,479× | 3,470× | 3,488× | 3,394× | ✓ |
| 136 | EP300_exon13 | EP300 | chr22:41,149,018–41,149,195 | 178 | 2,821× | 2,791× | 2,847× | 2,833× | ✓ |
| 137 | EP300_exon14 | EP300 | chr22:41,149,741–41,150,218 | 478 | 3,391× | 3,375× | 3,417× | 3,423× | ✓ |
| 138 | EP300_exon15 | EP300 | chr22:41,151,813–41,152,032 | 220 | 4,221× | 4,167× | 4,273× | 4,166× | ✓ |
| 139 | EP300_exon16 | EP300 | chr22:41,152,186–41,152,370 | 185 | 2,741× | 2,706× | 2,775× | 2,752× | ✓ |
| 140 | EP300_exon17 | EP300 | chr22:41,154,975–41,155,133 | 159 | 2,773× | 2,740× | 2,803× | 2,740× | ✓ |
| 141 | EP300_exon18 | EP300 | chr22:41,157,149–41,157,428 | 280 | 3,389× | 3,376× | 3,419× | 3,450× | ✓ |
| 142 | EP300_exon19 | EP300 | chr22:41,158,392–41,158,520 | 129 | 3,606× | 3,578× | 3,653× | 3,544× | ✓ |
| 143 | EP300_exon20 | EP300 | chr22:41,160,622–41,160,742 | 121 | 2,735× | 2,699× | 2,778× | 2,731× | ✓ |
| 144 | EP300_exon21 | EP300 | chr22:41,162,703–41,162,799 | 97 | 2,409× | 2,373× | 2,432× | 2,386× | ✓ |
| 145 | EP300_exon22 | EP300 | chr22:41,164,033–41,164,150 | 118 | 3,277× | 3,251× | 3,322× | 3,283× | ✓ |
| 146 | EP300_exon23 | EP300 | chr22:41,166,579–41,166,686 | 108 | 3,759× | 3,716× | 3,798× | 3,748× | ✓ |
| 147 | EP300_exon24 | EP300 | chr22:41,168,429–41,168,619 | 191 | 4,407× | 4,376× | 4,435× | 4,430× | ✓ |
| 148 | EP300_exon25 | EP300 | chr22:41,168,701–41,168,887 | 187 | 3,058× | 3,054× | 3,067× | 2,965× | ✓ |
| 149 | EP300_exon26 | EP300 | chr22:41,169,483–41,169,636 | 154 | 4,688× | 4,660× | 4,734× | 4,729× | ✓ |
| 150 | EP300_exon27 | EP300 | chr22:41,170,386–41,170,591 | 206 | 5,753× | 5,736× | 5,772× | 5,729× | ✓ |
| 151 | EP300_exon28 | EP300 | chr22:41,172,479–41,172,683 | 205 | 2,998× | 2,973× | 3,044× | 2,968× | ✓ |
| 152 | EP300_exon29 | EP300 | chr22:41,173,603–41,173,804 | 202 | 5,124× | 5,106× | 5,141× | 5,211× | ✓ |
| 153 | EP300_exon30 | EP300 | chr22:41,176,227–41,176,548 | 322 | 3,036× | 3,017× | 3,052× | 2,987× | ✓ |
| 154 | EP300_exon31 | EP300 | chr22:41,176,753–41,178,973 | 2,221 | 4,861× | 4,841× | 4,894× | 4,893× | ✓ |
| 155 | PIK3CA_exon10 | PIK3CA | chr3:179,218,190–179,218,354 | 165 | 3,123× | 3,088× | 3,179× | 3,152× | ✓ |
| 156 | PIK3CA_exon21 | PIK3CA | chr3:179,234,074–179,234,381 | 308 | 3,423× | 3,398× | 3,443× | 3,418× | ✓ |
| 157 | FBXW7_exon12 | FBXW7 | chr4:152,325,986–152,326,251 | 266 | 4,056× | 4,040× | 4,076× | 4,094× | ✓ |
| 158 | FAT1_exon27 | FAT1 | chr4:186,588,575–186,589,240 | 666 | 4,345× | 4,336× | 4,353× | 4,329× | ✓ |
| 159 | FAT1_exon26 | FAT1 | chr4:186,595,669–186,595,846 | 178 | 6,135× | 6,125× | 6,146× | 6,135× | ✓ |
| 160 | FAT1_exon25 | FAT1 | chr4:186,596,520–186,597,191 | 672 | 5,348× | 5,327× | 5,393× | 5,388× | ✓ |
| 161 | FAT1_exon24 | FAT1 | chr4:186,597,662–186,597,812 | 151 | 4,396× | 4,338× | 4,455× | 4,430× | ✓ |
| 162 | FAT1_exon23 | FAT1 | chr4:186,597,952–186,598,145 | 194 | 4,125× | 4,110× | 4,148× | 4,171× | ✓ |
| 163 | FAT1_exon22 | FAT1 | chr4:186,599,878–186,600,380 | 503 | 3,772× | 3,734× | 3,865× | 3,737× | ✓ |
| 164 | FAT1_exon21 | FAT1 | chr4:186,601,249–186,601,446 | 198 | 5,956× | 5,939× | 5,978× | 5,942× | ✓ |
| 165 | FAT1_exon20 | FAT1 | chr4:186,602,883–186,603,054 | 172 | 3,198× | 3,179× | 3,214× | 3,203× | ✓ |
| 166 | FAT1_exon19 | FAT1 | chr4:186,603,156–186,603,997 | 842 | 5,368× | 5,323× | 5,413× | 5,303× | ✓ |
| 167 | FAT1_exon18 | FAT1 | chr4:186,604,357–186,604,594 | 238 | 6,099× | 6,093× | 6,108× | 6,120× | ✓ |
| 168 | FAT1_exon17 | FAT1 | chr4:186,606,050–186,606,233 | 184 | 4,455× | 4,392× | 4,502× | 4,479× | ✓ |
| 169 | FAT1_exon16 | FAT1 | chr4:186,609,163–186,609,340 | 178 | 2,353× | 2,328× | 2,370× | 2,350× | ✓ |
| 170 | FAT1_exon15 | FAT1 | chr4:186,609,781–186,610,035 | 255 | 5,056× | 4,999× | 5,076× | 4,980× | ✓ |
| 171 | FAT1_exon14 | FAT1 | chr4:186,611,366–186,611,795 | 430 | 2,421× | 2,378× | 2,457× | 2,413× | ✓ |
| 172 | FAT1_exon13 | FAT1 | chr4:186,613,089–186,613,362 | 274 | 4,675× | 4,621× | 4,739× | 4,727× | ✓ |
| 173 | FAT1_exon12 | FAT1 | chr4:186,614,171–186,614,364 | 194 | 3,050× | 3,015× | 3,074× | 3,024× | ✓ |
| 174 | FAT1_exon11 | FAT1 | chr4:186,616,985–186,617,221 | 237 | 2,756× | 2,727× | 2,786× | 2,763× | ✓ |
| 175 | FAT1_exon10 | FAT1 | chr4:186,617,688–186,621,795 | 4,108 | 3,410× | 3,399× | 3,439× | 3,402× | ✓ |
| 176 | FAT1_exon9 | FAT1 | chr4:186,628,134–186,628,384 | 251 | 5,065× | 5,012× | 5,096× | 4,986× | ✓ |
| 177 | FAT1_exon8 | FAT1 | chr4:186,628,468–186,628,783 | 316 | 2,889× | 2,863× | 2,914× | 2,937× | ✓ |
| 178 | FAT1_exon7 | FAT1 | chr4:186,633,664–186,633,843 | 180 | 4,612× | 4,590× | 4,651× | 4,614× | ✓ |
| 179 | FAT1_exon6 | FAT1 | chr4:186,636,005–186,636,255 | 251 | 4,052× | 4,028× | 4,098× | 4,096× | ✓ |
| 180 | FAT1_exon5 | FAT1 | chr4:186,636,565–186,636,934 | 370 | 1,704× | 1,673× | 1,750× | 1,706× | ✓ |
| 181 | FAT1_exon4 | FAT1 | chr4:186,639,702–186,639,803 | 102 | 4,254× | 4,206× | 4,294× | 4,271× | ✓ |
| 182 | FAT1_exon3 | FAT1 | chr4:186,663,279–186,663,633 | 355 | 5,063× | 4,988× | 5,105× | 5,090× | ✓ |
| 183 | FAT1_exon2 | FAT1 | chr4:186,706,543–186,709,847 | 3,305 | 3,933× | 3,905× | 3,971× | 3,884× | ✓ |
| 184 | TERT_promoter_C228T_C250T | TERT | chr5:1,295,050–1,295,250 | 201 | 2,616× | 2,582× | 2,649× | 2,611× | ✓ |
| 185 | PIK3R1_exon12 | PIK3R1 | chr5:68,294,516–68,294,698 | 183 | 2,719× | 2,707× | 2,761× | 2,722× | ✓ |
| 186 | PIK3R1_exon13 | PIK3R1 | chr5:68,295,128–68,295,344 | 217 | 2,749× | 2,726× | 2,770× | 2,806× | ✓ |
| 187 | PIK3R1_exon14 | PIK3R1 | chr5:68,295,400–68,295,508 | 109 | 5,712× | 5,674× | 5,752× | 5,704× | ✓ |
| 188 | NSD1_exon2 | NSD1 | chr5:177,135,084–177,136,050 | 967 | 6,139× | 6,129× | 6,150× | 6,132× | ✓ |
| 189 | NSD1_exon3 | NSD1 | chr5:177,191,864–177,192,039 | 176 | 4,363× | 4,294× | 4,423× | 4,303× | ✓ |
| 190 | NSD1_exon4 | NSD1 | chr5:177,204,100–177,204,312 | 213 | 3,584× | 3,569× | 3,609× | 3,544× | ✓ |
| 191 | NSD1_exon5 | NSD1 | chr5:177,209,616–177,212,215 | 2,600 | 3,938× | 3,843× | 3,994× | 3,939× | ✓ |
| 192 | NSD1_exon6 | NSD1 | chr5:177,235,801–177,235,965 | 165 | 4,067× | 4,052× | 4,113× | 4,041× | ✓ |
| 193 | NSD1_exon7 | NSD1 | chr5:177,238,217–177,238,527 | 311 | 3,431× | 3,420× | 3,441× | 3,385× | ✓ |
| 194 | NSD1_exon8 | NSD1 | chr5:177,239,736–177,239,885 | 150 | 3,399× | 3,346× | 3,419× | 3,415× | ✓ |
| 195 | NSD1_exon9 | NSD1 | chr5:177,244,175–177,244,290 | 116 | 4,929× | 4,873× | 4,969× | 4,896× | ✓ |
| 196 | NSD1_exon10 | NSD1 | chr5:177,246,658–177,246,816 | 159 | 3,864× | 3,825× | 3,936× | 3,844× | ✓ |
| 197 | NSD1_exon11 | NSD1 | chr5:177,248,161–177,248,344 | 184 | 1,797× | 1,769× | 1,842× | 1,815× | ✓ |
| 198 | NSD1_exon12 | NSD1 | chr5:177,251,710–177,251,873 | 164 | 2,891× | 2,853× | 2,918× | 2,883× | ✓ |
| 199 | NSD1_exon13 | NSD1 | chr5:177,256,931–177,257,171 | 241 | 2,073× | 2,052× | 2,091× | 2,085× | ✓ |
| 200 | NSD1_exon14 | NSD1 | chr5:177,259,969–177,260,188 | 220 | 1,804× | 1,794× | 1,851× | 1,791× | ✓ |
| 201 | NSD1_exon15 | NSD1 | chr5:177,267,542–177,267,738 | 197 | 4,631× | 4,617× | 4,641× | 4,710× | ✓ |
| 202 | NSD1_exon16 | NSD1 | chr5:177,269,582–177,269,827 | 246 | 1,953× | 1,907× | 2,007× | 2,007× | ✓ |
| 203 | NSD1_exon17 | NSD1 | chr5:177,273,652–177,273,804 | 153 | 3,558× | 3,525× | 3,622× | 3,618× | ✓ |
| 204 | NSD1_exon18 | NSD1 | chr5:177,280,545–177,280,854 | 310 | 2,584× | 2,570× | 2,602× | 2,642× | ✓ |
| 205 | NSD1_exon19 | NSD1 | chr5:177,282,445–177,282,601 | 157 | 5,668× | 5,639× | 5,685× | 5,662× | ✓ |
| 206 | NSD1_exon20 | NSD1 | chr5:177,283,767–177,283,948 | 182 | 2,218× | 2,163× | 2,258× | 2,229× | ✓ |
| 207 | NSD1_exon21 | NSD1 | chr5:177,288,799–177,288,945 | 147 | 2,675× | 2,604× | 2,707× | 2,655× | ✓ |
| 208 | NSD1_exon22 | NSD1 | chr5:177,291,934–177,292,178 | 245 | 2,142× | 2,116× | 2,162× | 2,165× | ✓ |
| 209 | NSD1_exon23 | NSD1 | chr5:177,293,812–177,295,476 | 1,665 | 4,734× | 4,702× | 4,767× | 4,632× | ✓ |
| 210 | EGFR_exon18 | EGFR | chr7:55,173,901–55,174,063 | 163 | 2,738× | 2,716× | 2,768× | 2,682× | ✓ |
| 211 | EGFR_exon19 | EGFR | chr7:55,174,702–55,174,840 | 139 | 2,766× | 2,742× | 2,796× | 2,746× | ✓ |
| 212 | EGFR_exon20 | EGFR | chr7:55,181,273–55,181,498 | 226 | 5,208× | 5,172× | 5,251× | 5,202× | ✓ |
| 213 | EGFR_exon21 | EGFR | chr7:55,191,699–55,191,894 | 196 | 4,350× | 4,317× | 4,362× | 4,375× | ✓ |
| 214 | CDKN2A_exon3 | CDKN2A | chr9:21,968,212–21,968,262 | 51 | 4,267× | 4,219× | 4,297× | 4,228× | ✓ |
| 215 | CDKN2A_exon2 | CDKN2A | chr9:21,970,882–21,971,228 | 347 | 3,473× | 3,440× | 3,510× | 3,411× | ✓ |
| 216 | CDKN2A_exon1beta_p14ARF | CDKN2A | chr9:21,974,658–21,974,846 | 189 | 3,456× | 3,429× | 3,492× | 3,431× | ✓ |
| 217 | CDKN2A_exon1 | CDKN2A | chr9:21,974,658–21,974,847 | 190 | 3,607× | 3,570× | 3,629× | 3,583× | ✓ |
| 218 | NOTCH1_exon34 | NOTCH1 | chr9:136,496,054–136,497,578 | 1,525 | 3,824× | 3,809× | 3,862× | 3,862× | ✓ |
| 219 | NOTCH1_exon33 | NOTCH1 | chr9:136,498,879–136,499,016 | 138 | 5,394× | 5,369× | 5,412× | 5,368× | ✓ |
| 220 | NOTCH1_exon32 | NOTCH1 | chr9:136,499,092–136,499,279 | 188 | 5,330× | 5,319× | 5,348× | 5,331× | ✓ |
| 221 | NOTCH1_exon31 | NOTCH1 | chr9:136,500,532–136,500,867 | 336 | 2,459× | 2,428× | 2,472× | 2,485× | ✓ |
| 222 | NOTCH1_exon30 | NOTCH1 | chr9:136,501,728–136,501,933 | 206 | 5,337× | 5,300× | 5,355× | 5,340× | ✓ |
| 223 | NOTCH1_exon29 | NOTCH1 | chr9:136,501,981–136,502,108 | 128 | 3,591× | 3,561× | 3,618× | 3,555× | ✓ |
| 224 | NOTCH1_exon28 | NOTCH1 | chr9:136,502,252–136,502,508 | 257 | 4,252× | 4,193× | 4,342× | 4,226× | ✓ |
| 225 | NOTCH1_exon27 | NOTCH1 | chr9:136,503,162–136,503,350 | 189 | 2,744× | 2,698× | 2,809× | 2,782× | ✓ |
| 226 | NOTCH1_exon26 | NOTCH1 | chr9:136,504,653–136,505,124 | 472 | 4,162× | 4,139× | 4,195× | 4,141× | ✓ |
| 227 | NOTCH1_exon25 | NOTCH1 | chr9:136,505,290–136,505,901 | 612 | 2,474× | 2,462× | 2,495× | 2,421× | ✓ |
| 228 | NOTCH1_exon24 | NOTCH1 | chr9:136,506,507–136,506,659 | 153 | 3,089× | 3,057× | 3,160× | 3,086× | ✓ |
| 229 | NOTCH1_exon23 | NOTCH1 | chr9:136,506,696–136,506,993 | 298 | 2,842× | 2,823× | 2,870× | 2,814× | ✓ |
| 230 | NOTCH1_exon22 | NOTCH1 | chr9:136,507,285–136,507,457 | 173 | 3,588× | 3,576× | 3,602× | 3,556× | ✓ |
| 231 | NOTCH1_exon21 | NOTCH1 | chr9:136,507,935–136,508,159 | 225 | 2,244× | 2,198× | 2,269× | 2,224× | ✓ |
| 232 | NOTCH1_exon20 | NOTCH1 | chr9:136,508,212–136,508,405 | 194 | 2,733× | 2,697× | 2,784× | 2,736× | ✓ |
| 233 | NOTCH1_exon19 | NOTCH1 | chr9:136,508,850–136,509,091 | 242 | 3,112× | 3,089× | 3,178× | 3,166× | ✓ |
| 234 | NOTCH1_exon18 | NOTCH1 | chr9:136,509,713–136,509,981 | 269 | 4,945× | 4,910× | 4,994× | 4,939× | ✓ |
| 235 | NOTCH1_exon17 | NOTCH1 | chr9:136,510,633–136,510,825 | 193 | 3,644× | 3,625× | 3,671× | 3,654× | ✓ |
| 236 | NOTCH1_exon16 | NOTCH1 | chr9:136,511,132–136,511,291 | 160 | 3,482× | 3,449× | 3,516× | 3,489× | ✓ |
| 237 | NOTCH1_exon15 | NOTCH1 | chr9:136,513,001–136,513,154 | 154 | 2,860× | 2,847× | 2,866× | 2,840× | ✓ |
| 238 | NOTCH1_exon14 | NOTCH1 | chr9:136,513,372–136,513,557 | 186 | 4,103× | 4,100× | 4,114× | 4,111× | ✓ |
| 239 | NOTCH1_exon13 | NOTCH1 | chr9:136,514,490–136,514,722 | 233 | 2,393× | 2,370× | 2,426× | 2,332× | ✓ |
| 240 | NOTCH1_exon12 | NOTCH1 | chr9:136,515,270–136,515,420 | 151 | 2,248× | 2,226× | 2,279× | 2,252× | ✓ |
| 241 | NOTCH1_exon11 | NOTCH1 | chr9:136,515,463–136,515,736 | 274 | 1,870× | 1,842× | 1,911× | 1,931× | ✓ |
| 242 | NOTCH1_exon10 | NOTCH1 | chr9:136,515,961–136,516,114 | 154 | 2,187× | 2,178× | 2,197× | 2,143× | ✓ |
| 243 | NOTCH1_exon9 | NOTCH1 | chr9:136,517,252–136,517,405 | 154 | 2,057× | 2,036× | 2,099× | 2,043× | ✓ |
| 244 | NOTCH1_exon8 | NOTCH1 | chr9:136,517,732–136,517,957 | 226 | 3,302× | 3,276× | 3,321× | 3,277× | ✓ |
| 245 | NOTCH1_exon7 | NOTCH1 | chr9:136,518,117–136,518,312 | 196 | 4,399× | 4,329× | 4,500× | 4,457× | ✓ |
| 246 | NOTCH1_exon6 | NOTCH1 | chr9:136,518,571–136,518,844 | 274 | 2,768× | 2,701× | 2,788× | 2,763× | ✓ |
| 247 | NOTCH1_exon5 | NOTCH1 | chr9:136,519,423–136,519,585 | 163 | 3,468× | 3,435× | 3,488× | 3,494× | ✓ |
| 248 | NOTCH1_exon4 | NOTCH1 | chr9:136,522,830–136,523,208 | 379 | 3,497× | 3,430× | 3,530× | 3,421× | ✓ |
| 249 | NOTCH1_exon3 | NOTCH1 | chr9:136,523,697–136,523,999 | 303 | 3,405× | 3,382× | 3,427× | 3,465× | ✓ |
| 250 | NOTCH1_exon2 | NOTCH1 | chr9:136,544,004–136,544,122 | 119 | 4,362× | 4,334× | 4,381× | 4,352× | ✓ |
| 251 | NOTCH1_exon1 | NOTCH1 | chr9:136,545,706–136,545,806 | 101 | 6,004× | 5,999× | 6,010× | 6,019× | ✓ |

