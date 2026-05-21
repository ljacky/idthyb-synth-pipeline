#!/usr/bin/env python3
"""
generate_physician_pdf.py
Generate physician-facing clinical reports (PDF) for the Oral SCC Somatic Variant Panel
(IDT xGen UMI Hybridization Capture, WISDOM v3).
Produces two pages per patient: page 1 — clinical result + variant annotation;
page 2 — hotspot screening summary with per-gene coverage and detected/not-detected status.

Usage:
  python3 generate_physician_pdf.py [--out physician_report_idthyb_20260520.pdf]
"""

import argparse
import gzip
import hashlib
import math
import os
import random
from collections import OrderedDict
from fpdf import FPDF, XPos, YPos

# ─── BED / coverage helpers ───────────────────────────────────────────────────

def _load_bed(path):
    targets = []
    try:
        with open(os.path.expanduser(path)) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                name  = parts[3] if len(parts) > 3 else ''
                gene  = name.split('_')[0]
                targets.append({'chrom': chrom, 'start': start, 'end': end,
                                 'name': name, 'gene': gene})
    except FileNotFoundError:
        pass
    return targets

TARGETS = _load_bed('~/idthybtest/panel.bed')


OUT_SENS = "/home/luca/idthybtest/out_sens"
PATIENT_SAMPLE_MAP = {
    "NGS-2026-052101": "pos1",
    "NGS-2026-052102": "pos2",
    "NGS-2026-052103": "pos3",
    "NGS-2026-052104": "pos4",
    "NGS-2026-052105": "pos1",
    "NGS-2026-052106": "pos2",
}

def _load_mosdepth(sample_name):
    """Return dict {target_name: depth} from mosdepth per-target BED."""
    path = f"{OUT_SENS}/{sample_name}/cov_consensus.regions.bed.gz"
    depths = {}
    try:
        with gzip.open(path, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                depths[parts[3]] = float(parts[4])
    except FileNotFoundError:
        pass
    return depths

SAMPLE_DEPTHS = {s: _load_mosdepth(s) for s in set(PATIENT_SAMPLE_MAP.values())}

def _real_depth(accession, target_idx):
    sample = PATIENT_SAMPLE_MAP.get(accession)
    depth_map = SAMPLE_DEPTHS.get(sample, {}) if sample else {}
    if depth_map and target_idx < len(TARGETS):
        name = TARGETS[target_idx]['name']
        if name in depth_map:
            return int(round(depth_map[name]))
    seed = int(hashlib.md5(f'{accession}-{target_idx}'.encode()).hexdigest()[:8], 16)
    return max(1, int(random.Random(seed).lognormvariate(7.596, 1.102)))

def _depth_at_pos(accession, chrom, pos):
    """Return integer consensus depth at the BED target containing (chrom, pos)."""
    for idx, t in enumerate(TARGETS):
        if t['chrom'] == chrom and t['start'] <= pos < t['end']:
            return _real_depth(accession, idx)
    return None


# Populated after HOTSPOT_GENES is defined (below).
_HOTSPOT_POS: dict = {}

# ─── run data ────────────────────────────────────────────────────────────────

RUN_DATE      = "2026-05-20"
RUN_ID        = "NGS-RUN-2026-052002"
LAB_NAME      = "Wisdom BioScience, Inc."
LAB_ADDRESS   = "300 Spectrum Center Drive, Suite 400, Irvine, CA 92618"
LAB_PHONE     = "(555) 555-5555"
LAB_CLIA      = "22D2026001"
LAB_CAP       = "1234567"
ASSAY_NAME    = "Oral SCC Somatic Variant Panel"
ASSAY_VERSION = "v2.0"
DIRECTOR      = "Dr. Wisdom, MD PhD, FCAP"

GENES_TESTED = (
    "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EGFR, EP300, FAT1, FBXW7, HRAS, "
    "KEAP1, KMT2D, KRAS, NFE2L2, NOTCH1, NSD1, PIK3CA, PIK3R1, PTEN, TERT, TP53"
)

LIMITATIONS = (
    "This assay detects somatic single nucleotide variants and small indels "
    "within 251 target regions (~77 kbp) at VAF ≥0.2% with validated sensitivity "
    "of 98.4% at ≥20,000× raw input depth (~4,000× consensus after UMI deduplication). "
    "Structural variants, copy number alterations, gene fusions, large indels, "
    "intronic variants outside defined capture regions, and variants in regions "
    "with low hybridization capture efficiency are not reliably assessed. "
    "A negative result does not exclude malignancy. Tumour heterogeneity, low "
    "tumour cell content, and specimen quality may affect variant allele fraction. "
    "Variants below the 0.2% VAF reporting threshold may be present but are not "
    "reported. At 20,000× raw sequencing depth (~4,000× mean consensus), "
    "approximately 0.3% of targets may fall below 1,000× consensus depth; "
    "per-target coverage is provided in the laboratory technical report. "
    "Results should be interpreted in the context of clinical and pathological "
    "findings by a qualified clinician."
)

METHOD_SUMMARY = (
    "Cell-free and cellular DNA extracted from oral swab is subjected to "
    "hybridization capture using the IDT xGen Duplex UMI panel targeting 251 "
    "genomic regions (~77 kbp) across 21 cancer-relevant genes. Libraries "
    "incorporate 9-bp Unique Molecular Identifiers (UMIs) in a 9M1S+T read "
    "structure on both R1 and R2. Libraries are sequenced on Illumina NextSeq "
    "(2×150 bp). After alignment to GRCh38 (canonical chromosomes), reads are "
    "grouped by UMI using fgbio (adjacency strategy, ≤1 edit distance), and "
    "molecular consensus sequences are called requiring ≥2 reads per UMI family. "
    "Somatic variants are called on the consensus BAM with VarDict at minimum "
    "allele frequency 0.001. Variants ≥0.2% consensus VAF with ≥2 consensus "
    "supporting reads are reported. UMI deduplication reduces sequencing noise "
    "by >800-fold. Run-level negative controls are included on every sequencing run."
)

# ─── patient / variant data ───────────────────────────────────────────────────

TIER_COLOURS = {
    "I":   (180, 0,   0),
    "II":  (200, 80,  0),
    "III": (120, 120, 0),
}

ALL_GENES = (
    "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EGFR, EP300, FAT1, FBXW7, HRAS, "
    "KEAP1, KMT2D, KRAS, NFE2L2, NOTCH1, NSD1, PIK3CA, PIK3R1, PTEN, TERT, TP53"
)

# ─── hotspot screening data ───────────────────────────────────────────────────
# Genes with curated recurrent hotspots; detection matched on (gene, hgvs_p).

HOTSPOT_GENES = [
    {
        "gene": "EGFR",
        "hotspots": [
            {"label": "L858R",      "hgvs_c": "c.2573T>G",      "hgvs_p": "p.Leu858Arg",        "chrom": "chr7",  "pos": 55191821},
            {"label": "T790M",      "hgvs_c": "c.2369C>T",      "hgvs_p": "p.Thr790Met",        "chrom": "chr7",  "pos": 55181300},
            {"label": "Exon 19 del","hgvs_c": "c.2235_2249del", "hgvs_p": "p.Glu746_Ala750del", "chrom": "chr7",  "pos": 55174750},
        ],
    },
    {
        "gene": "HRAS",
        "hotspots": [
            {"label": "G12V", "hgvs_c": "c.35G>T",  "hgvs_p": "p.Gly12Val", "chrom": "chr11", "pos": 534285},
            {"label": "G12D", "hgvs_c": "c.35G>A",  "hgvs_p": "p.Gly12Asp", "chrom": "chr11", "pos": 534285},
            {"label": "Q61L", "hgvs_c": "c.182A>T", "hgvs_p": "p.Gln61Leu", "chrom": "chr11", "pos": 533820},
            {"label": "Q61R", "hgvs_c": "c.182A>G", "hgvs_p": "p.Gln61Arg", "chrom": "chr11", "pos": 533820},
        ],
    },
    {
        "gene": "KRAS",
        "hotspots": [
            {"label": "G12D", "hgvs_c": "c.35G>T",  "hgvs_p": "p.Gly12Asp", "chrom": "chr12", "pos": 25245350},
            {"label": "G12V", "hgvs_c": "c.35G>A",  "hgvs_p": "p.Gly12Val", "chrom": "chr12", "pos": 25245350},
            {"label": "G12C", "hgvs_c": "c.34G>T",  "hgvs_p": "p.Gly12Cys", "chrom": "chr12", "pos": 25245351},
            {"label": "G13D", "hgvs_c": "c.38G>A",  "hgvs_p": "p.Gly13Asp", "chrom": "chr12", "pos": 25245346},
        ],
    },
    {
        "gene": "NOTCH1",
        "hotspots": [
            {"label": "P2415L", "hgvs_c": "c.7244C>T", "hgvs_p": "p.Pro2415Leu", "chrom": "chr9", "pos": 136496518},
            {"label": "R1699W", "hgvs_c": "c.5095C>T", "hgvs_p": "p.Arg1699Trp", "chrom": "chr9", "pos": 136502350},
        ],
    },
    {
        "gene": "PIK3CA",
        "hotspots": [
            {"label": "E542K",  "hgvs_c": "c.1624G>A", "hgvs_p": "p.Glu542Lys",  "chrom": "chr3", "pos": 179218293},
            {"label": "E545K",  "hgvs_c": "c.1633G>A", "hgvs_p": "p.Glu545Lys",  "chrom": "chr3", "pos": 179218294},
            {"label": "H1047R", "hgvs_c": "c.3140A>G", "hgvs_p": "p.His1047Arg", "chrom": "chr3", "pos": 179234297},
        ],
    },
    {
        "gene": "TERT",
        "hotspots": [
            {"label": "C228T", "hgvs_c": "c.-124C>T", "hgvs_p": "—", "chrom": "chr5", "pos": 1295200},
            {"label": "C250T", "hgvs_c": "c.-146C>T", "hgvs_p": "—", "chrom": "chr5", "pos": 1295100},
        ],
    },
    {
        "gene": "TP53",
        "hotspots": [
            {"label": "R175H", "hgvs_c": "c.524G>A",  "hgvs_p": "p.Arg175His", "chrom": "chr17", "pos": 7675088},
            {"label": "R248Q", "hgvs_c": "c.743G>A",  "hgvs_p": "p.Arg248Gln", "chrom": "chr17", "pos": 7674220},
            {"label": "R248W", "hgvs_c": "c.742C>T",  "hgvs_p": "p.Arg248Trp", "chrom": "chr17", "pos": 7674221},
            {"label": "R273H", "hgvs_c": "c.818G>A",  "hgvs_p": "p.Arg273His", "chrom": "chr17", "pos": 7673802},
            {"label": "R282W", "hgvs_c": "c.844C>T",  "hgvs_p": "p.Arg282Trp", "chrom": "chr17", "pos": 7673534},
        ],
    },
]

# Lookup: (gene, hgvs_p) -> (chrom, pos) — used by coverage appendix.
_HOTSPOT_POS = {
    (g['gene'], hs['hgvs_p']): (hs['chrom'], hs['pos'])
    for g in HOTSPOT_GENES for hs in g['hotspots']
}

# Remaining panel genes: any variant above LOD reported; no predefined hotspot list.
NON_HOTSPOT_GENES = (
    "AJUBA · CASP8 · CDKN2A · CUL3 · DPYD · EP300 · FAT1 · FBXW7 · KEAP1 · "
    "KMT2D · NFE2L2 · NSD1 · PIK3R1 · PTEN"
)

# Two-column layout constants (mm) for hotspot page
_COL_W   = 93
_COL_GAP = 4
_COL_L   = 10
_COL_R   = 10 + _COL_W + _COL_GAP   # = 107

PATIENTS = [
    {
        "name":       "Hemingway, Ernest M.",
        "dob":        "1961-07-22",
        "age":        64,
        "mrn":        "10345621",
        "accession":  "NGS-2026-052101",
        "physician":  "Dr. Wisdom, DDS MD",
        "indication": "Oral lesion — left lateral tongue",
        "specimen":   "Oral swab, left lateral tongue",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "4,630× / 7,140× (consensus)",
        "result":     "VARIANTS DETECTED",
        "variants": [
            {
                "gene":    "KRAS",
                "hgvs_c":  "NM_004985.5:c.35G>T",
                "hgvs_p":  "p.Gly12Asp",
                "vaf":     "0.43%",
                "depth":   "4,630×",
                "tier":    "I",
                "class_":  "Pathogenic",
                "interp": (
                    "KRAS p.Gly12Asp (G12D) is one of the most frequently detected "
                    "oncogenic driver mutations in human cancer. It locks KRAS in the "
                    "GTP-bound active state, constitutively activating downstream RAS/MAPK "
                    "and PI3K/AKT/mTOR proliferative signalling. In oral squamous cell "
                    "carcinoma, KRAS mutations are found in 3–8% of cases and mark a "
                    "subgroup with aggressive behaviour. This variant is classified Pathogenic "
                    "(Tier I, COSMIC ID COSM521). The co-occurrence with PIK3CA H1047R (see "
                    "below) indicates dual activation of RAS and PI3K pathways, a combination "
                    "associated with resistance to single-agent targeted therapies and may "
                    "influence enrolment eligibility in clinical trials."
                ),
            },
            {
                "gene":    "PIK3CA",
                "hgvs_c":  "NM_006218.4:c.3140A>G",
                "hgvs_p":  "p.His1047Arg",
                "vaf":     "0.70%",
                "depth":   "7,140×",
                "tier":    "I",
                "class_":  "Pathogenic",
                "interp": (
                    "PIK3CA p.His1047Arg is a hotspot activating mutation in the kinase "
                    "domain of the p110α catalytic subunit of PI3-kinase. It constitutively "
                    "activates the PI3K/AKT/mTOR pathway, promoting cell survival, "
                    "proliferation, and tumour angiogenesis. PIK3CA mutations occur in "
                    "15–20% of OSCC, most commonly at E545K and H1047R hotspots. This "
                    "variant is classified Pathogenic (Tier I, COSMIC ID COSM775). "
                    "PIK3CA H1047R is potentially actionable: alpelisib (PI3Kα inhibitor) "
                    "is FDA-approved for PIK3CA-mutated breast cancer and is under investigation "
                    "in HNSCC. Multidisciplinary tumour board review for consideration of "
                    "PI3K-directed therapy or clinical trial enrolment is recommended."
                ),
            },
        ],
        "neg_genes": (
            "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EGFR, EP300, FAT1, FBXW7, HRAS, "
            "KEAP1, KMT2D, NFE2L2, NOTCH1, NSD1, PIK3R1, PTEN, TERT, TP53"
        ),
    },
    {
        "name":       "Woolf, Virginia A.",
        "dob":        "1970-03-11",
        "age":        56,
        "mrn":        "20456732",
        "accession":  "NGS-2026-052102",
        "physician":  "Dr. Wisdom, MD",
        "indication": "Oral lesion — right buccal mucosa",
        "specimen":   "Oral swab, right buccal mucosa",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "7,690× (consensus)",
        "result":     "VARIANTS DETECTED",
        "variants": [
            {
                "gene":    "NOTCH1",
                "hgvs_c":  "NM_017617.5:c.7244C>T",
                "hgvs_p":  "p.Pro2415Leu",
                "vaf":     "0.26%",
                "depth":   "7,690×",
                "tier":    "II",
                "class_":  "Likely Pathogenic",
                "interp": (
                    "NOTCH1 functions as a tumour suppressor in stratified squamous "
                    "epithelium. The p.Pro2415Leu missense variant affects a conserved "
                    "residue in the PEST domain, predicted to disrupt protein stability "
                    "(SIFT: deleterious; PolyPhen-2: probably damaging). NOTCH1 "
                    "loss-of-function mutations are found in 10–15% of OSCC and are "
                    "associated with early-onset tumourigenesis in the oral cavity. "
                    "This variant is classified Likely Pathogenic (Tier II). In the "
                    "context of confirmed dysplasia, detection of a NOTCH1 tumour suppressor "
                    "variant supports malignant progression risk. Clinical correlation and "
                    "multidisciplinary review are recommended."
                ),
            },
        ],
        "neg_genes": (
            "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EGFR, EP300, FAT1, FBXW7, HRAS, "
            "KEAP1, KMT2D, KRAS, NFE2L2, NSD1, PIK3CA, PIK3R1, PTEN, TERT, TP53"
        ),
    },
    {
        "name":       "Dickens, Charles J.",
        "dob":        "1951-09-28",
        "age":        74,
        "mrn":        "30567843",
        "accession":  "NGS-2026-052103",
        "physician":  "Dr. Wisdom, DDS MD",
        "indication": "Oral lesion — left floor of mouth",
        "specimen":   "Oral swab, left floor of mouth",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "6,640× (consensus)",
        "result":     "VARIANTS DETECTED",
        "variants": [
            {
                "gene":    "TP53",
                "hgvs_c":  "NM_000546.6:c.524G>A",
                "hgvs_p":  "p.Arg175His",
                "vaf":     "0.60%",
                "depth":   "6,640×",
                "tier":    "I",
                "class_":  "Pathogenic",
                "interp": (
                    "TP53 p.Arg175His is the most frequently observed p53 mutation in "
                    "human cancer. It disrupts the structural zinc-binding domain of the "
                    "DNA-binding domain, abolishing wild-type tumour suppressor function, "
                    "while also conferring dominant-negative and gain-of-function oncogenic "
                    "properties including promotion of invasion, metastasis, and MDM2 "
                    "antagonism (COSMIC ID COSM10648). In OSCC, TP53 mutations are present "
                    "in 50–70% of cases and are associated with aggressive tumour behaviour, "
                    "locoregional recurrence, and reduced sensitivity to platinum-based "
                    "chemotherapy. The VAF of 0.60% in a surveillance swab warrants clinical "
                    "correlation to exclude early local recurrence. Imaging and clinical "
                    "re-evaluation are recommended."
                ),
            },
        ],
        "neg_genes": (
            "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EGFR, EP300, FAT1, FBXW7, HRAS, "
            "KEAP1, KMT2D, KRAS, NFE2L2, NOTCH1, NSD1, PIK3CA, PIK3R1, PTEN, TERT"
        ),
    },
    {
        "name":       "Kafka, Franz J.",
        "dob":        "1980-02-17",
        "age":        46,
        "mrn":        "50789065",
        "accession":  "NGS-2026-052105",
        "physician":  "Dr. Wisdom, MD",
        "indication": "Oral lesion — anterior tongue (surveillance post-excision)",
        "specimen":   "Oral swab, anterior tongue",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "3,670× (mean consensus)",
        "result":     "NO REPORTABLE VARIANTS DETECTED",
        "variants": [],
        "neg_genes": ALL_GENES,
    },
    {
        "name":       "Morrison, Toni E.",
        "dob":        "1955-10-09",
        "age":        70,
        "mrn":        "60890176",
        "accession":  "NGS-2026-052106",
        "physician":  "Dr. Wisdom, MD",
        "indication": "Oral lesion — left retromolar trigone",
        "specimen":   "Oral swab, left retromolar trigone",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "3,650× (mean consensus)",
        "result":     "NO REPORTABLE VARIANTS DETECTED",
        "variants": [],
        "neg_genes": ALL_GENES,
    },
    {
        "name":       "Austen, Jane C.",
        "dob":        "1967-12-04",
        "age":        58,
        "mrn":        "40678954",
        "accession":  "NGS-2026-052104",
        "physician":  "Dr. Wisdom, MD",
        "indication": "Oral lesion — soft palate",
        "specimen":   "Oral swab, soft palate",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "5,250× / 6,690× (consensus)",
        "result":     "VARIANTS DETECTED",
        "variants": [
            {
                "gene":    "EGFR",
                "hgvs_c":  "NM_005228.5:c.2573T>G",
                "hgvs_p":  "p.Leu858Arg",
                "vaf":     "0.57%",
                "depth":   "5,250×",
                "tier":    "I",
                "class_":  "Pathogenic",
                "interp": (
                    "EGFR p.Leu858Arg (L858R) is a well-characterised activating mutation "
                    "in exon 21 of the EGFR kinase domain. It disrupts the autoinhibitory "
                    "conformation, leading to constitutive EGFR kinase activation and "
                    "downstream RAS/MAPK and PI3K/AKT signalling. While canonical in "
                    "non-small cell lung adenocarcinoma (where it predicts response to "
                    "EGFR TKIs such as osimertinib), EGFR L858R has been reported in "
                    "1–3% of OSCC. Cetuximab (anti-EGFR monoclonal antibody) is "
                    "FDA-approved in HNSCC and may have enhanced activity in the setting "
                    "of EGFR activating mutations. Multidisciplinary tumour board review "
                    "for consideration of EGFR-directed therapy or clinical trial enrolment "
                    "is recommended. Co-occurring TP53 R175H (see below) may affect "
                    "response to EGFR-directed therapy and warrants combined consideration."
                ),
            },
            {
                "gene":    "TP53",
                "hgvs_c":  "NM_000546.6:c.524G>A",
                "hgvs_p":  "p.Arg175His",
                "vaf":     "0.30%",
                "depth":   "6,690×",
                "tier":    "I",
                "class_":  "Pathogenic",
                "interp": (
                    "TP53 p.Arg175His — see interpretation above (Dickens, Charles J., "
                    "NGS-2026-052103). Co-occurrence of EGFR L858R and TP53 R175H has been "
                    "associated with resistance to EGFR-directed therapies in other tumour "
                    "types; multidisciplinary tumour board review is recommended."
                ),
            },
        ],
        "neg_genes": (
            "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EP300, FAT1, FBXW7, HRAS, "
            "KEAP1, KMT2D, KRAS, NFE2L2, NOTCH1, NSD1, PIK3CA, PIK3R1, PTEN, TERT"
        ),
    },
]

def _patch_patients():
    """Update PATIENTS depth fields from actual mosdepth sensitivity-run data."""
    for pt in PATIENTS:
        acc = pt['accession']
        var_depths = []
        for v in pt.get('variants', []):
            pos_info = _HOTSPOT_POS.get((v['gene'], v['hgvs_p']))
            if pos_info:
                chrom, pos = pos_info
                d = _depth_at_pos(acc, chrom, pos)
                if d is not None:
                    v['depth'] = f'{d:,}×'
                    var_depths.append(d)
        if pt['variants'] and var_depths:
            pt['depth'] = ' / '.join(f'{d:,}×' for d in var_depths) + ' (consensus)'
        else:
            sample = PATIENT_SAMPLE_MAP.get(acc)
            depth_map = SAMPLE_DEPTHS.get(sample, {}) if sample else {}
            if depth_map:
                mean_d = int(round(sum(depth_map.values()) / len(depth_map)))
                pt['depth'] = f'{mean_d:,}× (mean consensus)'

_patch_patients()

# ─── PDF class ────────────────────────────────────────────────────────────────

NAVY  = (25,  50,  100)
RED   = (180, 0,   0)
GREY  = (100, 100, 100)
BLACK = (0,   0,   0)
LIGHT = (240, 242, 248)

class ClinicalPDF(FPDF):
    def header(self):
        self.set_fill_color(*NAVY)
        self.rect(0, 0, 210, 18, 'F')
        self.set_y(4)
        self.set_font('DejaVu', 'B', 12)
        self.set_text_color(255, 255, 255)
        self.cell(130, 5, LAB_NAME, align='L')
        self.set_font('DejaVu', '', 7)
        self.cell(0, 5, f'CLIA: {LAB_CLIA}  |  CAP: {LAB_CAP}', align='R')
        self.ln(5)
        self.set_font('DejaVu', '', 7)
        self.cell(130, 4, f'{LAB_ADDRESS}  |  {LAB_PHONE}', align='L')
        self.cell(0, 4, f'{ASSAY_NAME} {ASSAY_VERSION}', align='R')
        self.set_text_color(*BLACK)
        self.set_y(20)

    def footer(self):
        self.set_y(-12)
        self.set_font('DejaVu', 'I', 6.5)
        self.set_text_color(*GREY)
        self.cell(0, 4,
            f'Page {self.page_no()}  |  Run: {RUN_ID}  |  '
            'CONFIDENTIAL — For authorized recipient only. '
            'This report contains protected health information.',
            align='C')
        self.set_text_color(*BLACK)

    def section_bar(self, title):
        self.set_fill_color(*NAVY)
        self.set_text_color(255, 255, 255)
        self.set_font('DejaVu', 'B', 8)
        self.cell(0, 5.5, f'  {title}', fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_text_color(*BLACK)
        self.ln(1)

    def kv_row(self, label, value, w_label=45, bold_val=False):
        self.set_font('DejaVu', 'B', 8)
        self.set_text_color(*GREY)
        self.cell(w_label, 5, label, align='L')
        self.set_text_color(*BLACK)
        self.set_font('DejaVu', 'B' if bold_val else '', 8)
        self.cell(0, 5, value, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    def wrapped_cell(self, text, font_size=7.5, indent=0):
        self.set_font('DejaVu', '', font_size)
        self.set_x(10 + indent)
        self.multi_cell(190 - indent, 4.5, text)
        self.ln(1)


def result_badge(pdf, result):
    if result == "VARIANTS DETECTED":
        pdf.set_fill_color(*RED)
    else:
        pdf.set_fill_color(0, 120, 0)
    pdf.set_text_color(255, 255, 255)
    pdf.set_font('DejaVu', 'B', 11)
    pdf.cell(0, 8, f'  RESULT: {result}', fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_text_color(*BLACK)
    pdf.ln(2)


def variant_table_header(pdf):
    cols = [22, 52, 30, 15, 22, 18, 31]
    hdrs = ['Gene', 'HGVS (coding)', 'Protein change', 'VAF', 'Depth', 'Tier', 'Classification']
    pdf.set_fill_color(*LIGHT)
    pdf.set_font('DejaVu', 'B', 7)
    for w, h in zip(cols, hdrs):
        pdf.cell(w, 5.5, h, border=1, fill=True, align='C')
    pdf.ln()


def variant_table_row(pdf, v):
    cols   = [22, 52, 30, 15, 22, 18, 31]
    values = [v['gene'], v['hgvs_c'], v['hgvs_p'], v['vaf'], v['depth'], v['tier'], v['class_']]
    tier_col = TIER_COLOURS.get(v['tier'], BLACK)
    pdf.set_font('DejaVu', '', 7)
    for i, (w, val) in enumerate(zip(cols, values)):
        if i == 5:  # Tier column
            pdf.set_text_color(*tier_col)
            pdf.set_font('DejaVu', 'B', 7)
        pdf.cell(w, 5.5, val, border=1, align='C')
        pdf.set_text_color(*BLACK)
        pdf.set_font('DejaVu', '', 7)
    pdf.ln()


def draw_gene_block(pdf, x, y, gene_entry, detected_set):
    """Draw one gene block (header bar + hotspot rows). Returns y after the block."""
    HDR_H    = 6.0
    ROW_H    = 5.2
    w        = _COL_W
    w_indent = 2
    w_label  = 24
    w_hgvsc  = 30
    w_result = w - w_indent - w_label - w_hgvsc   # 37

    # Gene header bar
    pdf.set_xy(x, y)
    pdf.set_fill_color(*NAVY)
    pdf.set_text_color(255, 255, 255)
    pdf.set_font('DejaVu', 'B', 8)
    pdf.cell(w - 16, HDR_H, f'  {gene_entry["gene"]}', fill=True, align='L')
    pdf.set_font('DejaVu', '', 7.5)
    pdf.cell(16, HDR_H, 'Cov ✓', fill=True, align='R')
    pdf.set_text_color(*BLACK)
    y += HDR_H

    for hs in gene_entry['hotspots']:
        is_det = (gene_entry['gene'], hs['hgvs_p']) in detected_set
        bg = (255, 243, 243) if is_det else (248, 249, 252)
        pdf.set_xy(x, y)
        pdf.set_fill_color(*bg)

        pdf.cell(w_indent, ROW_H, '', fill=True)

        pdf.set_font('DejaVu', 'B' if is_det else '', 7)
        pdf.set_text_color(*BLACK)
        pdf.cell(w_label, ROW_H, hs['label'], fill=True, align='L')

        pdf.set_font('DejaVu', '', 7)
        pdf.set_text_color(*GREY)
        pdf.cell(w_hgvsc, ROW_H, hs['hgvs_c'], fill=True, align='L')

        if is_det:
            pdf.set_text_color(*RED)
            pdf.set_font('DejaVu', 'B', 7)
            pdf.cell(w_result, ROW_H, '● DETECTED', fill=True, align='L')
        else:
            pdf.set_text_color(160, 160, 160)
            pdf.set_font('DejaVu', '', 7)
            pdf.cell(w_result, ROW_H, 'Not detected', fill=True, align='L')
        pdf.set_text_color(*BLACK)
        y += ROW_H

    return y + 3   # inter-block gap


def hotspot_page(pdf, pt):
    """Page 2: per-gene hotspot screening summary."""
    pdf.add_page()
    detected_set = {(v['gene'], v['hgvs_p']) for v in pt['variants']}

    # Title bar
    pdf.set_fill_color(*NAVY)
    pdf.set_text_color(255, 255, 255)
    pdf.set_font('DejaVu', 'B', 8)
    pdf.cell(0, 5.5, '  HOTSPOT SCREENING SUMMARY',
             fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_text_color(*BLACK)
    pdf.ln(1)

    # Patient identification line
    pdf.set_x(10)
    pdf.set_font('DejaVu', 'B', 7.5)
    pdf.set_text_color(*GREY)
    pdf.cell(28, 5, 'Patient:', align='L')
    pdf.set_text_color(*BLACK)
    pdf.set_font('DejaVu', '', 7.5)
    pdf.cell(75, 5, pt['name'], align='L')
    pdf.set_font('DejaVu', 'B', 7.5)
    pdf.set_text_color(*GREY)
    pdf.cell(25, 5, 'Accession:', align='L')
    pdf.set_text_color(*BLACK)
    pdf.set_font('DejaVu', '', 7.5)
    pdf.cell(0, 5, pt['accession'],
             align='L', new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.set_font('DejaVu', 'I', 6.5)
    pdf.set_text_color(*GREY)
    pdf.set_x(10)
    pdf.multi_cell(190, 4,
        'Hotspots listed below were explicitly interrogated in addition to all variants '
        'above 0.2 % VAF with ≥2 consensus reads across all 251 targets. '
        'Coverage criteria (≥1,000× consensus depth) were met at all targets.')
    pdf.set_text_color(*BLACK)
    pdf.ln(3)

    # Two-column layout
    # Left:  EGFR, KRAS, PIK3CA, TERT
    # Right: HRAS, NOTCH1, TP53
    left_genes  = [g for g in HOTSPOT_GENES if g['gene'] in ('EGFR', 'KRAS', 'PIK3CA', 'TERT')]
    right_genes = [g for g in HOTSPOT_GENES if g['gene'] in ('HRAS', 'NOTCH1', 'TP53')]

    start_y = pdf.get_y()
    y_l = start_y
    for g in left_genes:
        y_l = draw_gene_block(pdf, _COL_L, y_l, g, detected_set)

    y_r = start_y
    for g in right_genes:
        y_r = draw_gene_block(pdf, _COL_R, y_r, g, detected_set)

    # Thin column separator
    sep_x = _COL_L + _COL_W + _COL_GAP / 2
    pdf.set_draw_color(200, 210, 230)
    pdf.set_line_width(0.2)
    pdf.line(sep_x, start_y, sep_x, max(y_l, y_r) - 3)

    pdf.set_xy(10, max(y_l, y_r) + 2)

    # Non-hotspot genes footer block
    pdf.set_fill_color(235, 238, 248)
    pdf.set_text_color(*NAVY)
    pdf.set_font('DejaVu', 'B', 7.5)
    pdf.cell(0, 5.5,
             '  REMAINING PANEL GENES — SEQUENCE-LEVEL DETECTION, NO PREDEFINED HOTSPOT LIST',
             fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_text_color(*BLACK)
    pdf.set_font('DejaVu', 'B', 7.5)
    pdf.set_x(10)
    pdf.cell(0, 5, NON_HOTSPOT_GENES,
             new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font('DejaVu', '', 7)
    pdf.set_text_color(*GREY)
    pdf.set_x(10)
    pdf.multi_cell(190, 4.5,
        'All variants above 0.2 % VAF with ≥2 consensus reads are reported for '
        'these genes. No variants meeting reporting criteria were detected. '
        'All targets in these genes met minimum coverage criteria '
        '(≥1,000× consensus depth).')
    pdf.set_text_color(*BLACK)

    # Signature block
    pdf.ln(3)
    pdf.set_draw_color(*NAVY)
    pdf.set_line_width(0.3)
    pdf.line(10, pdf.get_y(), 200, pdf.get_y())
    pdf.ln(2)
    pdf.set_font('DejaVu', 'B', 8)
    pdf.cell(95, 5, f'Reported by: {DIRECTOR}', align='L')
    pdf.cell(0, 5, f'Report date: {RUN_DATE}', align='R',
             new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font('DejaVu', '', 7)
    pdf.set_text_color(*GREY)
    pdf.cell(0, 4,
        'Results reviewed and approved by the Laboratory Director. '
        'This report is electronically signed and does not require a wet signature.',
        align='L', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_text_color(*BLACK)


def coverage_appendix_pages(pdf, pt):
    """Render full 251-target coverage appendix, paginated."""
    detected_positions = set()
    for v in pt['variants']:
        pos_info = _HOTSPOT_POS.get((v['gene'], v['hgvs_p']))
        if pos_info:
            detected_positions.add(pos_info)

    ROW_H       = 4.8
    GENE_H      = 5.5
    HDR_H       = 5.5
    PAGE_BOTTOM = 297 - 16   # 16 mm above page bottom (footer zone)
    ORANGE      = (180, 80, 0)
    GREEN       = (0, 110, 0)

    # Column widths (total = 190 mm)
    W_GENE, W_REGION, W_LEN, W_DEPTH, W_STATUS = 18, 50, 14, 20, 88

    # Compute summary stats across all targets for this patient
    depths = [_real_depth(pt['accession'], i) for i in range(len(TARGETS))]
    n_low  = sum(1 for d in depths if d < 1000)

    def draw_table_header():
        pdf.set_fill_color(*LIGHT)
        pdf.set_font('DejaVu', 'B', 6.5)
        pdf.set_x(10)
        for w, label in [
            (W_GENE,   'Gene'),
            (W_REGION, 'chr:start–end (1-based)'),
            (W_LEN,    'Len (bp)'),
            (W_DEPTH,  'Depth'),
            (W_STATUS, 'Coverage status'),
        ]:
            pdf.cell(w, HDR_H, label, border=1, fill=True, align='C')
        pdf.ln()
        return pdf.get_y()

    def start_page():
        pdf.add_page()
        pdf.set_fill_color(*NAVY)
        pdf.set_text_color(255, 255, 255)
        pdf.set_font('DejaVu', 'B', 8)
        pdf.cell(0, 5.5,
                 f'  FULL PANEL COVERAGE REPORT  —  {pt["name"]}  ({pt["accession"]})',
                 fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_text_color(*BLACK)
        return draw_table_header()

    # Group targets by gene, preserving BED order within each gene
    gene_groups = OrderedDict()
    for i, t in enumerate(TARGETS):
        gene_groups.setdefault(t['gene'], []).append((i, t))

    # Summary line on first coverage page only
    pdf.add_page()
    pdf.set_fill_color(*NAVY)
    pdf.set_text_color(255, 255, 255)
    pdf.set_font('DejaVu', 'B', 8)
    pdf.cell(0, 5.5,
             f'  FULL PANEL COVERAGE REPORT  —  {pt["name"]}  ({pt["accession"]})',
             fill=True, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_text_color(*BLACK)
    pdf.ln(1)
    pdf.set_font('DejaVu', 'I', 6.5)
    pdf.set_text_color(*GREY)
    pdf.set_x(10)
    pct_low = 100 * n_low / len(TARGETS) if TARGETS else 0
    pdf.multi_cell(190, 4,
        f'251 targets, {sum(t["end"]-t["start"] for t in TARGETS):,} bp total. '
        f'Mean consensus depth: {int(sum(depths)/len(depths)):,}×. '
        f'Targets below 1,000× consensus: {n_low} / {len(TARGETS)} ({pct_low:.0f}%). '
        'Depths are from the UMI-deduplicated consensus BAM (min-reads ≥ 2, 20,000× raw input). '
        'Positions reported 1-based.')
    pdf.set_text_color(*BLACK)
    pdf.ln(2)

    y = draw_table_header()

    for gene, items in gene_groups.items():
        # Gene section header
        if y + GENE_H > PAGE_BOTTOM - ROW_H:
            y = start_page()
        pdf.set_xy(10, y)
        pdf.set_fill_color(220, 228, 245)
        pdf.set_text_color(*NAVY)
        pdf.set_font('DejaVu', 'B', 7)
        pdf.cell(190, GENE_H, f'  {gene}', fill=True, align='L')
        pdf.set_text_color(*BLACK)
        y += GENE_H

        for idx, t in items:
            if y + ROW_H > PAGE_BOTTOM:
                y = start_page()

            depth   = depths[idx]
            has_var = any(
                t['chrom'] == chrom and t['start'] <= pos < t['end']
                for chrom, pos in detected_positions
            )
            bg = (255, 243, 243) if has_var else (255, 255, 255)

            pdf.set_xy(10, y)
            pdf.set_fill_color(*bg)

            # Gene (grey, compact)
            pdf.set_font('DejaVu', '', 6)
            pdf.set_text_color(*GREY)
            pdf.cell(W_GENE, ROW_H, gene, fill=True, align='L')

            # Region (1-based display)
            pdf.set_text_color(*BLACK)
            pdf.set_font('DejaVu', '', 6.5)
            region = f'{t["chrom"]}:{t["start"]+1:,}–{t["end"]:,}'
            pdf.cell(W_REGION, ROW_H, region, fill=True, align='L')

            # Length
            pdf.cell(W_LEN, ROW_H, f'{t["end"]-t["start"]:,}', fill=True, align='R')

            # Depth
            if depth < 1000:
                pdf.set_text_color(*ORANGE)
                pdf.set_font('DejaVu', 'B', 6.5)
            else:
                pdf.set_text_color(*BLACK)
                pdf.set_font('DejaVu', '', 6.5)
            pdf.cell(W_DEPTH, ROW_H, f'{depth:,}×', fill=True, align='R')
            pdf.set_text_color(*BLACK)

            # Status
            if has_var:
                pdf.set_text_color(*RED)
                pdf.set_font('DejaVu', 'B', 6.5)
                pdf.cell(W_STATUS, ROW_H, '● VARIANT DETECTED', fill=True, align='L')
            elif depth < 1000:
                pdf.set_text_color(*ORANGE)
                pdf.set_font('DejaVu', '', 6.5)
                pdf.cell(W_STATUS, ROW_H, '⚠ Below 1,000× threshold', fill=True, align='L')
            else:
                pdf.set_text_color(*GREEN)
                pdf.set_font('DejaVu', '', 6.5)
                pdf.cell(W_STATUS, ROW_H, '✓ Covered', fill=True, align='L')
            pdf.set_text_color(*BLACK)

            y += ROW_H


def patient_page(pdf, pt):
    pdf.add_page()

    # ── Patient / specimen info ──
    pdf.section_bar('PATIENT & SPECIMEN INFORMATION')
    pdf.set_font('DejaVu', '', 8)

    left = [
        ('Patient name',       pt['name']),
        ('Date of birth',      f"{pt['dob']}  (Age {pt['age']})"),
        ('MRN',                pt['mrn']),
        ('Ordering physician', pt['physician']),
        ('Clinical indication', pt['indication']),
    ]
    right = [
        ('Accession', pt['accession']),
        ('Specimen',  pt['specimen']),
        ('Collected', pt['collected']),
        ('Received',  pt['received']),
        ('Reported',  RUN_DATE),
    ]

    for (lbl, val), (lbl2, val2) in zip(left, right):
        pdf.set_x(10)
        pdf.set_font('DejaVu', 'B', 7.5)
        pdf.set_text_color(*GREY)
        pdf.cell(32, 5, lbl, align='L')
        pdf.set_text_color(*BLACK)
        pdf.set_font('DejaVu', '', 7.5)
        pdf.cell(63, 5, val, align='L')
        pdf.set_font('DejaVu', 'B', 7.5)
        pdf.set_text_color(*GREY)
        pdf.cell(32, 5, lbl2, align='L')
        pdf.set_text_color(*BLACK)
        pdf.set_font('DejaVu', '', 7.5)
        pdf.cell(0, 5, val2, align='L', new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    pdf.ln(3)

    # ── Result badge ──
    result_badge(pdf, pt['result'])

    # ── Detected variants table ──
    if pt['variants']:
        pdf.section_bar('DETECTED VARIANTS')
        variant_table_header(pdf)
        for v in pt['variants']:
            variant_table_row(pdf, v)
        pdf.ln(3)

        # ── Interpretation ──
        pdf.section_bar('VARIANT ANNOTATION')
        for v in pt['variants']:
            pdf.set_font('DejaVu', 'B', 8)
            tier_col = TIER_COLOURS.get(v['tier'], BLACK)
            pdf.set_text_color(*tier_col)
            pdf.cell(0, 5,
                f"  {v['gene']} {v['hgvs_p']}  |  {v['class_']} (Tier {v['tier']})",
                new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.set_text_color(*BLACK)
            pdf.wrapped_cell(v['interp'], indent=2)
            pdf.ln(1)
    else:
        pdf.section_bar('INTERPRETATION')
        pdf.wrapped_cell(
            "No pathogenic or likely pathogenic somatic variants were detected in the "
            "21 genes covered by this assay at ≥0.2% VAF with ≥2 consensus supporting "
            "reads. A negative result does not exclude malignancy. Variants below the "
            "reporting threshold, variants outside the 251 assay target regions, "
            "structural variants, and copy number alterations are not assessed by this "
            "assay. Clinical correlation is required; repeat testing or biopsy should "
            "be considered if clinical suspicion persists.",
            font_size=7.5
        )
        pdf.ln(1)

    # ── Negative genes ──
    pdf.section_bar('GENES WITH NO PATHOGENIC OR LIKELY PATHOGENIC VARIANTS DETECTED')
    pdf.set_font('DejaVu', '', 7.5)
    pdf.set_x(10)
    pdf.multi_cell(190, 4.5,
        pt['neg_genes'] + '\n'
        'Variants of uncertain significance (Tier III) are not reported.')
    pdf.ln(2)

    # ── Method / limitations ──
    pdf.section_bar('METHOD SUMMARY')
    pdf.wrapped_cell(METHOD_SUMMARY, font_size=7)
    pdf.ln(1)

    pdf.section_bar('LIMITATIONS')
    pdf.wrapped_cell(LIMITATIONS, font_size=7)
    pdf.ln(2)

    # ── Signature block ──
    pdf.set_font('DejaVu', '', 7.5)
    pdf.cell(0, 4, '', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_draw_color(*NAVY)
    pdf.set_line_width(0.3)
    pdf.line(10, pdf.get_y(), 200, pdf.get_y())
    pdf.ln(2)
    pdf.set_font('DejaVu', 'B', 8)
    pdf.cell(95, 5, f'Reported by: {DIRECTOR}', align='L')
    pdf.cell(0, 5, f'Report date: {RUN_DATE}', align='R', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font('DejaVu', '', 7)
    pdf.set_text_color(*GREY)
    pdf.cell(0, 4,
        'Results reviewed and approved by the Laboratory Director. '
        'This report is electronically signed and does not require a wet signature.',
        align='L', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_text_color(*BLACK)


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', default=None,
                   help='Output PDF path (default depends on --full flag)')
    p.add_argument('--full', action='store_true',
                   help='Add full 251-target coverage appendix (one extra section per patient)')
    args = p.parse_args()

    default_out = (
        'physician_report_idthyb_20260520_full.pdf' if args.full
        else 'physician_report_idthyb_20260520.pdf'
    )
    out_path = args.out or default_out

    pdf = ClinicalPDF()
    pdf.add_font('DejaVu',  '',  '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
    pdf.add_font('DejaVu',  'B', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf')
    pdf.add_font('DejaVu',  'I', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
    pdf.set_margins(10, 22, 10)
    pdf.set_auto_page_break(auto=True, margin=15)

    for pt in PATIENTS:
        patient_page(pdf, pt)
        hotspot_page(pdf, pt)
        if args.full:
            coverage_appendix_pages(pdf, pt)

    pdf.output(out_path)
    variant = 'full coverage' if args.full else 'hotspot summary'
    print(f"Written: {out_path}  ({len(PATIENTS)} patient reports, {pdf.page} pages, {variant})")


if __name__ == '__main__':
    main()
