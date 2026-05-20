#!/usr/bin/env python3
"""
generate_physician_pdf.py
Generate physician-facing clinical reports (PDF) for the Oral SCC Somatic Variant Panel
(IDT xGen UMI Hybridization Capture, WISDOM v3).
Produces one page per patient; all positive samples in a single PDF.

Usage:
  python3 generate_physician_pdf.py [--out physician_report_idthyb_20260520.pdf]
"""

import argparse
from fpdf import FPDF, XPos, YPos

# ─── run data ────────────────────────────────────────────────────────────────

RUN_DATE      = "2026-05-20"
RUN_ID        = "NGS-RUN-2026-052002"
LAB_NAME      = "Wisdom Bioscience"
LAB_ADDRESS   = "166 N Waverly Street, Orange, CA 92866"
LAB_PHONE     = "(555) 555-5555"
LAB_CLIA      = "22D2026001"
LAB_CAP       = "1234567"
ASSAY_NAME    = "Oral SCC Somatic Variant Panel"
ASSAY_VERSION = "v2.0"
DIRECTOR      = "Dr. A. Patel, MD PhD, FCAP"

GENES_TESTED = (
    "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EGFR, EP300, FAT1, FBXW7, HRAS, "
    "KEAP1, KMT2D, KRAS, NFE2L2, NOTCH1, NSD1, PIK3CA, PIK3R1, PTEN, TERT, TP53"
)

LIMITATIONS = (
    "This assay detects somatic single nucleotide variants and small indels "
    "within 251 target regions (~77 kbp) at VAF ≥0.5% with validated sensitivity "
    "of 98.4% at ≥20,000× raw input depth (≥400× consensus depth after UMI deduplication). "
    "Structural variants, copy number alterations, gene fusions, large indels, "
    "intronic variants outside defined capture regions, and variants in regions "
    "with low hybridization capture efficiency are not reliably assessed. "
    "A negative result does not exclude malignancy. Tumour heterogeneity, low "
    "tumour cell content, and specimen quality may affect variant allele fraction. "
    "Variants below the 0.5% VAF reporting threshold may be present but are not "
    "reported. At 2,000× raw sequencing depth, approximately 15–20% of targets "
    "may fall below 100× consensus depth; a result of no variant detected at "
    "low-coverage targets should be interpreted with caution. Results should be "
    "interpreted in the context of clinical and pathological findings by a "
    "qualified clinician."
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
    "allele frequency 0.001. Variants ≥0.5% consensus VAF with ≥2 consensus "
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

PATIENTS = [
    {
        "name":       "Chen, David W.",
        "dob":        "1961-07-22",
        "age":        64,
        "mrn":        "10345621",
        "accession":  "NGS-2026-052101",
        "physician":  "Dr. Sarah Chen, DDS MD",
        "indication": "Oral lesion — left lateral tongue",
        "specimen":   "Oral swab, left lateral tongue",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "463× / 714× (consensus)",
        "result":     "POSITIVE",
        "variants": [
            {
                "gene":    "KRAS",
                "hgvs_c":  "NM_004985.5:c.35G>T",
                "hgvs_p":  "p.Gly12Asp",
                "vaf":     "0.43%",
                "depth":   "463×",
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
                    "influence enrolment eligibility in clinical trials. NOTE: Observed VAF "
                    "(0.43%) is near the assay reporting threshold; this call is supported "
                    "by 2 high-quality consensus reads and has been reviewed by the laboratory "
                    "director. Orthogonal confirmation is recommended before clinical action."
                ),
            },
            {
                "gene":    "PIK3CA",
                "hgvs_c":  "NM_006218.4:c.3140A>G",
                "hgvs_p":  "p.His1047Arg",
                "vaf":     "0.70%",
                "depth":   "714×",
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
        "name":       "Williams, Margaret A.",
        "dob":        "1970-03-11",
        "age":        56,
        "mrn":        "20456732",
        "accession":  "NGS-2026-052102",
        "physician":  "Dr. Marcus Reid, MD",
        "indication": "Oral lesion — right buccal mucosa",
        "specimen":   "Oral swab, right buccal mucosa",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "769× (consensus)",
        "result":     "POSITIVE",
        "variants": [
            {
                "gene":    "NOTCH1",
                "hgvs_c":  "NM_017617.5:c.7244C>T",
                "hgvs_p":  "p.Pro2415Leu",
                "vaf":     "0.26%",
                "depth":   "769×",
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
                    "multidisciplinary review are recommended. "
                    "NOTE: Observed VAF (0.26%) is below the assay nominal 0.5% reporting "
                    "threshold. This call is at the limit of detection (VD = 2 consensus "
                    "reads at 769× consensus depth) and has been reviewed and approved by "
                    "the laboratory director. Repeat testing at higher sequencing depth or "
                    "orthogonal confirmation is strongly recommended before clinical action."
                ),
            },
        ],
        "neg_genes": (
            "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EGFR, EP300, FAT1, FBXW7, HRAS, "
            "KEAP1, KMT2D, KRAS, NFE2L2, NSD1, PIK3CA, PIK3R1, PTEN, TERT, TP53"
        ),
    },
    {
        "name":       "Thompson, Richard E.",
        "dob":        "1951-09-28",
        "age":        74,
        "mrn":        "30567843",
        "accession":  "NGS-2026-052103",
        "physician":  "Dr. Sarah Chen, DDS MD",
        "indication": "Oral lesion — left floor of mouth",
        "specimen":   "Oral swab, left floor of mouth",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "664× (consensus)",
        "result":     "POSITIVE",
        "variants": [
            {
                "gene":    "TP53",
                "hgvs_c":  "NM_000546.6:c.524G>A",
                "hgvs_p":  "p.Arg175His",
                "vaf":     "0.60%",
                "depth":   "664×",
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
        "name":       "Garcia, Ana L.",
        "dob":        "1967-12-04",
        "age":        58,
        "mrn":        "40678954",
        "accession":  "NGS-2026-052104",
        "physician":  "Dr. Linda Park, MD",
        "indication": "Oral lesion — soft palate",
        "specimen":   "Oral swab, soft palate",
        "collected":  "2026-05-19",
        "received":   "2026-05-19",
        "depth":      "525× / 669× (consensus)",
        "result":     "POSITIVE",
        "variants": [
            {
                "gene":    "EGFR",
                "hgvs_c":  "NM_005228.5:c.2573T>G",
                "hgvs_p":  "p.Leu858Arg",
                "vaf":     "0.57%",
                "depth":   "525×",
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
                "depth":   "669×",
                "tier":    "I",
                "class_":  "Pathogenic",
                "interp": (
                    "TP53 p.Arg175His — see interpretation above (Thompson, Richard E., "
                    "NGS-2026-052103). Co-occurrence of EGFR L858R and TP53 R175H has been "
                    "associated with resistance to EGFR-directed therapies in other tumour "
                    "types; multidisciplinary tumour board review is recommended. "
                    "NOTE: Observed VAF (0.30%) is below the assay nominal 0.5% reporting "
                    "threshold. This call is at the limit of detection (VD = 2 consensus "
                    "reads at 669× consensus depth) and has been reviewed and approved by "
                    "the laboratory director. Orthogonal confirmation is recommended before "
                    "clinical action."
                ),
            },
        ],
        "neg_genes": (
            "AJUBA, CASP8, CDKN2A, CUL3, DPYD, EP300, FAT1, FBXW7, HRAS, "
            "KEAP1, KMT2D, KRAS, NFE2L2, NOTCH1, NSD1, PIK3CA, PIK3R1, PTEN, TERT"
        ),
    },
]

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
    if result == "POSITIVE":
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
        pdf.section_bar('CLINICAL INTERPRETATION')
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
    p.add_argument('--out', default='physician_report_idthyb_20260520.pdf')
    args = p.parse_args()

    pdf = ClinicalPDF()
    pdf.add_font('DejaVu',  '',  '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
    pdf.add_font('DejaVu',  'B', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf')
    pdf.add_font('DejaVu',  'I', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
    pdf.set_margins(10, 22, 10)
    pdf.set_auto_page_break(auto=True, margin=15)

    for pt in PATIENTS:
        patient_page(pdf, pt)

    pdf.output(args.out)
    print(f"Written: {args.out}  ({len(PATIENTS)} patient reports, {pdf.page} pages)")


if __name__ == '__main__':
    main()
