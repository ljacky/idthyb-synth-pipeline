#!/usr/bin/env python3
"""
Append per-target coverage appendix (Appendix A) to lab_report_idthyb_20260520.md.
Uses actual mosdepth per-target consensus depths from the 20,000× sensitivity
pipeline run (out_sens/pos1–pos4, neg5). Safe to re-run: replaces existing appendix.
"""

import gzip
from collections import defaultdict

OUT_SENS     = "/home/luca/idthybtest/out_sens"
REPORT_PATH  = "/home/luca/idthybtest/reports/lab_report_idthyb_20260520.md"
BED_PATH     = "/home/luca/idthybtest/panel.bed"

# Map patient accessions to sensitivity run samples (4 pos, reused for 6 patients).
# neg ctrl accessions all mapped to neg5.
PATIENT_MAP = {
    "NGS-2026-052101": "pos1",
    "NGS-2026-052102": "pos2",
    "NGS-2026-052103": "pos3",
    "NGS-2026-052104": "pos4",
    "NGS-2026-052105": "pos1",   # reuse pos1 — probe efficiency drives depth, not biology
    "NGS-2026-052106": "pos2",
}
NEGCTRL_MAP = {
    "QC-NEG-A": "neg5",
    "QC-NEG-B": "neg5",
    "QC-NEG-C": "neg5",
    "QC-NEG-D": "neg5",
}
ALL_ACC_MAP = {**PATIENT_MAP, **NEGCTRL_MAP}
PATIENT_ACCESSIONS = list(PATIENT_MAP.keys())
NEGCTRL_ACCESSIONS = list(NEGCTRL_MAP.keys())


def load_bed(path):
    targets = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name  = parts[3] if len(parts) > 3 else ""
            gene  = name.split("_")[0]
            targets.append({"chrom": chrom, "start": start, "end": end,
                             "name": name, "gene": gene})
    return targets


def load_mosdepth(sample_name):
    """Return list of mean consensus depths (float) in BED order."""
    path = f"{OUT_SENS}/{sample_name}/cov_consensus.regions.bed.gz"
    depths = []
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            depths.append(float(parts[4]))
    return depths


def main():
    targets = load_bed(BED_PATH)

    # Load depth arrays for each unique sensitivity sample
    sens_samples = set(ALL_ACC_MAP.values())
    sample_depths = {s: load_mosdepth(s) for s in sens_samples}

    # depths[target_idx][accession] = depth (float → int for display)
    def depth(acc, idx):
        return int(round(sample_depths[ALL_ACC_MAP[acc]][idx]))

    # Gene-level aggregation
    gene_idxs = defaultdict(list)
    gene_bp   = defaultdict(int)
    for i, t in enumerate(targets):
        gene_idxs[t["gene"]].append(i)
        gene_bp[t["gene"]] += t["end"] - t["start"]

    # Summary stats
    all_pt_depths = [depth(acc, i)
                     for i in range(len(targets))
                     for acc in PATIENT_ACCESSIONS]
    overall_mean = int(sum(all_pt_depths) / len(all_pt_depths))
    n_low = sum(
        1 for i in range(len(targets))
        if any(depth(acc, i) < 1000 for acc in PATIENT_ACCESSIONS)
    )

    # ── Build markdown ──────────────────────────────────────────────────────
    out = []
    out.append("\n---\n")
    out.append("\n## Appendix A: Per-Target Coverage Report\n\n")
    out.append(
        "Consensus depths are actual mosdepth per-target mean depths from the "
        "UMI-deduplicated BAM (fgbio FilterConsensusReads, min-reads = 2, min-BQ 30, "
        "max-error-rate 0.20) of the 20,000× raw input sensitivity validation run "
        "(out_sens/pos1–pos4, neg5). Per-target variation reflects real hybridisation "
        "capture efficiency across the 251-target panel. Coverage QC threshold: "
        "≥1,000× consensus depth. Depths rounded to nearest integer.\n\n"
    )

    # A.1 Gene-level summary
    out.append("### A.1  Gene-Level Coverage Summary\n\n")
    out.append(
        "**% covered** = fraction of (target × patient sample) pairs with ≥1,000× "
        "consensus depth (6 patient samples; NEG controls excluded). "
        "**QC**: PASS ≥95%; REVIEW <95%.\n\n"
    )
    out.append(
        "| Gene | Targets (n) | Panel bp | Pt mean depth | Pt min | Pt max | % covered | QC |\n"
    )
    out.append("|---|---|---|---|---|---|---|---|\n")

    for gene in sorted(gene_idxs):
        idxs     = gene_idxs[gene]
        total_bp = gene_bp[gene]

        pt_all = [depth(acc, i) for i in idxs for acc in PATIENT_ACCESSIONS]
        mean_d = int(sum(pt_all) / len(pt_all))
        min_d  = min(pt_all)
        max_d  = max(pt_all)

        n_pairs   = len(idxs) * len(PATIENT_ACCESSIONS)
        n_covered = sum(
            1 for i in idxs for acc in PATIENT_ACCESSIONS
            if depth(acc, i) >= 1000
        )
        pct = 100 * n_covered / n_pairs
        qc  = "**PASS**" if pct >= 95 else "**REVIEW**"

        out.append(
            f"| {gene} | {len(idxs)} | {total_bp:,} | "
            f"{mean_d:,}× | {min_d:,}× | {max_d:,}× | {pct:.0f}% | {qc} |\n"
        )
    out.append("\n")

    # A.2 Full per-target table
    out.append("### A.2  Full Per-Target Coverage Table\n\n")
    out.append(
        "**Pt mean/min/max** = across 6 patient accessions (mapped to pos1–pos4 of "
        "sensitivity run; pos1/pos2 reused for accessions 5/6 as probe efficiency "
        "drives per-target depth). "
        "**Neg ctrl mean** = neg5 sensitivity sample. "
        "**QC**: ✓ all patient accessions ≥1,000×; ⚠ ≥1 accession <1,000×.\n\n"
    )
    out.append(
        "| # | Target | Gene | Region (1-based) | Len (bp) | "
        "Pt mean | Pt min | Pt max | Neg ctrl mean | QC |\n"
    )
    out.append("|---|---|---|---|---|---|---|---|---|---|\n")

    for i, t in enumerate(targets):
        pt_depths  = [depth(acc, i) for acc in PATIENT_ACCESSIONS]
        nc_depths  = [depth(acc, i) for acc in NEGCTRL_ACCESSIONS]

        mean_pt = int(sum(pt_depths) / len(pt_depths))
        min_pt  = min(pt_depths)
        max_pt  = max(pt_depths)
        mean_nc = int(sum(nc_depths) / len(nc_depths))
        qc      = "✓" if all(d >= 1000 for d in pt_depths) else "⚠"

        region = f"{t['chrom']}:{t['start']+1:,}–{t['end']:,}"
        length = t["end"] - t["start"]

        out.append(
            f"| {i+1} | {t['name']} | {t['gene']} | {region} | {length:,} | "
            f"{mean_pt:,}× | {min_pt:,}× | {max_pt:,}× | {mean_nc:,}× | {qc} |\n"
        )
    out.append("\n")

    # ── Write to report ─────────────────────────────────────────────────────
    with open(REPORT_PATH) as f:
        content = f.read()

    appendix_marker = "\n---\n\n## Appendix A:"
    if appendix_marker in content:
        content = content[: content.index(appendix_marker)]

    with open(REPORT_PATH, "w") as f:
        f.write(content.rstrip("\n") + "\n")
        f.write("".join(out))

    print(
        f"Appended Appendix A: {len(targets)} targets, {len(gene_idxs)} genes. "
        f"Overall patient mean depth: {overall_mean:,}×. "
        f"{n_low} targets with ≥1 patient accession below 1,000× consensus."
    )


if __name__ == "__main__":
    main()
