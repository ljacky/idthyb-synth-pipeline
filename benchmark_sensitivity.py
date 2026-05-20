#!/usr/bin/env python3
"""
benchmark_sensitivity.py — report 0.5% VAF sensitivity panel results.

Reads manifest from fastq_sens/manifest.json for ground truth.
Queries consensus VCFs in out_sens/{sample}/{sample}.vardict.consensus.vcf.gz.
Reports per-gene detection summary and AF distribution.

Usage:
  cd ~/idthybtest
  python3 benchmark_sensitivity.py
"""

import argparse
import gzip
import json
import os
import statistics
import subprocess
from collections import defaultdict


def bcftools_query(vcf, chrom, pos):
    region = f"{chrom}:{pos}-{pos}"
    try:
        out = subprocess.check_output(
            ["bcftools", "view", "-H", "-r", region, vcf],
            stderr=subprocess.DEVNULL,
        ).decode()
    except subprocess.CalledProcessError:
        return []
    rows = []
    for line in out.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True)
                    for kv in fields[7].split(";"))
        rows.append({
            "ref":    fields[3],
            "alt":    fields[4],
            "filter": fields[6],
            "af":     float(info.get("AF", 0)),
            "dp":     int(info.get("DP", 0)),
            "vd":     int(info.get("VD", 0)),
            "mq":     info.get("MQ", "?"),
        })
    return rows


def count_vcf_rows(vcf, pass_only=False):
    cmd = ["bcftools", "view", "-H"]
    if pass_only:
        cmd += ["-f", "PASS"]
    cmd.append(vcf)
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
        return len([l for l in out.decode().splitlines() if l and not l.startswith("#")])
    except subprocess.CalledProcessError:
        return None


def mean_depth(cov_gz):
    if not os.path.exists(cov_gz):
        return None
    total, n = 0.0, 0
    with gzip.open(cov_gz, "rt") as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    total += float(parts[4])
                    n += 1
                except ValueError:
                    pass
    return total / n if n else None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fastq-dir", default="fastq_sens")
    p.add_argument("--out-dir",   default="out_sens")
    args = p.parse_args()

    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    manifest_path = os.path.join(args.fastq_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        raise SystemExit(f"Manifest not found: {manifest_path}")

    with open(manifest_path) as f:
        manifest = json.load(f)

    variants = manifest["variants"]   # {vname: {chrom, pos_1based, ref, alt, gene, sample}}

    pos_samples = ["pos1", "pos2", "pos3", "pos4"]
    neg_sample  = "neg5"

    sep = "-" * 110

    # ── Coverage ──────────────────────────────────────────────────────────────
    print("\n" + sep)
    print("COVERAGE SUMMARY (mean over panel.bed, raw and consensus)")
    print(sep)
    print(f"{'Sample':<8}  {'Raw depth':>12}  {'Consensus depth':>16}  {'Fold drop':>10}  {'Label'}")
    for sname in pos_samples + [neg_sample]:
        raw_cov = os.path.join(args.out_dir, sname, "cov_raw.regions.bed.gz")
        con_cov = os.path.join(args.out_dir, sname, "cov_consensus.regions.bed.gz")
        rd = mean_depth(raw_cov)
        cd = mean_depth(con_cov)
        rd_str = f"{rd:.0f}x" if rd else "pending"
        cd_str = f"{cd:.0f}x" if cd else "pending"
        fd_str = f"{rd/cd:.1f}x" if rd and cd else ""
        vcf = os.path.join(args.out_dir, sname, f"{sname}.vardict.consensus.vcf.gz")
        label = "(pending)" if not os.path.exists(vcf) else ""
        print(f"{sname:<8}  {rd_str:>12}  {cd_str:>16}  {fd_str:>10}  {label}")

    # ── Per-gene detection summary ─────────────────────────────────────────────
    print("\n" + sep)
    print("DETECTION SUMMARY BY GENE  (consensus VCF; MISS = not called)")
    print(sep)

    gene_stats = defaultdict(lambda: {
        "total": 0, "pass": 0, "fail_filter": defaultdict(int), "miss": 0, "pending": 0,
    })
    all_results = {}

    # Group variants by sample
    sample_variants = defaultdict(dict)
    for vname, vinfo in variants.items():
        s = vinfo.get("sample", "pos4")
        if s in pos_samples:
            sample_variants[s][vname] = vinfo

    for sname in pos_samples:
        vcf = os.path.join(args.out_dir, sname, f"{sname}.vardict.consensus.vcf.gz")
        vcf_exists = os.path.exists(vcf)

        for vname, vinfo in sample_variants[sname].items():
            gene  = vinfo["gene"]
            chrom = vinfo["chrom"]
            pos   = vinfo["pos_1based"]
            gs    = gene_stats[gene]
            gs["total"] += 1

            if not vcf_exists:
                gs["pending"] += 1
                all_results[vname] = {"status": "pending", "gene": gene, "calls": []}
                continue

            calls = bcftools_query(vcf, chrom, pos)
            if not calls:
                gs["miss"] += 1
                all_results[vname] = {"status": "MISS", "gene": gene, "calls": []}
            else:
                flt = calls[0]["filter"]
                if flt == "PASS":
                    gs["pass"] += 1
                else:
                    gs["fail_filter"][flt] += 1
                all_results[vname] = {"status": flt, "gene": gene, "calls": calls}

    print(f"{'Gene':<18}  {'n':>4}  {'PASS':>6}  {'MISS':>6}  {'det%':>7}  {'Filters'}")
    print("-" * 80)
    total_all = pass_all = miss_all = 0
    for gene in sorted(gene_stats):
        gs   = gene_stats[gene]
        n    = gs["total"]
        pas  = gs["pass"]
        miss = gs["miss"]
        pend = gs["pending"]
        pct  = 100 * pas / n if n > 0 else 0
        filt_str = "  ".join(f"{k}={v}" for k, v in sorted(gs["fail_filter"].items()))
        if pend:
            filt_str += f"  pending={pend}"
        total_all += n; pass_all += pas; miss_all += miss
        print(f"{gene:<18}  {n:>4}  {pas:>6}  {miss:>6}  {pct:>6.1f}%  {filt_str}")
    print("-" * 80)
    pct_all = 100 * pass_all / total_all if total_all else 0
    print(f"{'TOTAL':<18}  {total_all:>4}  {pass_all:>6}  {miss_all:>6}  {pct_all:>6.1f}%")

    # ── False positives ────────────────────────────────────────────────────────
    print("\n" + sep)
    print("FALSE POSITIVES — NEGATIVE CONTROL (neg5, consensus VCF)")
    print(sep)
    vcf = os.path.join(args.out_dir, neg_sample, f"{neg_sample}.vardict.consensus.vcf.gz")
    if os.path.exists(vcf):
        n_total = count_vcf_rows(vcf)
        n_pass  = count_vcf_rows(vcf, pass_only=True)
        print(f"  Total calls: {n_total}   PASS: {n_pass}")
    else:
        print("  pending")

    # ── Non-PASS / MISS details ────────────────────────────────────────────────
    non_pass = [(vn, r) for vn, r in all_results.items()
                if r["status"] not in ("PASS", "pending")]
    if non_pass:
        print("\n" + sep)
        print(f"NON-PASS / MISS DETAILS  ({len(non_pass)} variants)")
        print(sep)
        print(f"{'Variant':<40}  {'Gene':<14}  {'Status':<12}  {'AF':>8}  {'VD/DP':>10}  {'MQ':>4}")
        print("-" * 100)
        for vn, r in sorted(non_pass, key=lambda x: (x[1]["gene"], x[0])):
            if r["calls"]:
                c = r["calls"][0]
                af_s = f"{c['af']:.4f}"; vd_s = f"{c['vd']}/{c['dp']}"; mq_s = str(c["mq"])
            else:
                af_s = vd_s = mq_s = "—"
            print(f"{vn:<40}  {r['gene']:<14}  {r['status']:<12}  {af_s:>8}  {vd_s:>10}  {mq_s:>4}")

    # ── AF distribution ────────────────────────────────────────────────────────
    pass_calls = [(vn, r) for vn, r in all_results.items() if r["status"] == "PASS"]
    if pass_calls:
        afs = sorted(r["calls"][0]["af"] for _, r in pass_calls)
        n = len(afs)
        print("\n" + sep)
        print(f"OBSERVED AF DISTRIBUTION FOR {n} PASS CALLS (target VAF=0.5%)")
        print(sep)
        for pct in [5, 25, 50, 75, 95]:
            idx = max(0, min(n - 1, int(pct / 100 * n)))
            print(f"  p{pct:2d}: {afs[idx]:.4f}")
        print(f"  mean: {statistics.mean(afs):.4f}  stdev: {statistics.stdev(afs):.4f}")
        print(f"  min:  {afs[0]:.4f}   max:  {afs[-1]:.4f}")


if __name__ == "__main__":
    main()
