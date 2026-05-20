#!/usr/bin/env python3
"""
benchmark.py — xGen UMI hybridization capture pipeline validation.

Reports: coverage summary, spike-in detection table, FP counts in negatives,
and simulator truth vs observed AF for positive controls.

Usage:
  cd ~/idthybtest
  python3 benchmark.py
"""

import gzip
import json
import os
import subprocess
import sys

SAMPLES  = [f"sample{i}" for i in range(1, 9)]
NEG_CTRL = [f"sample{i}" for i in range(5, 9)]

SPIKEIN_TRUTH = [
    {"sample": "sample1", "name": "KRAS_G12D",     "chrom": "chr12", "pos": 25245351,  "vaf": 0.010},
    {"sample": "sample1", "name": "PIK3CA_H1047R", "chrom": "chr3",  "pos": 179234298, "vaf": 0.010},
    {"sample": "sample2", "name": "PIK3CA_E545K",   "chrom": "chr3",  "pos": 179218303, "vaf": 0.005},
    {"sample": "sample2", "name": "NOTCH1_P2415L", "chrom": "chr9",  "pos": 136496519, "vaf": 0.005},
    {"sample": "sample3", "name": "TP53_R175H",    "chrom": "chr17", "pos": 7675089,   "vaf": 0.020},
    {"sample": "sample4", "name": "EGFR_L858R",    "chrom": "chr7",  "pos": 55191822,  "vaf": 0.003},
    {"sample": "sample4", "name": "TP53_R175H",    "chrom": "chr17", "pos": 7675089,   "vaf": 0.010},
]


def bcftools_query_pos(vcf, chrom, pos):
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
            "chrom": fields[0], "pos": int(fields[1]),
            "ref": fields[3], "alt": fields[4],
            "filter": fields[6],
            "af": float(info.get("AF", 0)),
            "dp": int(info.get("DP", 0)),
            "vd": int(info.get("VD", 0)),
            "mq": info.get("MQ", "?"),
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


def load_manifest():
    path = "fastq/manifest.json"
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def collect_results():
    manifest = load_manifest()
    results = {}
    for s in SAMPLES:
        raw_vcf = f"out/{s}/{s}.vardict.raw.vcf.gz"
        con_vcf = f"out/{s}/{s}.vardict.consensus.vcf.gz"
        raw_cov = f"out/{s}/cov_raw.regions.bed.gz"
        con_cov = f"out/{s}/cov_consensus.regions.bed.gz"

        results[s] = {
            "raw_vcf_exists":  os.path.exists(raw_vcf),
            "con_vcf_exists":  os.path.exists(con_vcf),
            "raw_mean_depth":  mean_depth(raw_cov),
            "con_mean_depth":  mean_depth(con_cov),
            "raw_n_total":     count_vcf_rows(raw_vcf) if os.path.exists(raw_vcf) else None,
            "raw_n_pass":      count_vcf_rows(raw_vcf, pass_only=True) if os.path.exists(raw_vcf) else None,
            "con_n_total":     count_vcf_rows(con_vcf) if os.path.exists(con_vcf) else None,
            "con_n_pass":      count_vcf_rows(con_vcf, pass_only=True) if os.path.exists(con_vcf) else None,
            "spikeins": {},
        }

        sim_truth = {}
        if manifest:
            for ms in manifest["samples"]:
                if ms["sample"] == s:
                    sim_truth = ms.get("variant_read_counts", {})
                    break

        for t in SPIKEIN_TRUTH:
            if t["sample"] != s:
                continue
            vname = t["name"]
            calls = bcftools_query_pos(con_vcf, t["chrom"], t["pos"]) \
                if os.path.exists(con_vcf) else []
            sim = sim_truth.get(vname, {})
            results[s]["spikeins"][vname] = {
                "target_vaf":    t["vaf"],
                "sim_obs_vaf":   sim.get("observed_vaf"),
                "sim_alt_reads": sim.get("n_alt_reads"),
                "sim_eff_depth": sim.get("effective_depth"),
                "sim_p_detect":  sim.get("p_detect_at_least_1"),
                "calls":         calls,
            }

    return results


def print_report(results):
    sep = "-" * 100

    print("\n" + sep)
    print("COVERAGE SUMMARY (mean depth over panel.bed)")
    print(sep)
    print(f"{'Sample':<10}  {'Raw depth':>12}  {'Consensus depth':>16}  {'Fold drop':>10}")
    for s in SAMPLES:
        r = results[s]
        rd = f"{r['raw_mean_depth']:.0f}x" if r["raw_mean_depth"] else "pending"
        cd = f"{r['con_mean_depth']:.0f}x" if r["con_mean_depth"] else "pending"
        fd = ""
        if r["raw_mean_depth"] and r["con_mean_depth"] and r["con_mean_depth"] > 0:
            fd = f"{r['raw_mean_depth']/r['con_mean_depth']:.1f}x"
        print(f"{s:<10}  {rd:>12}  {cd:>16}  {fd:>10}")

    print("\n" + sep)
    print("SPIKE-IN DETECTION (consensus VCF)")
    print(sep)
    print(f"{'Variant (sample, VAF)':<30}  {'Result'}")
    for t in SPIKEIN_TRUTH:
        s, vname, vaf = t["sample"], t["name"], t["vaf"]
        label = f"{vname} ({s}, {vaf:.1%})"
        info  = results[s]["spikeins"].get(vname, {})
        calls = info.get("calls", [])
        if not results[s]["con_vcf_exists"]:
            cell = "pending"
        elif not calls:
            sim_alt = info.get("sim_alt_reads")
            cell = f"MISS  (sim alt={sim_alt})"
        else:
            c = calls[0]
            status = "PASS" if c["filter"] == "PASS" else f"FAIL:{c['filter']}"
            cell = f"{status}  AF={c['af']:.4f}  VD={c['vd']}/{c['dp']}  MQ={c['mq']}"
        print(f"{label:<30}  {cell}")

    print("\n" + sep)
    print("FALSE POSITIVES IN NEGATIVE CONTROLS (samples 5-8)")
    print(sep)
    print(f"{'Sample':<10}  {'Raw total/PASS':>18}  {'Consensus total/PASS':>22}")
    for s in NEG_CTRL:
        r = results[s]
        if r["raw_n_total"] is None:
            raw_str = "pending"
        else:
            raw_str = f"{r['raw_n_total']} / {r['raw_n_pass']}"
        if r["con_n_total"] is None:
            con_str = "pending"
        else:
            con_str = f"{r['con_n_total']} / {r['con_n_pass']}"
        print(f"{s:<10}  {raw_str:>18}  {con_str:>22}")

    print("\n" + sep)
    print("SIMULATOR TRUTH vs PIPELINE — POSITIVE CONTROLS")
    print(sep)
    print(f"{'Variant':<20}  {'Sim VAF':>8}  {'Sim alt/eff_d':>14}  "
          f"{'P(det)':>7}  {'Obs AF':>8}  {'VD/DP':>12}  {'Filter'}")
    for t in SPIKEIN_TRUTH:
        s, vname = t["sample"], t["name"]
        info  = results[s]["spikeins"].get(vname, {})
        calls = info.get("calls", [])
        sim_v = f"{info['sim_obs_vaf']:.4f}" if info.get("sim_obs_vaf") is not None else "—"
        sim_a = f"{info.get('sim_alt_reads','?')}/{info.get('sim_eff_depth','?')}"
        p_det = f"{info['sim_p_detect']:.3f}" if info.get("sim_p_detect") is not None else "—"
        if calls:
            c = calls[0]
            obs_af = f"{c['af']:.4f}"
            vd_dp  = f"{c['vd']}/{c['dp']}"
            filt   = c["filter"]
        else:
            obs_af, vd_dp, filt = "—", "—", "MISS"
        print(f"{vname:<20}  {sim_v:>8}  {sim_a:>14}  {p_det:>7}  "
              f"{obs_af:>8}  {vd_dp:>12}  {filt}")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    results = collect_results()
    print_report(results)


if __name__ == "__main__":
    main()
