#!/usr/bin/env python3
"""
generate_sensitivity_xgen.py

Sensitivity benchmark: spike one variant per target region at 0.5% VAF.
Tests all 251 xGen capture targets across 4 positive samples + 1 negative control.

Sample design (by gene cluster):
  pos1: KMT2D                            (54 targets)
  pos2: NOTCH1 + FAT1                    (~60 targets)
  pos3: EP300 + NSD1 + CUL3             (~69 targets)
  pos4: TP53 + PTEN + CASP8 + AJUBA + KEAP1 + KRAS + HRAS + EGFR +
        PIK3R1 + CDKN2A + PIK3CA + TERT + NFE2L2 + FBXW7 + DPYD + other
  neg5: no spike-ins

Variant selection: one C>T (or G>A) transition per target at the position
nearest the target center, avoiding homopolymer runs ≥ 3.
Falls back to A>G or T>C if no C/G is found near center.

Chemistry: IDT xGen Duplex UMI 9M1S+T — variable-length fragments, Lander-Waterman
duplicate model. Each copy of a unique molecule gets independent sequencing errors.
"""

import argparse
import gzip
import json
import math
import os
import random
import statistics
import sys
import time
from collections import defaultdict

import pysam

# ─────────────────────────────────── constants ───────────────────────────────

READ_LEN    = 150
UMI_LEN     = 9
LINKER_LEN  = 1
INSERT_LEN  = READ_LEN - UMI_LEN - LINKER_LEN  # 140 bp

FRAG_MEAN   = 220
FRAG_MIN    = INSERT_LEN
FRAG_MAX    = 400
FRAG_WINDOW = 50   # extension beyond target edges for fragment start sampling

ADAPTER_FWD = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_REV = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

ERROR_PROFILE = [
    (110, "F", 0.0002),
    (25,  "?", 0.001),
    (15,  "7", 0.0063),
]

INSTRUMENT = "NB551234"
RUN_ID     = 55
FLOWCELL   = "HSENSCGXY"
LANE       = 1

TARGET_VAF = 0.005   # 0.5% for all positive spike-ins

# ─────────────────────────────────── gene / sample assignment ─────────────────

# gene label → sample index (0-based)
SAMPLE_GENE_MAP = {
    "KMT2D":  0,
    "NOTCH1": 1, "FAT1": 1,
    "EP300":  2, "NSD1": 2, "CUL3": 2,
    # everything else → 3 (pos4)
}

SAMPLE_CONFIGS = [
    {"name": "pos1", "index": "ATCACG", "label": "KMT2D"},
    {"name": "pos2", "index": "CGATGT", "label": "NOTCH1+FAT1"},
    {"name": "pos3", "index": "TTAGGC", "label": "EP300+NSD1+CUL3"},
    {"name": "pos4", "index": "TGACCA",
     "label": "TP53+PTEN+CASP8+AJUBA+KEAP1+KRAS+HRAS+EGFR+PIK3R1+CDKN2A+PIK3CA+TERT+NFE2L2+FBXW7+DPYD+other"},
    {"name": "neg5", "index": "ACAGTG", "label": "negative_control"},
]

# ─────────────────────────────────── helpers ─────────────────────────────────

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def rev_comp(seq):
    return seq.translate(COMP)[::-1]

def load_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
            rows.append((chrom, int(f[1]), int(f[2]), f[3] if len(f) > 3 else "."))
    return rows

def gene_from_name(name):
    """Extract gene label from BED name like 'KMT2D_exon5' or 'DPYD_2A_rs3918290_IVS14plus1'."""
    # Strip trailing _exonN or _exonNNN
    import re
    m = re.match(r'^([A-Za-z0-9]+(?:_[A-Za-z]+\d*)?)', name)
    if not m:
        return "other"
    prefix = m.group(1)
    # Known multi-word gene names that should be kept together
    known = {"CDKN2Abeta": "CDKN2A", "TERT_promoter": "TERT"}
    for k, v in known.items():
        if name.startswith(k):
            return v
    # Single-word gene names: first underscore-delimited token
    return name.split("_")[0]

def random_umi(rng, length=UMI_LEN):
    return "".join(rng.choice("ACGT") for _ in range(length))

def intro_errors(seq, rng):
    bases = list(seq)
    qual  = []
    pos   = 0
    for n_cyc, qchar, p_err in ERROR_PROFILE:
        for _ in range(n_cyc):
            if pos >= len(bases):
                break
            qual.append(qchar)
            if bases[pos] not in "Nn" and rng.random() < p_err:
                bases[pos] = rng.choice([b for b in "ACGT" if b != bases[pos]])
            pos += 1
    return "".join(bases), "".join(qual)

def pad_adapter(seq, is_r2):
    if len(seq) >= INSERT_LEN:
        return seq[:INSERT_LEN]
    pad = (ADAPTER_REV if is_r2 else ADAPTER_FWD)[:INSERT_LEN - len(seq)]
    return seq + pad

def make_qname(tile, x, y):
    return f"{INSTRUMENT}:{RUN_ID}:{FLOWCELL}:{LANE}:{tile}:{x}:{y}"

class RNG:
    def __init__(self, seed):
        self._r = random.Random(seed)

    def random(self):
        return self._r.random()

    def choice(self, seq):
        return self._r.choice(seq)

    def randint(self, a, b):
        return self._r.randint(a, b)

    def poisson(self, lam):
        if lam > 30:
            return max(0, round(self._r.gauss(lam, math.sqrt(lam))))
        L = math.exp(-lam)
        k, p = 0, 1.0
        while p > L:
            k += 1
            p *= self._r.random()
        return k - 1

    def geometric(self, p_success):
        if p_success <= 0:
            return 1
        return math.ceil(math.log(self._r.random()) / math.log(1.0 - p_success))

# ─────────────────────────────────── variant selection ───────────────────────

def pick_variant_for_target(chrom, tstart, tend, fasta):
    """
    Pick one SNV near the target center, preferring C>T or G>A.
    Avoids homopolymer runs ≥ 3 and Ns. Falls back to A>G or T>C.
    Returns (pos_0based, ref_base, alt_base) or raises RuntimeError.
    """
    tlen   = tend - tstart
    seq    = fasta.fetch(chrom, tstart, tend).upper()
    center = tlen // 2

    preferred_fallback = None

    for offset in range(tlen):
        for sign in ([0] if offset == 0 else [1, -1]):
            idx = center + sign * offset
            if idx < 2 or idx >= tlen - 2:
                continue
            ref_base = seq[idx]
            if ref_base == 'N':
                continue
            window = seq[max(0, idx - 2): idx + 3]
            if any(window.count(b * 3) > 0 for b in 'ACGT'):
                continue
            pos = tstart + idx
            if ref_base == 'C':
                return pos, 'C', 'T'
            if ref_base == 'G':
                return pos, 'G', 'A'
            if preferred_fallback is None:
                alt = {'A': 'G', 'T': 'C'}.get(ref_base, 'G')
                preferred_fallback = (pos, ref_base, alt)

    if preferred_fallback:
        return preferred_fallback

    raise RuntimeError(
        f"Cannot find suitable variant in {chrom}:{tstart}-{tend} "
        "(all N or homopolymer?)")

# ─────────────────────────────────── variant assignment ──────────────────────

def build_variant_assignments(targets, fasta):
    """
    One variant per target, assigned to a sample by gene.
    Returns:
      variants        – {vname: {chrom, pos, ref, alt, gene}}
      target_variants – {target_idx: [{name, vaf, pos, ref, alt}]}
      sample_variants – list of 5 dicts {vname: TARGET_VAF}
      gene_labels     – [gene_label per target]
    """
    variants        = {}
    target_variants = defaultdict(list)
    sample_variants = [{} for _ in SAMPLE_CONFIGS]
    gene_labels     = []

    for tidx, (chrom, tstart, tend, tname) in enumerate(targets):
        gene = gene_from_name(tname)
        gene_labels.append(gene)
        sample_idx = SAMPLE_GENE_MAP.get(gene, 3)

        pos, ref, alt = pick_variant_for_target(chrom, tstart, tend, fasta)
        vname = f"{gene}_{chrom}:{pos + 1}"
        if vname in variants:
            vname = f"{vname}_{tidx}"

        variants[vname] = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                           "gene": gene, "target_idx": tidx,
                           "sample": SAMPLE_CONFIGS[sample_idx]["name"]}
        target_variants[tidx].append({"name": vname, "vaf": TARGET_VAF,
                                      "pos": pos, "ref": ref, "alt": alt})
        sample_variants[sample_idx][vname] = TARGET_VAF

    return variants, target_variants, sample_variants, gene_labels

# ─────────────────────────────────── simulation ──────────────────────────────

def simulate_sample(sample_cfg, sample_idx, targets, target_variants,
                    fasta, target_depth, dup_rate, rng, outdir):
    t0    = time.time()
    name  = sample_cfg["name"]
    index = sample_cfg["index"]

    r1_path = os.path.join(outdir, f"{name}_S{sample_idx + 1}_L001_R1_001.fastq.gz")
    r2_path = os.path.join(outdir, f"{name}_S{sample_idx + 1}_L001_R2_001.fastq.gz")

    pairs_written  = 0
    variant_counts = defaultdict(lambda: {"n_alt": 0, "n_total": 0})
    per_target_raw = []  # raw read pairs written per target
    tile, x, y = 1101, 0, 0

    with gzip.open(r1_path, "wt", compresslevel=4) as f1, \
         gzip.open(r2_path, "wt", compresslevel=4) as f2:

        for tidx, (chrom, tstart, tend, _tname) in enumerate(targets):
            target_len = tend - tstart
            n_reads_target = target_depth * target_len / INSERT_LEN
            n_unique = max(1, round(n_reads_target * (1.0 - dup_rate)))
            variants_here = target_variants.get(tidx, [])

            win_start = max(0, tstart - FRAG_WINDOW)
            win_end   = tend

            reads_this_target = 0

            for _ in range(n_unique):
                for _ in range(20):
                    frag_start = rng.randint(win_start, win_end - 1)
                    frag_len   = max(FRAG_MIN, min(FRAG_MAX, rng.poisson(FRAG_MEAN)))
                    frag_end   = frag_start + frag_len
                    if frag_start < tend and frag_end > tstart:
                        break
                else:
                    frag_start = tstart
                    frag_len   = FRAG_MIN
                    frag_end   = frag_start + frag_len

                umi = random_umi(rng)

                fetch_start = frag_start
                fetch_end   = min(frag_end, fasta.get_reference_length(chrom))
                template = fasta.fetch(chrom, fetch_start, fetch_end).upper()

                for v in variants_here:
                    if fetch_start <= v["pos"] < fetch_end:
                        # This molecule covers the variant position
                        variant_counts[v["name"]]["n_total"] += 1
                        if rng.random() < v["vaf"]:
                            off = v["pos"] - fetch_start
                            template = template[:off] + v["alt"] + template[off + 1:]
                            variant_counts[v["name"]]["n_alt"] += 1

                r1_template = pad_adapter(template[:INSERT_LEN], False)
                r2_template = pad_adapter(rev_comp(template[-INSERT_LEN:]), True)

                n_copies = min(rng.geometric(1.0 - dup_rate), 20)
                for _ in range(n_copies):
                    r1_seq, r1_qual = intro_errors(r1_template, rng)
                    r2_seq, r2_qual = intro_errors(r2_template, rng)
                    r1_raw = umi + "T" + r1_seq
                    r2_raw = umi + "T" + r2_seq
                    r1_qual_raw = "F" * (UMI_LEN + LINKER_LEN) + r1_qual
                    r2_qual_raw = "F" * (UMI_LEN + LINKER_LEN) + r2_qual
                    qname = make_qname(tile, x, y)
                    f1.write(f"@{qname} 1:N:0:{index}\n{r1_raw}\n+\n{r1_qual_raw}\n")
                    f2.write(f"@{qname} 2:N:0:{index}\n{r2_raw}\n+\n{r2_qual_raw}\n")
                    pairs_written += 1
                    reads_this_target += 1
                    y += 1
                    if y > 9999:
                        y = 0; x += 1
                    if x > 9999:
                        x = 0; tile += 1

            per_target_raw.append(reads_this_target)

    # Coverage stats
    def depth_stats(counts):
        n  = len(counts)
        mn = statistics.mean(counts) if counts else 0
        sd = statistics.pstdev(counts) if counts else 0
        return {
            "n_targets": n,
            "mean":      round(mn, 1),
            "median":    round(float(statistics.median(counts)), 1) if counts else 0,
            "cv":        round(sd / mn, 3) if mn > 0 else 0,
            "min":       min(counts) if counts else 0,
            "max":       max(counts) if counts else 0,
        }

    # Variant summary — n_total counts only molecules that actually cover the position
    vcounts_out = {}
    for vname, vc in variant_counts.items():
        n_alt = vc["n_alt"]
        n_tot = vc["n_total"]
        obs_vaf = n_alt / n_tot if n_tot > 0 else 0.0
        vcounts_out[vname] = {
            "vaf_target":          TARGET_VAF,
            "n_alt_reads":         n_alt,
            "n_covering_reads":    n_tot,
            "observed_vaf":        round(obs_vaf, 4),
            "expected_alt_reads":  round(n_tot * TARGET_VAF, 1),
        }

    return {
        "sample":              name,
        "label":               sample_cfg["label"],
        "index":               index,
        "r1":                  os.path.basename(r1_path),
        "r2":                  os.path.basename(r2_path),
        "pairs_written":       pairs_written,
        "runtime_sec":         round(time.time() - t0, 1),
        "coverage_stats":      depth_stats(per_target_raw),
        "n_spikeins":          len(variant_counts),
        "variant_read_counts": vcounts_out,
    }

# ─────────────────────────────────── main ────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed",
        default=os.path.expanduser("~/idthybtest/panel.bed"))
    p.add_argument("--ref",
        default=os.path.expanduser("~/ref/hg38_canonical/hg38_canonical.fa"))
    p.add_argument("--outdir",       default="fastq_sens")
    p.add_argument("--target-depth", type=int, default=20000,
        help="Raw reads per base per target region (default 20000)")
    p.add_argument("--dup-rate",     type=float, default=0.72)
    p.add_argument("--seed",         type=int,   default=123)
    p.add_argument("--samples",      nargs="*",
        help="Subset of sample names (default: all 5)")
    return p.parse_args()


def main():
    args = parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading target BED: {args.bed}", flush=True)
    targets = load_bed(args.bed)
    total_bp = sum(e - s for _, s, e, _ in targets)
    print(f"  {len(targets)} targets, {total_bp:,} bp total", flush=True)

    print(f"Opening reference: {args.ref}", flush=True)
    fasta = pysam.FastaFile(args.ref)

    print("Selecting variants (one per target, 0.5% VAF)...", flush=True)
    variants, target_variants, sample_variants, gene_labels = \
        build_variant_assignments(targets, fasta)

    total_variants = sum(len(sv) for sv in sample_variants)
    print(f"  {total_variants} spike-in variants across {len(targets)} targets", flush=True)

    # Per-sample breakdown by gene
    gene_counts = defaultdict(lambda: defaultdict(int))
    for tidx, gene in enumerate(gene_labels):
        sample_idx = SAMPLE_GENE_MAP.get(gene, 3)
        gene_counts[SAMPLE_CONFIGS[sample_idx]["name"]][gene] += 1

    for sc in SAMPLE_CONFIGS:
        gc = gene_counts[sc["name"]]
        if gc:
            breakdown = ", ".join(f"{g}={n}" for g, n in sorted(gc.items()))
            print(f"  {sc['name']} ({sc['label']}): {sum(gc.values())} — {breakdown}",
                  flush=True)
        else:
            print(f"  {sc['name']} ({sc['label']}): 0 (negative control)", flush=True)

    samples_to_run = set(args.samples or [sc["name"] for sc in SAMPLE_CONFIGS])
    manifest_samples = []

    for s_idx, s_cfg in enumerate(SAMPLE_CONFIGS):
        if s_cfg["name"] not in samples_to_run:
            continue
        # Per-sample seed for independence
        rng = RNG(args.seed * 1000 + s_idx)
        n_sv = len(sample_variants[s_idx])
        print(f"\n[{s_cfg['name']}] Generating "
              f"(target_depth={args.target_depth}x, dup_rate={args.dup_rate}, "
              f"{n_sv} spike-ins)...", flush=True)
        # Build per-sample target_variants: only spike in this sample's assigned variants.
        # The global target_variants maps every target to a variant; without filtering,
        # neg5 would get all 251 spike-ins and each positive sample would get all 251
        # instead of its own gene cluster.
        sv = sample_variants[s_idx]   # {vname: TARGET_VAF}, empty for neg5
        sample_tvars = defaultdict(list)
        for tidx, var_list in target_variants.items():
            for v in var_list:
                if v["name"] in sv:
                    sample_tvars[tidx].append(v)

        result = simulate_sample(
            s_cfg, s_idx, targets, sample_tvars,
            fasta, args.target_depth, args.dup_rate, rng, args.outdir,
        )
        manifest_samples.append(result)

        cs = result["coverage_stats"]
        print(f"  {result['pairs_written']:,} pairs  "
              f"mean_depth={cs['mean']}  CV={cs['cv']}  "
              f"{result['runtime_sec']:.0f}s", flush=True)

    manifest = {
        "run": {
            "instrument":      INSTRUMENT,
            "run_id":          RUN_ID,
            "flowcell":        FLOWCELL,
            "read_length":     READ_LEN,
            "chemistry":       "xGen Duplex UMI (9M1S+T single-strand grouping)",
            "panel":           "WISDOM v3 xGen final green — 251 targets",
            "panel_targets":   len(targets),
            "panel_bp":        total_bp,
            "target_raw_depth": args.target_depth,
            "dup_rate":         args.dup_rate,
            "target_vaf":       TARGET_VAF,
            "seed":             args.seed,
            "reference": "hg38_canonical (GRCh38, chr-prefixed, no ALT/random)",
        },
        "variants": {
            vn: {
                "chrom":      v["chrom"],
                "pos_0based": v["pos"],
                "pos_1based": v["pos"] + 1,
                "ref":        v["ref"],
                "alt":        v["alt"],
                "gene":       v["gene"],
                "sample":     v["sample"],
            }
            for vn, v in variants.items()
        },
        "samples": manifest_samples,
    }

    manifest_path = os.path.join(args.outdir, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\nManifest written to {manifest_path}")
    print("Done.")


if __name__ == "__main__":
    main()
