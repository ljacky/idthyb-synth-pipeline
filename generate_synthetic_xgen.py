#!/usr/bin/env python3
"""
generate_synthetic_xgen.py

Simulate paired-end IDT xGen Duplex UMI hybridization capture sequencing on the
WISDOM v3 xGen panel (251 targets, 14 chromosomes, ~77 kbp).

Chemistry model:
  - UMI structure: 9M1S+T (9 bp molecular barcode + 1 bp T-linker + insert).
    Both R1 and R2 carry the same 9 bp UMI — matches fgbio --strategy Adjacency.
  - Fragments: variable length (Poisson mean 220 bp, clamped 140–400 bp),
    starting uniformly within the extended target window.
  - Library complexity (Lander-Waterman): N_unique unique molecules per target,
    each duplicated geometrically at rate --dup-rate. Expected raw depth =
    N_unique / (1 - dup_rate).
  - Sequencing errors: position-dependent Q-score model (same as rhampseqtest).
  - Somatic variants spiked in at molecule level: a VAF fraction of molecules
    whose fragment overlaps the variant position carry the alt allele.
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
LINKER_LEN  = 1           # the T-linker skipped by fgbio 9M1S+T
INSERT_LEN  = READ_LEN - UMI_LEN - LINKER_LEN   # 140 bp of insert per read end
FRAG_MEAN   = 220         # mean fragment length (bp)
FRAG_MIN    = INSERT_LEN  # minimum usable fragment
FRAG_MAX    = 400
FRAG_WINDOW = 50          # extension beyond target edges for fragment start sampling

ADAPTER_FWD = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_REV = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

ERROR_PROFILE = [
    (110, "F", 0.0002),
    (25,  "?", 0.001),
    (15,  "7", 0.0063),
]

INSTRUMENT = "NB551234"
RUN_ID     = 44
FLOWCELL   = "HYYYYCGXY"
LANE       = 1

DROPOUT_THRESHOLDS = [100, 500, 1000]

# ─────────────────────────────────── variant / sample definitions ─────────────

# 0-based genomic coordinates on the + strand.
# ref is filled in by verify_refs(); alt must differ from actual ref.
VARIANTS = {
    "KRAS_G12D":    {"chrom": "chr12", "pos": 25245350,  "ref": None, "alt": "T"},
    "PIK3CA_H1047R":{"chrom": "chr3",  "pos": 179234297, "ref": None, "alt": "G"},
    "PIK3CA_E545K": {"chrom": "chr3",  "pos": 179218302, "ref": None, "alt": "A"},
    "TP53_R175H":   {"chrom": "chr17", "pos": 7675088,   "ref": None, "alt": "T"},
    "EGFR_L858R":   {"chrom": "chr7",  "pos": 55191821,  "ref": None, "alt": "G"},
    "NOTCH1_P2415L":{"chrom": "chr9",  "pos": 136496518, "ref": None, "alt": "T"},
}

SAMPLES = [
    {"name": "sample1", "index": "ATCACG",
     "variants": {"KRAS_G12D": 0.010, "PIK3CA_H1047R": 0.010}},
    {"name": "sample2", "index": "CGATGT",
     "variants": {"PIK3CA_E545K": 0.005, "NOTCH1_P2415L": 0.005}},
    {"name": "sample3", "index": "TTAGGC",
     "variants": {"TP53_R175H": 0.020}},
    {"name": "sample4", "index": "TGACCA",
     "variants": {"EGFR_L858R": 0.003, "TP53_R175H": 0.010}},
    {"name": "sample5", "index": "ACAGTG", "variants": {}},
    {"name": "sample6", "index": "GCCAAT", "variants": {}},
    {"name": "sample7", "index": "CAGATC", "variants": {}},
    {"name": "sample8", "index": "ACTTGA", "variants": {}},
]

# ─────────────────────────────────── helpers ─────────────────────────────────

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def rev_comp(seq):
    return seq.translate(COMP)[::-1]

def load_bed(path):
    """Load BED; prepend 'chr' if chromosome name lacks it."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
            rows.append((chrom, int(f[1]), int(f[2]), f[3] if len(f) > 3 else "."))
    return rows

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

# ─────────────────────────────────── reference verification ──────────────────

def verify_refs(variants, fasta):
    print("Verifying variant positions against reference...", flush=True)
    for vname, v in variants.items():
        ref_base = fasta.fetch(v["chrom"], v["pos"], v["pos"] + 1).upper()
        v["ref"] = ref_base
        if ref_base == v["alt"]:
            raise ValueError(
                f"{vname}: ref == alt == {ref_base!r} at {v['chrom']}:{v['pos']+1}. "
                "Choose a different alt allele.")
        print(f"  {vname}: {v['chrom']}:{v['pos']+1}  {ref_base}>{v['alt']}  OK", flush=True)

# ─────────────────────────────────── variant → target mapping ────────────────

def assign_variants_to_targets(variants, targets):
    mapping = defaultdict(list)
    for vname, v in variants.items():
        for i, (chrom, start, end, _name) in enumerate(targets):
            if chrom == v["chrom"] and start <= v["pos"] < end:
                mapping[vname].append(i)
        if not mapping[vname]:
            print(f"  WARNING: {vname} at {v['chrom']}:{v['pos']+1} "
                  "not covered by any target region", file=sys.stderr)
    return mapping

# ─────────────────────────────────── RNG (stdlib, no numpy) ──────────────────

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
        """Number of trials until first success (min 1)."""
        if p_success <= 0:
            return 1
        return math.ceil(math.log(self._r.random()) / math.log(1.0 - p_success))

# ─────────────────────────────────── simulation ──────────────────────────────

def simulate_sample(sample_cfg, targets, variant_target_map,
                    fasta, target_depth, dup_rate, rng, outdir):
    t0 = time.time()
    name    = sample_cfg["name"]
    index   = sample_cfg["index"]
    s_num   = int(name.replace("sample", ""))
    vaf_map = sample_cfg["variants"]

    # Build per-target variant list: target_idx → [{name, vaf, pos, ref, alt}]
    target_variants = defaultdict(list)
    for vname, vaf in vaf_map.items():
        for tidx in variant_target_map.get(vname, []):
            v = VARIANTS[vname]
            target_variants[tidx].append({
                "name": vname, "vaf": vaf,
                "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
            })

    r1_path = os.path.join(outdir, f"{name}_S{s_num}_L001_R1_001.fastq.gz")
    r2_path = os.path.join(outdir, f"{name}_S{s_num}_L001_R2_001.fastq.gz")

    pairs_written = 0
    variant_counts = defaultdict(lambda: {"n_alt_reads": 0, "n_total_reads": 0})
    per_target_depth = []
    tile, x, y = 1101, 0, 0

    with gzip.open(r1_path, "wt", compresslevel=4) as f1, \
         gzip.open(r2_path, "wt", compresslevel=4) as f2:

        for tidx, (chrom, tstart, tend, tname) in enumerate(targets):
            target_len = tend - tstart
            n_reads_target = target_depth * target_len / INSERT_LEN
            n_unique = max(1, round(n_reads_target * (1.0 - dup_rate)))
            variants_here = target_variants.get(tidx, [])

            # Window for fragment start sampling: extend by FRAG_WINDOW on each side
            win_start = max(0, tstart - FRAG_WINDOW)
            win_end   = tend  # fragment start can be up to tend-1 and still overlap

            reads_this_target = 0

            for _ in range(n_unique):
                # Draw fragment geometry
                for _ in range(20):  # retry until fragment overlaps target
                    frag_start = rng.randint(win_start, win_end - 1)
                    frag_len = max(FRAG_MIN,
                                   min(FRAG_MAX, rng.poisson(FRAG_MEAN)))
                    frag_end = frag_start + frag_len
                    if frag_start < tend and frag_end > tstart:
                        break
                else:
                    frag_start = tstart
                    frag_len   = FRAG_MIN
                    frag_end   = frag_start + frag_len

                umi = random_umi(rng)

                # Fetch reference template
                fetch_start = frag_start
                fetch_end   = min(frag_end, fasta.get_reference_length(chrom))
                template = fasta.fetch(chrom, fetch_start, fetch_end).upper()

                # Apply variant alt alleles
                for v in variants_here:
                    variant_counts[v["name"]]["n_total_reads"] += 1
                    if fetch_start <= v["pos"] < fetch_end:
                        if rng.random() < v["vaf"]:
                            off = v["pos"] - fetch_start
                            template = template[:off] + v["alt"] + template[off+1:]
                            variant_counts[v["name"]]["n_alt_reads"] += 1

                # Template carries the correct molecule sequence (ref or alt).
                # Sequencing errors are introduced independently per copy so that
                # UMI consensus correctly suppresses read errors (each copy has
                # different random errors; only the true alt allele is consistent).
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

            per_target_depth.append(reads_this_target)

    # Coverage stats across targets
    def depth_stats(counts):
        n  = len(counts)
        mn = statistics.mean(counts) if counts else 0
        md = float(statistics.median(counts)) if counts else 0
        sd = statistics.pstdev(counts) if counts else 0
        cv = sd / mn if mn > 0 else 0.0
        return {
            "n_targets":               n,
            "mean":                    round(mn, 1),
            "median":                  round(md, 1),
            "stdev":                   round(sd, 1),
            "cv":                      round(cv, 3),
            "min":                     min(counts) if counts else 0,
            "max":                     max(counts) if counts else 0,
            "uniformity_0_2x_mean_pct": round(
                100 * sum(1 for c in counts if c >= 0.2 * mn) / n, 1) if n else 0,
            **{f"n_below_{t}x": sum(1 for c in counts if c < t)
               for t in DROPOUT_THRESHOLDS},
        }

    # Variant summary
    vcounts_out = {}
    for vname, vaf in vaf_map.items():
        n_alt = variant_counts[vname]["n_alt_reads"]
        n_tot = variant_counts[vname]["n_total_reads"]
        obs_vaf = n_alt / n_tot if n_tot > 0 else 0.0
        # Effective depth = reads overlapping all covering targets
        eff_depth = sum(
            per_target_depth[tidx]
            for tidx in variant_target_map.get(vname, [])
        )
        p_detect = 1.0 - (1.0 - vaf) ** eff_depth if eff_depth > 0 else 0.0
        vcounts_out[vname] = {
            "vaf_target":          vaf,
            "n_alt_reads":         n_alt,
            "n_total_reads":       n_tot,
            "observed_vaf":        round(obs_vaf, 4),
            "n_covering_targets":  len(variant_target_map.get(vname, [])),
            "effective_depth":     eff_depth,
            "expected_alt_reads":  round(eff_depth * vaf, 1),
            "p_detect_at_least_1": round(p_detect, 4),
        }

    return {
        "sample":         name,
        "index":          index,
        "r1":             os.path.basename(r1_path),
        "r2":             os.path.basename(r2_path),
        "pairs_written":  pairs_written,
        "runtime_sec":    round(time.time() - t0, 1),
        "coverage_stats": depth_stats(per_target_depth),
        "variant_read_counts": vcounts_out,
    }

# ─────────────────────────────────── main ────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed",
        default=os.path.expanduser(
            "~/seqtest/idtdesigns/Design Output wisdom_v3_final green/"
            "Targets-XGEN.7ECA75CAFC214C2E982211ED4AF76C73.g.bed"))
    p.add_argument("--ref",
        default=os.path.expanduser("~/ref/hg38_canonical/hg38_canonical.fa"))
    p.add_argument("--outdir",   default="fastq")
    p.add_argument("--target-depth", type=int, default=2000,
        help="Raw read pairs per base in target regions (default 2000)")
    p.add_argument("--dup-rate",  type=float, default=0.72,
        help="PCR duplicate fraction — Lander-Waterman model (default 0.72)")
    p.add_argument("--seed",      type=int,   default=42)
    p.add_argument("--samples",   nargs="*",
        help="Subset of sample names to generate (default: all 8)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading target BED: {args.bed}", flush=True)
    targets = load_bed(args.bed)
    total_bp = sum(e - s for _, s, e, _ in targets)
    print(f"  {len(targets)} targets, {total_bp:,} bp total", flush=True)

    print(f"Opening reference: {args.ref}", flush=True)
    fasta = pysam.FastaFile(args.ref)

    verify_refs(VARIANTS, fasta)

    print("Assigning variants to targets...", flush=True)
    variant_target_map = assign_variants_to_targets(VARIANTS, targets)
    for vname, tidxs in variant_target_map.items():
        print(f"  {vname}: {len(tidxs)} target region(s)", flush=True)

    samples_to_run = args.samples or [s["name"] for s in SAMPLES]
    manifest_samples = []

    for s_cfg in SAMPLES:
        if s_cfg["name"] not in samples_to_run:
            continue
        # Per-sample seed: each sample is fully independent of the others.
        s_num = int(s_cfg["name"].replace("sample", ""))
        rng = RNG(args.seed * 1000 + s_num)
        print(f"\n[{s_cfg['name']}] Generating "
              f"(raw depth ~{args.target_depth}x, dup_rate={args.dup_rate}, "
              f"seed={args.seed * 1000 + s_num})...",
              flush=True)
        result = simulate_sample(
            s_cfg, targets, variant_target_map,
            fasta, args.target_depth, args.dup_rate, rng, args.outdir,
        )
        manifest_samples.append(result)

        cs  = result["coverage_stats"]
        vrc = result["variant_read_counts"]
        vstr = ", ".join(
            f"{vn} obs_vaf={d['observed_vaf']:.4f} eff_d={d['effective_depth']} "
            f"P(det)={d['p_detect_at_least_1']:.3f}"
            for vn, d in vrc.items()
        ) or "no spike-ins"
        print(f"  {result['pairs_written']:,} pairs  "
              f"mean_depth={cs['mean']}  CV={cs['cv']}  {result['runtime_sec']:.0f}s",
              flush=True)
        if vstr != "no spike-ins":
            print(f"  variants: {vstr}", flush=True)

    manifest = {
        "run": {
            "instrument":    INSTRUMENT,
            "run_id":        RUN_ID,
            "flowcell":      FLOWCELL,
            "kit":           "NextSeq 550 High Output 300-cycle (2x150 PE)",
            "lane":          LANE,
            "read_length":   READ_LEN,
            "chemistry":     "xGen Duplex UMI (9M1S+T single-strand grouping)",
            "panel":         "WISDOM v3 xGen final green — 251 targets",
            "panel_targets": len(targets),
            "panel_bp":      total_bp,
            "target_raw_depth": args.target_depth,
            "dup_rate":          args.dup_rate,
            "frag_mean_bp":      FRAG_MEAN,
            "frag_range_bp":     [FRAG_MIN, FRAG_MAX],
            "umi_len":           UMI_LEN,
            "seed":          args.seed,
            "reference":     "hg38_canonical (GRCh38, chr-prefixed, no ALT/random)",
            "error_profile": [{"cycles": n, "qchar": q, "p_err": p}
                              for n, q, p in ERROR_PROFILE],
            "adapter_r1_readthrough": ADAPTER_FWD,
            "adapter_r2_readthrough": ADAPTER_REV,
        },
        "variants": {
            vn: {
                "chrom":              v["chrom"],
                "pos_0based":         v["pos"],
                "pos_1based":         v["pos"] + 1,
                "ref":                v["ref"],
                "alt":                v["alt"],
                "n_covering_targets": len(variant_target_map.get(vn, [])),
            }
            for vn, v in VARIANTS.items()
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
