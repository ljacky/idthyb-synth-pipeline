# idthybtest — xGen UMI hybridization capture synthetic data

Synthetic NGS dataset for the WISDOM v3 xGen final green panel (251 targets, ~77 kbp),
modelling IDT xGen Duplex UMI hybridization capture chemistry.

## Key design differences from ~/rhampseqtest

| Aspect | ~/rhampseqtest (RHAMPseq) | ~/idthybtest (xGen UMI) |
|--------|--------------------------|-------------------------|
| Chemistry | Amplicon PCR, no UMI | Hybridization capture + 9M1S+T UMI |
| Dedup | None — all reads retained | fgbio GroupReadsByUmi + CallMolecularConsensusReads |
| Duplicates modelled | n/a | Lander-Waterman ~72% dup rate |
| Fragment | Fixed amplicon length | Variable 140–400 bp overlapping target |
| Read structure | Insert only | [9 bp UMI][T][140 bp insert] per R1 and R2 |
| Pipeline | fastp → align → ampliconclip → VarDict | UMI extract → align → GroupReadsByUmi → consensus → VarDict |

## Working environment

| Item | Path |
|------|------|
| Project root | `~/idthybtest` |
| Conda env | `mamba activate ngs` |
| Pipeline script | `~/idthybtest/run_sample.sh` |
| Generator | `~/idthybtest/generate_synthetic_xgen.py` |
| Benchmark | `~/idthybtest/benchmark.py` |
| FASTQs | `~/idthybtest/fastq/` |
| Panel BED (chr-prefixed) | `~/idthybtest/panel.bed` (251 lines) |
| Source BED | `~/seqtest/idtdesigns/Design Output wisdom_v3_final green/Targets-XGEN.7ECA75CAFC214C2E982211ED4AF76C73.g.bed` |
| Reference | `~/ref/hg38_canonical/hg38_canonical.fa` |
| Output | `~/idthybtest/out/$SAMPLE/` |

## Spike-in design

| Sample | Variant | VAF | Difficulty |
|--------|---------|-----|------------|
| sample1 | KRAS G12D (chr12:25245351) | 1% | easy |
| sample1 | PIK3CA H1047R (chr3:179234298) | 1% | easy |
| sample2 | PIK3CA E545K (chr3:179218303) | 0.5% | medium |
| sample2 | NOTCH1 P2415L (chr9:136496519) | 0.5% | medium |
| sample3 | TP53 R175H (chr17:7675089) | 2% | easy — primary positive |
| sample4 | EGFR L858R (chr7:55191822) | 0.3% | LOD edge (~1.7 expected alt molecules) |
| sample4 | TP53 R175H (chr17:7675089) | 1% | easy — positive control alongside LOD edge |
| sample5–8 | — | — | negative controls |

All positions 1-based. EGFR L858R at 0.3% with ~560 unique molecules → ~1.7 expected alt
molecules at min-reads=2 consensus — likely MISS, which characterises the LOD for this
chemistry and depth (same function as HRAS Q61L in seqtest/rhampseqtest).

## Expected results

| Metric | Expected |
|--------|----------|
| Raw depth | ~2,000x per base (target depth) |
| Consensus depth | ~560x (72% dup → ~28% unique) |
| Raw VarDict FP count | 30–80 per sample |
| Consensus FP count | 0–5 for negatives |
| KRAS G12D, PIK3CA H1047R, PIK3CA E545K, TP53 R175H | HIT |
| EGFR L858R (0.3%) | MISS expected (LOD) |

## Running

```bash
cd ~/idthybtest
mamba activate ngs

# Generate FASTQs (~15 min, 8 samples)
python3 generate_synthetic_xgen.py

# Single sample validation
time ./run_sample.sh sample3 8 2>&1 | tee out/sample3.log

# Batch (3 concurrent × 8 threads)
parallel -j 3 './run_sample.sh {} 8 > out/{}.log 2>&1' ::: sample{1..8}

# Benchmark
python3 benchmark.py
```

## Pipeline notes

The run_sample.sh is identical to ~/seqtest/run_sample.sh except for FQDIR, BED, and OUT paths.
All seqtest gotchas apply: FIFO read renaming, SPANPAIR sed fix, bcftools reheader, etc.
See ~/seqtest/CLAUDE.md for the full gotcha list.

The BED note: source BED uses non-chr-prefixed chromosome names ("1", "2" …). The
panel.bed in this directory has "chr" prepended and is what run_sample.sh uses.
