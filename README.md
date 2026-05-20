# idthyb-synth-pipeline

Synthetic NGS dataset generator and UMI consensus pipeline for the **WISDOM v3 xGen final green** hybridization capture panel (251 targets, ~77 kbp, 14 chromosomes, 21 genes: AJUBA CASP8 CDKN2A CUL3 DPYD EGFR EP300 FAT1 FBXW7 HRAS KEAP1 KMT2D KRAS NFE2L2 NOTCH1 NSD1 PIK3CA PIK3R1 PTEN TERT TP53).

Models IDT xGen Duplex UMI hybridization capture chemistry (9M1S+T read structure, variable-length fragments, Lander-Waterman ~72% duplicate rate).

## Requirements

- Reference genome: `hg38_canonical.fa` (canonical chromosomes only — chr1–22, X, Y, M; no ALT/random contigs) with bwa-mem2 index
- Conda environment: see `environment.yml`

```bash
mamba env create -f environment.yml
mamba activate ngs
```

## Repository contents

| File | Description |
|------|-------------|
| `panel.bed` | Chr-prefixed 251-target BED (source: IDT xGen design, chr added) |
| `generate_synthetic_xgen.py` | 8-sample spike-in FASTQ generator (6 hotspot variants, 0.3–2% VAF) |
| `generate_sensitivity_xgen.py` | 251-target sensitivity panel generator (one variant per target, 0.5% VAF) |
| `run_sample.sh` | fgbio UMI consensus pipeline per sample |
| `run_sensitivity_pipeline.sh` | Parallel runner for 5-sample sensitivity panel |
| `benchmark.py` | Detection + FP summary for spike-in dataset |
| `benchmark_sensitivity.py` | Per-gene detection summary for sensitivity panel |
| `CLAUDE.md` | Full project handoff and pipeline notes |

## Spike-in dataset (8 samples, 2000× raw depth)

| Sample | Variant | VAF | Difficulty |
|--------|---------|-----|------------|
| sample1 | KRAS G12D + PIK3CA H1047R | 1% | easy |
| sample2 | PIK3CA E545K + NOTCH1 P2415L | 0.5% | medium |
| sample3 | TP53 R175H | 2% | easy — primary positive |
| sample4 | EGFR L858R + TP53 R175H | 0.3% / 1% | LOD edge + positive |
| sample5–8 | — | — | negative controls |

```bash
# Generate FASTQs (~2.5 hr, 5 samples at 20k depth)
python3 generate_synthetic_xgen.py

# Run pipeline (3 concurrent × 8 threads)
parallel -j 3 './run_sample.sh {} 8 > out/{}.log 2>&1' ::: sample{1..8}

# Benchmark
python3 benchmark.py
```

## Sensitivity panel (5 samples, 20 000× raw depth)

One C>T or G>A SNV per target at 0.5% VAF; 4 positive samples by gene cluster + 1 negative control.

| Sample | Gene cluster | Targets |
|--------|-------------|---------|
| pos1 | KMT2D | 54 |
| pos2 | NOTCH1 + FAT1 | 60 |
| pos3 | EP300 + NSD1 + CUL3 | 69 |
| pos4 | TP53 + PTEN + CASP8 + AJUBA + KEAP1 + KRAS + HRAS + EGFR + PIK3R1 + CDKN2A + PIK3CA + TERT + NFE2L2 + FBXW7 + DPYD | 68 |
| neg5 | negative control | 0 |

```bash
python3 generate_sensitivity_xgen.py
./run_sensitivity_pipeline.sh 8 3
python3 benchmark_sensitivity.py
```

## Pipeline overview

```
FastqToBam (fgbio, 9M1S+T UMI extraction)
  → bwa-mem2 align + ZipperBams
  → GroupReadsByUmi (Adjacency strategy, edits=1)
  → CallMolecularConsensusReads (min-reads=2)
  → remap + FilterConsensusReads (min-reads=2, min-BQ=30, max-error=0.2)
  → mosdepth (raw + consensus coverage)
  → VarDict (raw mapped.bam + consensus.bam, -f 0.001)
```

See `CLAUDE.md` for known gotchas (FIFO read renaming, SPANPAIR VCF fix, bcftools reheader).
