#!/usr/bin/env bash
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh && conda activate ngs
START_TIME=$(date +%s)
trap 'echo ""; echo "=== TOTAL: $(( ($(date +%s) - START_TIME) / 60 )) min ($(( $(date +%s) - START_TIME )) sec) ==="' EXIT
SAMPLE="$1"
THREADS="${2:-8}"
FQDIR="${3:-$HOME/idthybtest/fastq}"
OUTBASE="${4:-$HOME/idthybtest/out}"

REF=~/ref/hg38_canonical/hg38_canonical.fa
BED=~/idthybtest/panel.bed
OUT="$OUTBASE/$SAMPLE"
mkdir -p "$OUT"

R1=$(ls $FQDIR/${SAMPLE}_S*_L001_R1_001.fastq.gz)
R2=$(ls $FQDIR/${SAMPLE}_S*_L001_R2_001.fastq.gz)

echo "[$SAMPLE] FastqToBam (UMI extraction; renaming reads via FIFOs)"
TMP=$(mktemp -d -t ngs-${SAMPLE}-XXXXXX)
trap 'rm -rf "$TMP"; echo ""; echo "=== TOTAL: $(( ($(date +%s) - START_TIME) / 60 )) min ($(( $(date +%s) - START_TIME )) sec) ==="' EXIT
mkfifo "$TMP/r1.fq" "$TMP/r2.fq"
python3 -c "
import fcntl, os
for p in ['$TMP/r1.fq', '$TMP/r2.fq']:
    fd = os.open(p, os.O_RDWR | os.O_NONBLOCK)
    fcntl.fcntl(fd, 1031, 1048576)
    os.close(fd)
"
( zcat "$R1" | awk 'BEGIN{FS=OFS=" "} NR%4==1 {$1=$1"_"(++c)} {print}' > "$TMP/r1.fq" ) &
PID1=$!
( zcat "$R2" | awk 'BEGIN{FS=OFS=" "} NR%4==1 {$1=$1"_"(++c)} {print}' > "$TMP/r2.fq" ) &
PID2=$!
fgbio -Xmx8g --compression 1 FastqToBam \
  --input "$TMP/r1.fq" "$TMP/r2.fq" \
  --read-structures 9M1S+T 9M1S+T \
  --sample "$SAMPLE" --library "$SAMPLE" \
  --output "$OUT/unmapped.bam"
wait $PID1 $PID2

echo "[$SAMPLE] Map + ZipperBams"
samtools fastq "$OUT/unmapped.bam" \
  | bwa-mem2 mem -t $THREADS -p -K 150000000 -Y "$REF" - \
  | fgbio -Xmx8g --compression 1 --async-io ZipperBams \
        --unmapped "$OUT/unmapped.bam" --ref "$REF" \
  | samtools sort -@ $THREADS -o "$OUT/mapped.bam" -
samtools index "$OUT/mapped.bam"

echo "[$SAMPLE] GroupReadsByUmi"
fgbio -Xmx8g GroupReadsByUmi \
  --input "$OUT/mapped.bam" \
  --strategy Adjacency \
  --edits 1 \
  --output "$OUT/grouped.bam" \
  --family-size-histogram "$OUT/family_sizes.tsv"

echo "[$SAMPLE] CallMolecularConsensusReads"
fgbio -Xmx8g CallMolecularConsensusReads \
  --input "$OUT/grouped.bam" \
  --min-reads 2 \
  --output "$OUT/consensus.unmapped.bam"

echo "[$SAMPLE] Remap consensus + FilterConsensusReads"
samtools fastq "$OUT/consensus.unmapped.bam" \
  | bwa-mem2 mem -t $THREADS -p -K 150000000 -Y "$REF" - \
  | fgbio -Xmx8g --compression 1 --async-io ZipperBams \
        --unmapped "$OUT/consensus.unmapped.bam" --ref "$REF" \
        --tags-to-reverse Consensus --tags-to-revcomp Consensus \
  | fgbio -Xmx8g --compression 1 FilterConsensusReads \
        --input /dev/stdin --output /dev/stdout --ref "$REF" \
        --min-reads 2 --min-base-quality 30 --max-base-error-rate 0.2 \
  | samtools sort -@ $THREADS -o "$OUT/consensus.bam" -
samtools index "$OUT/consensus.bam"

echo "[$SAMPLE] Coverage"
mosdepth --by "$BED" --thresholds 100,500,1000,5000 --no-per-base \
  "$OUT/cov_raw" "$OUT/mapped.bam"
mosdepth --by "$BED" --thresholds 50,100,500,1000 --no-per-base \
  "$OUT/cov_consensus" "$OUT/consensus.bam"

echo "[$SAMPLE] VarDict (raw)"
vardict-java -G "$REF" -f 0.001 -N "$SAMPLE" \
  -b "$OUT/mapped.bam" -c 1 -S 2 -E 3 -g 4 "$BED" \
  | teststrandbias.R \
  | var2vcf_valid.pl -N "$SAMPLE" -E -f 0.001 \
  | sed 's/SPANPAIR=\([0-9]*\)GT:/SPANPAIR=\1\tGT:/' \
  | bcftools reheader --fai "$REF.fai" \
  | bcftools sort -Ov 2>/dev/null \
  | bcftools norm -f "$REF" -m -both -Oz -o "$OUT/$SAMPLE.vardict.raw.vcf.gz" 2>/dev/null
tabix -f -p vcf "$OUT/$SAMPLE.vardict.raw.vcf.gz"

echo "[$SAMPLE] VarDict (consensus)"
vardict-java -G "$REF" -f 0.001 -N "$SAMPLE" \
  -b "$OUT/consensus.bam" -c 1 -S 2 -E 3 -g 4 "$BED" \
  | teststrandbias.R \
  | var2vcf_valid.pl -N "$SAMPLE" -E -f 0.001 \
  | sed 's/SPANPAIR=\([0-9]*\)GT:/SPANPAIR=\1\tGT:/' \
  | bcftools reheader --fai "$REF.fai" \
  | bcftools sort -Ov 2>/dev/null \
  | bcftools norm -f "$REF" -m -both -Oz -o "$OUT/$SAMPLE.vardict.consensus.vcf.gz" 2>/dev/null
tabix -f -p vcf "$OUT/$SAMPLE.vardict.consensus.vcf.gz"

echo "[$SAMPLE] DONE"
