#!/usr/bin/env bash
# ============================================================
# STEP 07a - Count replicate-level reads over D/V union peaks
# Peak regions:
#   D/V union confident peaks generated in Step 06
#
# Purpose:
#   Count reads/fragments from each biological replicate over
#   exactly the same genomic regions.
#
# Output:
#   count matrix for downstream DESeq2 analysis
#
# NOTE:
#   This step does NOT use merged BAM files.
#   Differential analysis requires individual biological
#   replicates.
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_confident_union_peaks
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths -----------------------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

bam_dir="$out/bam"

# D/V union peaks generated in Step 06
union_peaks="$out/peaks/06_confident_signal_deeptools/LNCD_CHIRP-DV_union_peaks.bed"

# Step 07 output
step7_dir="$out/differential_binding/07_D_vs_V"
count_dir="$step7_dir/counts"

mkdir -p "$count_dir"

threads=16

# ============================================================
# Check required programs
# ============================================================

for cmd in bedtools samtools; do

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found."
        exit 1
    fi

done

# ============================================================
# Check union peak file
# ============================================================

if [[ ! -s "$union_peaks" ]]; then

    echo "ERROR: Cannot find union peak file:"
    echo "$union_peaks"
    exit 1

fi

echo ""
echo "========================================"
echo "STEP 07A"
echo "Count reads over D/V union peaks"
echo "========================================"

echo ""
echo "Union peak file:"
echo "$union_peaks"

echo ""
echo "Number of union peaks:"
wc -l "$union_peaks"


# ============================================================
# A) Create numbered union peak BED
#
# BED4:
#
# chr    start    end    peak_ID
#
# Example:
# chr1   1000     1500   peak_1
# ============================================================

awk 'BEGIN{OFS="\t"}
     {
         print $1,$2,$3,"peak_"NR
     }' \
    "$union_peaks" \
    > "$count_dir/DV_union_peaks_numbered.bed"

numbered_peaks="$count_dir/DV_union_peaks_numbered.bed"

echo ""
echo "Numbered union peaks:"
echo "$numbered_peaks"

# ============================================================
# B) Define individual biological replicate BAM files
#
# IMPORTANT:
# Use FILTERED individual BAMs from Step 3.
#
# Do NOT use merged BAMs here.
# ============================================================

D1="$bam_dir/CHIRP-D1_S10_L002_filt.bam"
D2="$bam_dir/CHIRP-D2_S11_L002_filt.bam"
D3="$bam_dir/CHIRP-D3_S12_L002_filt.bam"

V1="$bam_dir/CHIRP-V1_S7_L002_filt.bam"
V2="$bam_dir/CHIRP-V2_S8_L002_filt.bam"
V3="$bam_dir/CHIRP-V3_S9_L002_filt.bam"


# ============================================================
# Check BAM files
# ============================================================

for bam in \
    "$D1" "$D2" "$D3" \
    "$V1" "$V2" "$V3"
do

    if [[ ! -s "$bam" ]]; then

        echo "ERROR: Cannot find BAM:"
        echo "$bam"
        exit 1

    fi

done

# ============================================================
# C) Count fragments overlapping each union peak
#
# bedtools multicov:
#
#   -bams = individual replicate BAM files
#   -bed  = common union peak regions
#
# Because the same BED is used for all six samples, every row
# corresponds to exactly the same genomic region.
# ============================================================

echo ""
echo "========================================"
echo "Counting reads/fragments"
echo "========================================"

bedtools multicov \
    -bams \
        "$D1" \
        "$D2" \
        "$D3" \
        "$V1" \
        "$V2" \
        "$V3" \
    -bed "$numbered_peaks" \
    > "$count_dir/DV_union_multicov_raw.tsv"


# ============================================================
# D) Add column names
#
# Raw multicov output:
#
# chr start end peak_ID D1 D2 D3 V1 V2 V3
# ============================================================

{
    printf "chr\tstart\tend\tpeak_ID\tD1\tD2\tD3\tV1\tV2\tV3\n"

    cat "$count_dir/DV_union_multicov_raw.tsv"

} > "$count_dir/DV_union_counts.tsv"


# ============================================================
# E) Generate count-only matrix for DESeq2
#
# Output:
#
# peak_ID   D1   D2   D3   V1   V2   V3
# ============================================================

awk 'BEGIN{OFS="\t"}
     NR==1 {
         print "peak_ID","D1","D2","D3","V1","V2","V3"
         next
     }
     {
         print $4,$5,$6,$7,$8,$9,$10
     }' \
    "$count_dir/DV_union_counts.tsv" \
    > "$count_dir/DV_union_DESeq2_counts.tsv"


# ============================================================
# F) Basic QC
# ============================================================

echo ""
echo "========================================"
echo "Count matrix QC"
echo "========================================"

echo ""
echo "Number of union regions:"
echo "$(($(wc -l < "$count_dir/DV_union_DESeq2_counts.tsv") - 1))"

echo ""
echo "Total counts in each sample:"

awk '
NR>1 {
    D1 += $2
    D2 += $3
    D3 += $4
    V1 += $5
    V2 += $6
    V3 += $7
}
END {
    print "D1:",D1
    print "D2:",D2
    print "D3:",D3
    print "V1:",V1
    print "V2:",V2
    print "V3:",V3
}' "$count_dir/DV_union_DESeq2_counts.tsv"


# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 07A COMPLETE"
echo "========================================"

echo ""
echo "Output directory:"
echo "$count_dir"

echo ""
echo "Main output for DESeq2:"
echo "$count_dir/DV_union_DESeq2_counts.tsv"

echo ""
echo "Counts with genomic coordinates:"
echo "$count_dir/DV_union_counts.tsv"

echo ""
echo "Numbered peak BED:"
echo "$numbered_peaks"