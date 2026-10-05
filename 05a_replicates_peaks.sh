#!/usr/bin/env bash
# =============================================================================
# STEP 5a - Call peaks for individual ChIRP-seq replicates
#
# Apply the same criteria as pooled peak calling:
#   MACS2 q < 0.01
#   FC >= 5
#   remove blacklist regions
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 05a_replicates_peaks.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_replicates
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"
bam_dir="$out/bam"
peak_dir="$out/peaks/replicates"
blacklist="$out/ref/hg38-blacklist.v2.bed"
fc_min=5

mkdir -p "$peak_dir"
# ============================================================
# D condition
# ============================================================

for rep in 1 2 3; do

    # Match the actual BAM filenames from Step 3
    case "$rep" in
        1)
            chirp_bam="$bam_dir/CHIRP-D1_S10_L002_filt.bam"
            input_bam="$bam_dir/Input-D1_S4_L002_filt.bam"
            ;;
        2)
            chirp_bam="$bam_dir/CHIRP-D2_S11_L002_filt.bam"
            input_bam="$bam_dir/Input-D2_S5_L002_filt.bam"
            ;;
        3)
            chirp_bam="$bam_dir/CHIRP-D3_S12_L002_filt.bam"
            input_bam="$bam_dir/Input-D3_S6_L002_filt.bam"
            ;;
    esac

    name="LNCD_CHIRP-D${rep}"

    echo "========================================"
    echo "MACS2: CHIRP-D${rep} vs Input-D${rep}"
    echo "========================================"

    macs2 callpeak \
        -t "$chirp_bam" \
        -c "$input_bam" \
        -f BAMPE \
        -g hs \
        -q 0.01 \
        --keep-dup all \
        --name "$name" \
        --outdir "$peak_dir" \
        2> "$peak_dir/${name}_macs2.log"

    np="$peak_dir/${name}_peaks.narrowPeak"

    # FC >= 5
    awk -v fc="$fc_min" \
        'BEGIN{OFS="\t"} $7 >= fc' \
        "$np" \
        > "$peak_dir/${name}_FC${fc_min}.narrowPeak"

    # Remove blacklist regions
    bedtools intersect -v \
        -a "$peak_dir/${name}_FC${fc_min}.narrowPeak" \
        -b "$blacklist" \
        > "$peak_dir/${name}_FC${fc_min}_noBL.narrowPeak"

    echo "D${rep} raw:     $(wc -l < "$np")"
    echo "D${rep} FC>=5:   $(wc -l < "$peak_dir/${name}_FC${fc_min}.narrowPeak")"
    echo "D${rep} noBL:    $(wc -l < "$peak_dir/${name}_FC${fc_min}_noBL.narrowPeak")"
    echo ""

done
# ============================================================
# V condition
# IMPORTANT:
# CHIRP uses uppercase V
# Input uses lowercase v
# ============================================================

for rep in 1 2 3; do

    case "$rep" in
        1)
            chirp_bam="$bam_dir/CHIRP-V1_S7_L002_filt.bam"
            input_bam="$bam_dir/Input-v1_S1_L002_filt.bam"
            ;;
        2)
            chirp_bam="$bam_dir/CHIRP-V2_S8_L002_filt.bam"
            input_bam="$bam_dir/Input-v2_S2_L002_filt.bam"
            ;;
        3)
            chirp_bam="$bam_dir/CHIRP-V3_S9_L002_filt.bam"
            input_bam="$bam_dir/Input-v3_S3_L002_filt.bam"
            ;;
    esac

    name="LNCD_CHIRP-V${rep}"

    echo "========================================"
    echo "MACS2: CHIRP-V${rep} vs Input-v${rep}"
    echo "========================================"

    macs2 callpeak \
        -t "$chirp_bam" \
        -c "$input_bam" \
        -f BAMPE \
        -g hs \
        -q 0.01 \
        --keep-dup all \
        --name "$name" \
        --outdir "$peak_dir" \
        2> "$peak_dir/${name}_macs2.log"

    np="$peak_dir/${name}_peaks.narrowPeak"

    # FC >= 5
    awk -v fc="$fc_min" \
        'BEGIN{OFS="\t"} $7 >= fc' \
        "$np" \
        > "$peak_dir/${name}_FC${fc_min}.narrowPeak"

    # Remove blacklist regions
    bedtools intersect -v \
        -a "$peak_dir/${name}_FC${fc_min}.narrowPeak" \
        -b "$blacklist" \
        > "$peak_dir/${name}_FC${fc_min}_noBL.narrowPeak"

    echo "V${rep} raw:     $(wc -l < "$np")"
    echo "V${rep} FC>=5:   $(wc -l < "$peak_dir/${name}_FC${fc_min}.narrowPeak")"
    echo "V${rep} noBL:    $(wc -l < "$peak_dir/${name}_FC${fc_min}_noBL.narrowPeak")"
    echo ""

done

echo "========================================"
echo "STEP 5a complete"
echo "========================================"
echo "Individual replicate peaks are in:"
echo "$peak_dir"
