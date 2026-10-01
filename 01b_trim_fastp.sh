#!/usr/bin/env bash
# =============================================================================
# STEP 1b - Clean up the reads with fastp (trim adapters and poly-G tails).
#
# The raw FastQC showed leftover sequencing adapters and long runs of "G"
# (a common artifact of Illumina machines). fastp removes both so they don't
# confuse the alignment step. Because these reads are short (38 bp), we throw
# away anything shorter than 20 bp after trimming.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 01b_trim_fastp.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_fastp
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"  # output directory for this project
fastq_dir="$out/IGM_rawdata"           # input: raw reads from step 00
trim_dir="$out/fastq_trimmed"    # output: cleaned reads (used by step 02)
rep_dir="$out/fastp_reports"     # fastp's html/json reports and logs
threads=8

mkdir -p "$trim_dir" "$rep_dir"

# Loop over each sample and clean up its two read files.
for sample in CHIRP-D1_S10_L002 \
    CHIRP-D2_S11_L002 \
    CHIRP-D3_S12_L002 \
    CHIRP-V1_S7_L002 \
    CHIRP-V2_S8_L002 \
    CHIRP-V3_S9_L002 \
    Input-D1_S4_L002 \
    Input-D2_S5_L002 \
    Input-D3_S6_L002 \
    Input-v1_S1_L002 \
    Input-v2_S2_L002 \
    Input-v3_S3_L002; do
    echo "=== trimming $sample ==="
    fastp \
        --in1  "$fastq_dir/${sample}_R1_001.fastq.gz" \
        --in2  "$fastq_dir/${sample}_R2_001.fastq.gz" \
        --out1 "$trim_dir/${sample}_R1.fastq.gz" \
        --out2 "$trim_dir/${sample}_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --length_required 20 \
        --thread "$threads" \
        --html "$rep_dir/${sample}_fastp.html" \
        --json "$rep_dir/${sample}_fastp.json" \
        2> "$rep_dir/${sample}_fastp.log"
    echo "Done: $sample"
done

echo "Cleaned reads are in: $trim_dir"
