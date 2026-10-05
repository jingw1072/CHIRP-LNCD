#!/usr/bin/env bash
# ============================================================
# STEP 07C - Replicate QC
#
# Purpose:
#   Evaluate D1/D2/D3 and V1/V2/V3 before considering
#   exclusion of any biological replicate.
#
# Analyses:
#   A. BAM statistics
#   B. FRiP
#   C. Genome-wide replicate correlation
#   D. PCA
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=07C_replicate_QC
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiw169@ucsd.edu

# Activate software
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths -----------------------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

bam_dir="$out/bam"
peak_dir="$out/peaks/replicates"

qc_dir="$out/07C_replicate_QC"

mkdir -p "$qc_dir"

# ============================================================
# Define samples
# ============================================================

samples=(
    "D1"
    "D2"
    "D3"
    "V1"
    "V2"
    "V3"
)

bams=(
    "$bam_dir/CHIRP-D1_S10_L002_filt.bam"
    "$bam_dir/CHIRP-D2_S11_L002_filt.bam"
    "$bam_dir/CHIRP-D3_S12_L002_filt.bam"
    "$bam_dir/CHIRP-V1_S7_L002_filt.bam"
    "$bam_dir/CHIRP-V2_S8_L002_filt.bam"
    "$bam_dir/CHIRP-V3_S9_L002_filt.bam"
)


# ============================================================
# A. BAM statistics
# ============================================================

echo ""
echo "========================================"
echo "A. BAM statistics"
echo "========================================"

printf "Sample\tMapped_reads\tProperly_paired\n" \
    > "$qc_dir/BAM_QC_summary.tsv"

for i in "${!samples[@]}"
do

    sample="${samples[$i]}"
    bam="${bams[$i]}"

    echo "Processing $sample"

    if [[ ! -s "$bam" ]]; then
        echo "ERROR: BAM not found:"
        echo "$bam"
        exit 1
    fi

    samtools flagstat "$bam" \
        > "$qc_dir/${sample}_flagstat.txt"

    mapped=$(samtools view -c -F 4 "$bam")

    proper=$(samtools view -c -f 2 "$bam")

    printf "%s\t%s\t%s\n" \
        "$sample" \
        "$mapped" \
        "$proper" \
        >> "$qc_dir/BAM_QC_summary.tsv"

done


# ============================================================
# B. FRiP
#
# FRiP = reads overlapping that replicate's own confident
# MACS2 peak set / mapped reads
#
# IMPORTANT:
# Change peak filenames below if your Step 05A names differ.
# ============================================================

echo ""
echo "========================================"
echo "B. FRiP"
echo "========================================"

printf "Sample\tMapped_reads\tReads_in_peaks\tFRiP\n" \
    > "$qc_dir/FRiP_summary.tsv"


for i in "${!samples[@]}"
do

    sample="${samples[$i]}"
    bam="${bams[$i]}"

    # Find replicate peak file
    case "$sample" in

        D1)
            peak="$peak_dir/LNCD_CHIRP-D1_FC5_noBL.narrowPeak"
            ;;

        D2)
            peak="$peak_dir/LNCD_CHIRP-D2_FC5_noBL.narrowPeak"
            ;;

        D3)
            peak="$peak_dir/LNCD_CHIRP-D3_FC5_noBL.narrowPeak"
            ;;

        V1)
            peak="$peak_dir/LNCD_CHIRP-V1_FC5_noBL.narrowPeak"
            ;;

        V2)
            peak="$peak_dir/LNCD_CHIRP-V2_FC5_noBL.narrowPeak"
            ;;

        V3)
            peak="$peak_dir/LNCD_CHIRP-V3_FC5_noBL.narrowPeak"
            ;;

    esac


    if [[ ! -s "$peak" ]]; then
        echo "ERROR: Peak file not found:"
        echo "$peak"
        exit 1
    fi


    mapped=$(samtools view -c -F 4 "$bam")

    in_peaks=$(bedtools intersect \
        -u \
        -abam "$bam" \
        -b "$peak" \
        | samtools view -c)


    frip=$(awk \
        -v x="$in_peaks" \
        -v y="$mapped" \
        'BEGIN {
            if (y > 0)
                printf "%.6f", x/y;
            else
                print "NA"
        }')


    printf "%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$mapped" \
        "$in_peaks" \
        "$frip" \
        >> "$qc_dir/FRiP_summary.tsv"

done


# ============================================================
# C. Genome-wide correlation
#
# multiBamSummary bins:
#   counts signal across genome-wide bins for all replicates.
#
# This avoids judging replicates only by MACS2 peak number.
# ============================================================

echo ""
echo "========================================"
echo "C. Genome-wide signal correlation"
echo "========================================"


multiBamSummary bins \
    --bamfiles \
        "${bams[@]}" \
    --labels \
        D1 D2 D3 V1 V2 V3 \
    --binSize 10000 \
    --numberOfProcessors 8 \
    -out "$qc_dir/replicate_signal.npz" \
    --outRawCounts "$qc_dir/replicate_signal_rawCounts.tsv"


# Pearson correlation

plotCorrelation \
    -in "$qc_dir/replicate_signal.npz" \
    --corMethod pearson \
    --skipZeros \
    --whatToPlot heatmap \
    --plotNumbers \
    --plotFile "$qc_dir/replicate_Pearson_heatmap.pdf" \
    --outFileCorMatrix "$qc_dir/replicate_Pearson_matrix.tsv"


# Spearman correlation

plotCorrelation \
    -in "$qc_dir/replicate_signal.npz" \
    --corMethod spearman \
    --skipZeros \
    --whatToPlot heatmap \
    --plotNumbers \
    --plotFile "$qc_dir/replicate_Spearman_heatmap.pdf" \
    --outFileCorMatrix "$qc_dir/replicate_Spearman_matrix.tsv"


# ============================================================
# D. PCA
# ============================================================

echo ""
echo "========================================"
echo "D. PCA"
echo "========================================"


plotPCA \
    -in "$qc_dir/replicate_signal.npz" \
    --plotFile "$qc_dir/replicate_PCA.pdf" \
    --outFileNameData "$qc_dir/replicate_PCA_data.tsv"


# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 07C COMPLETE"
echo "========================================"

echo ""
echo "Output directory:"
echo "$qc_dir"

echo ""
echo "Important outputs:"
echo "  BAM_QC_summary.tsv"
echo "  FRiP_summary.tsv"
echo "  replicate_Pearson_heatmap.pdf"
echo "  replicate_Spearman_heatmap.pdf"
echo "  replicate_PCA.pdf"