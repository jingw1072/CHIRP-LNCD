#!/usr/bin/env bash
# ============================================================
# STEP 08 - confident Peak annotation
#
# Purpose:
#   Annotate reproducible LNC-D ChIRP peaks with:
#
#   - genomic annotation
#   - nearest gene
#   - nearest TSS
#   - distance to TSS
#
# Input:
#   confident_D.bed generated from replicate-supported
#   D peaks (>=2/3 replicates)
#
# Genome:
#   hg38
#
# Tool:
#   HOMER annotatePeaks.pl
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=08_confidentpeak_annotation
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiw169@ucsd.edu
# Activate software
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths -----------------------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

peak_dir="$out/peaks/confident"

step8_dir="$out/confident_peak_annotation/08_confident_D"

mkdir -p "$step8_dir"


# ============================================================
# Input peak file
# ============================================================
D_peaks="$peak_dir/LNCD_CHIRP-D_confident_peaks.bed"

# ============================================================
# Check required program
# ============================================================

if ! command -v annotatePeaks.pl >/dev/null 2>&1; then

    echo "ERROR: HOMER annotatePeaks.pl not found."
    echo "Check whether HOMER is installed in CHIRP-SEQ."

    exit 1

fi

# ============================================================
# Check input
# ============================================================

if [[ ! -s "$D_peaks" ]]; then

    echo "ERROR: Cannot find confident D peak file:"
    echo "$D_peaks"

    exit 1

fi

# ============================================================
# Report input
# ============================================================

echo ""
echo "========================================"
echo "STEP 08 - Peak annotation"
echo "========================================"

echo ""
echo "Input:"
echo "$D_peaks"

echo ""
echo "Number of confident D peaks:"
wc -l "$D_peaks"

echo ""
echo "Genome:"
echo "hg38"


# ============================================================
# A) HOMER peak annotation
#
# Output includes:
#
#   Annotation
#   Detailed Annotation
#   Distance to TSS
#   Nearest PromoterID
#   Entrez ID
#   Nearest Unigene
#   Nearest Refseq
#   Nearest Ensembl
#   Gene Name
#   Gene Alias
#   Gene Description
#   Gene Type
# ============================================================

echo ""
echo "========================================"
echo "Running HOMER annotatePeaks.pl"
echo "========================================"

annotatePeaks.pl \
    "$D_peaks" \
    hg38 \
    > "$step8_dir/confident_D_HOMER_annotation.tsv"


# ============================================================
# B) Basic QC
# ============================================================

echo ""
echo "========================================"
echo "Annotation QC"
echo "========================================"

echo ""
echo "Input peaks:"
wc -l "$D_peaks"

echo ""
echo "Annotation table rows:"
echo "$(($(wc -l < "$step8_dir/confident_D_HOMER_annotation.tsv") - 1))"


# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 08 COMPLETE"
echo "========================================"

echo ""
echo "Output directory:"
echo "$step8_dir"

echo ""
echo "Main annotation file:"
echo "$step8_dir/confident_D_HOMER_annotation.tsv"

echo ""
echo "Next step:"
echo "STEP 09 - Genomic distribution"