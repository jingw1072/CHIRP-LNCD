#!/usr/bin/env bash
# ============================================================
# STEP 08 - Annotate D-specific pooled ChIRP-seq peaks
#
# D = treatment
# V = control
#
# Input:
#   D-specific pooled peaks generated in Step 07
#
#     D_specific_peaks_ID.bed
#
# Purpose:
#   Annotate D-specific LNC-D binding regions with:
#
#     - genomic location
#     - promoter / exon / intron / intergenic annotation
#     - nearest transcription start site (TSS)
#     - distance to TSS
#     - nearest gene
#     - gene symbol
#
# Genome:
#   hg38 / GRCh38
#
# Main output:
#   D_specific_peaks_annotation.txt
#
# NOTE:
#   These are annotations of D-specific POOLED peaks.
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=Annotate_pooled_Dspecific_peaks
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

# ---- Activate software ------------------------------------------------------

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ
# Check HOMER
# ============================================================

if ! command -v annotatePeaks.pl >/dev/null 2>&1; then
    echo "ERROR: annotatePeaks.pl not found."
    echo "HOMER is not installed in the current environment."
    exit 1
fi
# Check/install HOMER hg38
# ============================================================

homer_config="$CONDA_PREFIX/share/homer/configureHomer.pl"
homer_hg38="$CONDA_PREFIX/share/homer/data/genomes/hg38"

echo "Checking HOMER hg38 genome..."

if [[ ! -d "$homer_hg38" ]]; then
    echo "HOMER hg38 is not installed. Installing..."
    perl "$homer_config" -install hg38
else
    echo "HOMER hg38 is already installed."
fi

# Verify
if [[ ! -d "$homer_hg38" ]]; then
    echo "ERROR: HOMER hg38 genome directory not found."
    exit 1
fi

echo "HOMER hg38 is ready."

# Paths
# ============================================================

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

step7_dir="$out/peaks/07_D_specific_pooled"
step8_dir="$out/annotation/08_D_specific_peaks"

mkdir -p "$step8_dir"

D_specific="$step7_dir/D_specific_peaks_ID.bed"

# Check input
# ============================================================

if [[ ! -s "$D_specific" ]]; then
    echo "ERROR: Cannot find:"
    echo "$D_specific"
    exit 1
fi

# Annotate peaks
# ============================================================

echo "Annotating D-specific peaks..."

annotatePeaks.pl \
    "$D_specific" \
    hg38 \
    > "$step8_dir/D_specific_peaks_annotation.txt"


# Check output
if [[ ! -s "$step8_dir/D_specific_peaks_annotation.txt" ]]; then
    echo "ERROR: Annotation output is empty."
    exit 1
fi


echo "STEP 08 complete."
echo "Output:"
echo "$step8_dir/D_specific_peaks_annotation.txt"