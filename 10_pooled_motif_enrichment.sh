#!/usr/bin/env bash
# ============================================================
# STEP 10 - Motif enrichment in D-specific ChIRP-seq peaks
#
# D = treatment
# V = control
#
# Foreground:
#   D-specific pooled peaks from Step 07
#
# Background:
#   V pooled peaks from Step 07
#
# Question:
#   Are particular DNA sequence motifs enriched in
#   D-specific LNC-D binding regions relative to
#   control-associated LNC-D binding regions?
#
# Genome:
#   hg38
#
# IMPORTANT:
#   There are only 6 D-specific peaks.
#   Therefore this analysis should be considered exploratory.
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=10_pooled_motif_enrichment
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu
# ============================================================
# Activate software
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# Paths
# ============================================================

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

step7_dir="$out/peaks/07_D_specific_pooled"

step10_dir="$out/motif/10_D_specific_motif"

mkdir -p "$step10_dir"

# Input files
# ============================================================

D_specific="$step7_dir/D_specific_peaks_ID.bed"

V_background="$step7_dir/V_pooled_peaks.bed"

# Check input files
# ============================================================

for f in "$D_specific" "$V_background"; do

    if [[ ! -s "$f" ]]; then

        echo "ERROR: Cannot find or empty file:"
        echo "$f"

        exit 1
    fi

done

echo ""
echo "========================================"
echo "STEP 10 - Motif enrichment"
echo "========================================"

echo ""
echo "Foreground:"
echo "$D_specific"

echo ""
echo "D-specific peaks:"
wc -l "$D_specific"

echo ""
echo "Background:"
echo "$V_background"

echo ""
echo "V background peaks:"
wc -l "$V_background"

# Check HOMER
# ============================================================

if ! command -v findMotifsGenome.pl >/dev/null 2>&1; then

    echo "ERROR: findMotifsGenome.pl not found."
    echo "HOMER is not available in CHIRP-SEQ."

    exit 1

fi

# Check HOMER hg38
# ============================================================

homer_hg38="$CONDA_PREFIX/share/homer/data/genomes/hg38"

if [[ ! -d "$homer_hg38" ]]; then

    echo "ERROR: HOMER hg38 genome is not installed."
    echo "Run Step 08 hg38 installation first."

    exit 1

fi


echo ""
echo "HOMER hg38 is ready."


# ============================================================
# Run motif enrichment
#
# -bg:
#   use V peaks as the background
#
# -size given:
#   use the actual peak sizes instead of resizing all peaks
#
# -mask:
#   mask repetitive sequence
#
# ============================================================

echo ""
echo "========================================"
echo "Running HOMER motif enrichment"
echo "========================================"


findMotifsGenome.pl \
    "$D_specific" \
    hg38 \
    "$step10_dir" \
    -bg "$V_background" \
    -size given \
    -mask \
    -p 8



# Check output
# ============================================================

if [[ ! -s "$step10_dir/knownResults.txt" ]]; then

    echo ""
    echo "WARNING:"
    echo "knownResults.txt was not generated or is empty."
    echo ""
    echo "Because only 6 D-specific peaks are available,"
    echo "motif enrichment may have insufficient power."

else

    echo ""
    echo "Known motif enrichment completed."

fi

# Print top known motifs
# ============================================================

if [[ -s "$step10_dir/knownResults.txt" ]]; then

    echo ""
    echo "========================================"
    echo "Top known motifs"
    echo "========================================"

    head -11 "$step10_dir/knownResults.txt"

fi

# Check de novo motifs
# ============================================================

if [[ -s "$step10_dir/homerResults.html" ]]; then

    echo ""
    echo "De novo motif results generated."

fi

# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 10 COMPLETE"
echo "========================================"

echo ""
echo "Output directory:"
echo "$step10_dir"

echo ""
echo "Main known-motif result:"
echo "$step10_dir/knownResults.txt"

echo ""
echo "Known motif HTML:"
echo "$step10_dir/knownResults.html"

echo ""
echo "De novo motif HTML:"
echo "$step10_dir/homerResults.html"