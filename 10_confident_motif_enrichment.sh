#!/usr/bin/env bash
# ============================================================
# STEP 10 - Motif enrichment of confident D peaks
#
# Purpose:
#   Identify enriched sequence motifs within reproducible
#   LNC-D ChIRP binding regions.
#
# Input:
#   LNCD_CHIRP-D_confident_peaks.bed
#
# Genome:
#   hg38
#
# Tool:
#   HOMER findMotifsGenome.pl
#
# Output:
#   Known motif enrichment
#   De novo motif enrichment
# ============================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=10_confident_D_motif
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

# ============================================================
# Activate software
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ============================================================
# Paths
# ============================================================

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

peak_dir="$out/peaks/confident"

D_peaks="$peak_dir/LNCD_CHIRP-D_confident_peaks.bed"

step10_dir="$out/motif_enrichment/10_confident_D"

mkdir -p "$step10_dir"

# ============================================================
# Check HOMER
# ============================================================

if ! command -v findMotifsGenome.pl >/dev/null 2>&1; then

    echo "ERROR: findMotifsGenome.pl not found."
    echo "Check HOMER installation in CHIRP-SEQ."

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
echo "STEP 10 - Motif enrichment"
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
# Run HOMER motif enrichment
#
# -size given
#   Use the actual confident peak regions.
#
# -len 8,10,12
#   Search de novo motifs of these lengths.
#
# -p 8
#   Use 8 CPU threads.
#
# HOMER automatically generates background regions matched
# for sequence properties.
# ============================================================

echo ""
echo "========================================"
echo "Running HOMER findMotifsGenome.pl"
echo "========================================"

findMotifsGenome.pl \
    "$D_peaks" \
    hg38 \
    "$step10_dir" \
    -size given \
    -len 8,10,12 \
    -p 8


# ============================================================
# Check output
# ============================================================

echo ""
echo "========================================"
echo "Motif analysis output"
echo "========================================"

echo ""

ls -lh "$step10_dir"


# ============================================================
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
echo "Important results:"
echo "$step10_dir/knownResults.txt"
echo "$step10_dir/homerResults.html"
echo "$step10_dir/knownResults.html"

echo ""
echo "Known motifs:"
echo "$step10_dir/knownResults.txt"

echo ""
echo "De novo motifs:"
echo "$step10_dir/homerResults.html"