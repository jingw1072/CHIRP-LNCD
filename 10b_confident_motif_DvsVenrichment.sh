#!/usr/bin/env bash
# ============================================================
# STEP 10b - D vs V motif enrichment
#
# Purpose:
#   Identify sequence motifs enriched in reproducible
#   LNC-D ChIRP binding regions relative to reproducible
#   V-control binding regions.
#
# Target:
#   confident D peaks
#
# Background/control:
#   confident V peaks
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
#SBATCH --job-name=10b_confident_D_vs_V_motif
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
V_peaks="$peak_dir/LNCD_CHIRP-V_confident_peaks.bed"

step10b_dir="$out/motif_enrichment/10B_D_vs_V"

mkdir -p "$step10b_dir"


# ============================================================
# Check HOMER
# ============================================================

if ! command -v findMotifsGenome.pl >/dev/null 2>&1; then

    echo "ERROR: findMotifsGenome.pl not found."
    echo "Check HOMER installation in CHIRP-SEQ."

    exit 1

fi

# ============================================================
# Check input files
# ============================================================

if [[ ! -s "$D_peaks" ]]; then

    echo "ERROR: Cannot find confident D peaks:"
    echo "$D_peaks"

    exit 1

fi

if [[ ! -s "$V_peaks" ]]; then

    echo "ERROR: Cannot find confident V peaks:"
    echo "$V_peaks"

    exit 1

fi

# ============================================================
# Input QC
# ============================================================

echo ""
echo "========================================"
echo "STEP 10B - D vs V motif enrichment"
echo "========================================"

echo ""
echo "Target:"
echo "$D_peaks"

echo ""
echo "Number of confident D peaks:"
wc -l "$D_peaks"

echo ""
echo "Background/control:"
echo "$V_peaks"

echo ""
echo "Number of confident V peaks:"
wc -l "$V_peaks"

echo ""
echo "Genome:"
echo "hg38"

# ============================================================
# Peak-width QC
#
# Because -size given is used, check whether D and V peak
# widths are substantially different.
# ============================================================

echo ""
echo "========================================"
echo "Peak width summary"
echo "========================================"

echo ""
echo "D peak widths:"
awk '
{
    width=$3-$2
    sum+=width
    n++
}
END {
    if(n>0)
        printf "N=%d\tMean=%.2f bp\n",n,sum/n
}' "$D_peaks"


echo "V peak widths:"
awk '
{
    width=$3-$2
    sum+=width
    n++
}
END {
    if(n>0)
        printf "N=%d\tMean=%.2f bp\n",n,sum/n
}' "$V_peaks"


# ============================================================
# Run HOMER
#
# Target:
#   confident D peaks
#
# Background:
#   confident V peaks
#
# -bg:
#   explicitly tells HOMER to use V peaks as background
#
# -size given:
#   use the original peak widths
#
# -len:
#   de novo motif lengths
#
# -p:
#   threads
# ============================================================

echo ""
echo "========================================"
echo "Running HOMER: confident D vs confident V"
echo "========================================"

findMotifsGenome.pl \
    "$D_peaks" \
    hg38 \
    "$step10b_dir" \
    -bg "$V_peaks" \
    -size given \
    -len 8,10,12 \
    -p 8


# ============================================================
# Check outputs
# ============================================================

echo ""
echo "========================================"
echo "STEP 10B results"
echo "========================================"

ls -lh "$step10b_dir"


# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 10B COMPLETE"
echo "========================================"

echo ""
echo "Target:"
echo "confident D peaks"

echo ""
echo "Background:"
echo "confident V peaks"

echo ""
echo "Output directory:"
echo "$step10b_dir"

echo ""
echo "Important outputs:"
echo "$step10b_dir/knownResults.txt"
echo "$step10b_dir/knownResults.html"
echo "$step10b_dir/homerResults.html"