#!/usr/bin/env bash
# ============================================================
# STEP 07B
# Run DESeq2 differential binding analysis
# Comparison:
#   D1/D2/D3 vs V1/V2/V3
#
# Input from Step 07A:
#   DV_union_DESeq2_counts.tsv
#   DV_union_counts.tsv
#
# R script:
#   07B_DESeq2.R
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=07b_DESeq2_D_vs_V
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiw169@ucsd.edu
# ============================================================
# Activate software
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-DESEQ2

# ---- Paths -----------------------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

script_dir="$out/scripts"

R_script="/tscc/nfs/home/jiw169/github/CHIRP-LNCD/07B_DESeq2.R"


# ============================================================
# Check R
# ============================================================

echo ""
echo "========================================"
echo "Software check"
echo "========================================"

echo "Conda environment:"
echo "$CONDA_DEFAULT_ENV"

echo ""
echo "Rscript:"
which Rscript

echo ""
Rscript --version

# ============================================================
# Check R script
# ============================================================

if [[ ! -s "$R_script" ]]; then

    echo "ERROR: Cannot find R script:"
    echo "$R_script"

    exit 1

fi

# ============================================================
# Check required R packages
# ============================================================

echo ""
echo "Checking required R packages..."

Rscript -e '
packages <- c("DESeq2", "ggplot2", "pheatmap")

missing <- packages[
    !sapply(packages, requireNamespace, quietly = TRUE)
]

if (length(missing) > 0) {
    stop(
        paste(
            "Missing R packages:",
            paste(missing, collapse = ", ")
        )
    )
}

cat("All required R packages are available.\n")
'

# ============================================================
# Run DESeq2
# ============================================================

echo ""
echo "========================================"
echo "STEP 07B - DESeq2 D vs V"
echo "========================================"

echo ""
echo "Start time:"
date

echo ""
echo "R script:"
echo "$R_script"

echo ""

Rscript "$R_script"


# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 07B COMPLETE"
echo "========================================"

echo ""
echo "End time:"
date