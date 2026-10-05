#!/usr/bin/env bash
# ============================================================
# STEP 09 - Genomic distribution of D-specific ChIRP-seq peaks
#
# D = treatment
# V = control
#
# Input:
#   HOMER annotation generated in Step 08:
#
#   D_specific_peaks_annotation.txt
#
# Purpose:
#   Classify D-specific pooled peaks into genomic categories:
#
#     Promoter
#     5' UTR
#     3' UTR
#     Exon
#     Intron
#     Intergenic
#     Other
#
# Outputs:
#   1) peak-level genomic category table
#   2) category counts
#   3) category percentages
#
# NOTE:
#   These are D-specific pooled peaks, not statistically
#   significant differential-binding peaks.
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=09_pooled_Dpeaks_Genomic_distribution
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
# ---- Activate software ------------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths -----------------------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

step8_dir="$out/annotation/08_D_specific_peaks"

step9_dir="$out/annotation/09_D_specific_genomic_distribution"

mkdir -p "$step9_dir"

# ---- Input -----------------------------------------------------------------

annotation="$step8_dir/D_specific_peaks_annotation.txt"

# A) Check input
# ============================================================

if [[ ! -s "$annotation" ]]; then

    echo "ERROR: Cannot find Step 08 annotation file:"
    echo "$annotation"
    exit 1

fi

echo ""
echo "========================================"
echo "STEP 09 - Genomic distribution"
echo "========================================"

echo ""
echo "Input:"
echo "$annotation"

echo ""
echo "Number of annotated peaks:"
echo "$(($(wc -l < "$annotation") - 1))"

# ============================================================
# B) Show HOMER header
#
# This is useful because HOMER versions can have slightly
# different output columns.
#
# The broad genomic annotation is normally column 8.
# ============================================================

echo ""
echo "HOMER annotation header:"
head -1 "$annotation"

# ============================================================
# C) Assign each peak to a simplified genomic category
#
# HOMER annotation examples:
#
#   promoter-TSS (...)
#   exon (...)
#   intron (...)
#   Intergenic
#   5' UTR
#   3' UTR
#
# Output:
#
# PeakID   HOMER_annotation   Category
# ============================================================

awk -F'\t' '
BEGIN {
    OFS="\t"
    print "PeakID","HOMER_annotation","Category"
}

NR > 1 {

    peak=$1
    ann=$8

    ann_lower=tolower(ann)

    if (ann_lower ~ /promoter/) {
        category="Promoter"
    }
    else if (ann_lower ~ /5'\'' utr/) {
        category="5UTR"
    }
    else if (ann_lower ~ /3'\'' utr/) {
        category="3UTR"
    }
    else if (ann_lower ~ /exon/) {
        category="Exon"
    }
    else if (ann_lower ~ /intron/) {
        category="Intron"
    }
    else if (ann_lower ~ /intergenic/) {
        category="Intergenic"
    }
    else {
        category="Other"
    }

    print peak,ann,category
}
' "$annotation" \
> "$step9_dir/D_specific_peak_categories.tsv"


# ============================================================
# D) Count peaks in each category
# ============================================================

awk -F'\t' '
BEGIN {
    OFS="\t"
}

NR > 1 {
    count[$3]++
    total++
}

END {

    print "Category","Count","Percentage"

    order[1]="Promoter"
    order[2]="5UTR"
    order[3]="3UTR"
    order[4]="Exon"
    order[5]="Intron"
    order[6]="Intergenic"
    order[7]="Other"

    for (i=1; i<=7; i++) {

        category=order[i]

        n=count[category]+0

        if (total > 0)
            pct=(n/total)*100
        else
            pct=0

        printf "%s\t%d\t%.2f\n", category,n,pct
    }
}
' "$step9_dir/D_specific_peak_categories.tsv" \
> "$step9_dir/D_specific_genomic_distribution.tsv"

# E) Remove zero-count categories
#
# This gives a compact table for plotting.
# ============================================================

awk -F'\t' '
BEGIN {
    OFS="\t"
}

NR==1 {
    print
    next
}

$2 > 0 {
    print
}
' "$step9_dir/D_specific_genomic_distribution.tsv" \
> "$step9_dir/D_specific_genomic_distribution_nonzero.tsv"

# F) Print results
# ============================================================

echo ""
echo "========================================"
echo "Genomic distribution"
echo "========================================"

echo ""

cat "$step9_dir/D_specific_genomic_distribution.tsv"

# G) Check total
# ============================================================

total=$(awk -F'\t' 'NR>1 {sum += $2} END {print sum+0}' \
    "$step9_dir/D_specific_genomic_distribution.tsv")


echo ""
echo "Total classified peaks: $total"

# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 09 COMPLETE"
echo "========================================"

echo ""
echo "Output directory:"
echo "$step9_dir"

echo ""
echo "Peak-level categories:"
echo "$step9_dir/D_specific_peak_categories.tsv"

echo ""
echo "Genomic distribution:"
echo "$step9_dir/D_specific_genomic_distribution.tsv"

echo ""
echo "Non-zero categories:"
echo "$step9_dir/D_specific_genomic_distribution_nonzero.tsv"