#!/usr/bin/env bash
# ============================================================
# STEP 07 - Identify D-specific pooled ChIRP-seq peaks
#
# D = treatment
# V = control
#
# Input:
#   D pooled peaks:
#     LNCD_CHIRP-D_FC5_noBL.narrowPeak
#
#   V pooled peaks:
#     LNCD_CHIRP-V_FC5_noBL.narrowPeak
#
# Both peak sets were generated in Step 4 using:
#   MACS2 vs matched Input
#   q < 0.01
#   FC >= 5
#   blacklist removal
#
# Definition:
#   A D peak is considered shared with V if it overlaps
#   a V peak by >= 100 bp.
#
#   D-specific peak =
#   D peak without >=100 bp overlap with any V peak.
#
# IMPORTANT:
#   These are D-specific POOLED peaks.
#   They are not statistically significant differential peaks.
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_specific_peaks
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu
# ---- Activate software ------------------------------------------------------

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths -----------------------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

peak_dir="$out/peaks"

step7_dir="$peak_dir/07_D_specific_pooled"

mkdir -p "$step7_dir"


# ---- Parameters -------------------------------------------------------------

min_overlap=100

# ---- Input peak files -------------------------------------------------------

D_np="$peak_dir/LNCD_CHIRP-D_FC5_noBL.narrowPeak"
V_np="$peak_dir/LNCD_CHIRP-V_FC5_noBL.narrowPeak"

# ============================================================
# Check input files
# ============================================================

for f in "$D_np" "$V_np"; do

    if [[ ! -s "$f" ]]; then
        echo "ERROR: Cannot find or empty file:"
        echo "$f"
        exit 1
    fi

done


if ! command -v bedtools >/dev/null 2>&1; then
    echo "ERROR: bedtools not found."
    exit 1
fi


echo ""
echo "========================================"
echo "STEP 07 - D-specific pooled peaks"
echo "========================================"

echo ""
echo "D = treatment"
echo "V = control"
echo "Minimum shared overlap = ${min_overlap} bp"

# ============================================================
# A) Convert narrowPeak to BED4
# ============================================================

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' "$D_np" \
    | sort -k1,1 -k2,2n \
    > "$step7_dir/D_pooled_peaks.bed"

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' "$V_np" \
    | sort -k1,1 -k2,2n \
    > "$step7_dir/V_pooled_peaks.bed"

D_bed="$step7_dir/D_pooled_peaks.bed"
V_bed="$step7_dir/V_pooled_peaks.bed"

echo ""
echo "D total peaks: $(wc -l < "$D_bed")"
echo "V total peaks: $(wc -l < "$V_bed")"

# ============================================================
# B) Find D peaks that overlap V by >=100 bp
#
# Output columns:
#
# D peak: 1-4
# V peak: 5-8
# overlap: 9
# ============================================================

echo ""
echo "Finding D peaks shared with V..."

bedtools intersect \
    -a "$D_bed" \
    -b "$V_bed" \
    -wo \
    | awk -v m="$min_overlap" '
        BEGIN{OFS="\t"}
        $9 >= m
    ' \
    > "$step7_dir/D_vs_V_overlap_pairs.tsv"

# ============================================================
# C) Get unique D peaks that overlap V
# ============================================================

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' \
    "$step7_dir/D_vs_V_overlap_pairs.tsv" \
    | sort -k1,1 -k2,2n -u \
    > "$step7_dir/D_shared_with_V.bed"

# ============================================================
# D) Extract D-specific peaks
#
# Remove D peaks that were identified as shared with V.
# ============================================================

cut -f4 "$step7_dir/D_shared_with_V.bed" \
    | sort -u \
    > "$step7_dir/D_shared_peak_names.txt"


awk '
    NR==FNR {
        shared[$1]=1
        next
    }

    !($4 in shared)
' \
    "$step7_dir/D_shared_peak_names.txt" \
    "$D_bed" \
    > "$step7_dir/D_specific_peaks.bed"


# ============================================================
# E) Add simple peak IDs
#
# This will be the main file for downstream annotation,
# motif enrichment, etc.
# ============================================================

awk 'BEGIN{OFS="\t"}
{
    print $1,$2,$3,"D_specific_"NR
}' "$step7_dir/D_specific_peaks.bed" \
    > "$step7_dir/D_specific_peaks_ID.bed"

# ============================================================
# F) Summary
# ============================================================

n_D=$(wc -l < "$D_bed")
n_V=$(wc -l < "$V_bed")
n_shared=$(wc -l < "$step7_dir/D_shared_with_V.bed")
n_D_specific=$(wc -l < "$step7_dir/D_specific_peaks.bed")

echo ""
echo "========================================"
echo "STEP 07 SUMMARY"
echo "========================================"

echo ""
echo "D total peaks:             $n_D"
echo "V total peaks:             $n_V"
echo "D peaks shared with V:     $n_shared"
echo "D-specific pooled peaks:   $n_D_specific"


# Save summary

{
    printf "category\tcount\n"
    printf "D_total\t%s\n" "$n_D"
    printf "V_total\t%s\n" "$n_V"
    printf "D_shared_with_V\t%s\n" "$n_shared"
    printf "D_specific\t%s\n" "$n_D_specific"

} > "$step7_dir/D_specific_summary.tsv"


# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 07 COMPLETE"
echo "========================================"

echo ""
echo "Main output:"
echo "$step7_dir/D_specific_peaks_ID.bed"

echo ""
echo "Summary:"
echo "$step7_dir/D_specific_summary.tsv"