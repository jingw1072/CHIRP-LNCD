#!/usr/bin/env bash
# =============================================================================
# STEP 5b - Identify confident ChIRP-seq peaks
#
# Input:
#   Individual replicate peak files generated in Step 5a.
#
# Each replicate peak set has already been:
#   1) called against its matched Input using MACS2 (q < 0.01)
#   2) filtered for fold enrichment (FC >= 5)
#   3) filtered to remove ENCODE blacklist regions
# Definition used here:
#
#   Confident peak = a genomic region supported by at least
#   2 of 3 biological replicates.
#
#   Peaks from two replicates must overlap by >= 100 bp.

# Final output:
#   LNCD_CHIRP-D_confident_peaks.bed
#   LNCD_CHIRP-V_confident_peaks.bed
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 05b_confident_peaks.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_confident
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ
# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"
rep_dir="$out/peaks/replicates"
conf_dir="$out/peaks/confident"

fc_min=5
min_overlap=100

mkdir -p "$conf_dir"

# ============================================================
# Function: identify >=2/3 replicate-supported peaks
# ============================================================

find_confident_peaks () {

    condition=$1

    echo ""
    echo "========================================"
    echo "Condition: $condition"
    echo "========================================"

    # --------------------------------------------------------
    # Individual replicate peaks from Step 5A
    # --------------------------------------------------------

    rep1="$rep_dir/LNCD_CHIRP-${condition}1_FC${fc_min}_noBL.narrowPeak"
    rep2="$rep_dir/LNCD_CHIRP-${condition}2_FC${fc_min}_noBL.narrowPeak"
    rep3="$rep_dir/LNCD_CHIRP-${condition}3_FC${fc_min}_noBL.narrowPeak"

    # --------------------------------------------------------
    # Check that all three input files exist
    # --------------------------------------------------------

    for f in "$rep1" "$rep2" "$rep3"; do

        if [[ ! -f "$f" ]]; then
            echo "ERROR: Cannot find:"
            echo "$f"
            exit 1
        fi

    done

    # --------------------------------------------------------
    # Convert narrowPeak files to BED3
    #
    # Keep only:
    #   chromosome
    #   start
    #   end
    #
    # These files are used only for replicate overlap analysis.
    # --------------------------------------------------------

    for rep in 1 2 3; do

        input="$rep_dir/LNCD_CHIRP-${condition}${rep}_FC${fc_min}_noBL.narrowPeak"

        awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$input" \
            | sort -k1,1 -k2,2n \
            > "$conf_dir/${condition}${rep}.bed"

    done

    # ========================================================
    # Find pairwise overlaps
    #
    # Because both A and B are BED3:
    #
    #   A: chr start end        columns 1-3
    #   B: chr start end        columns 4-6
    #   overlap length          column 7
    #
    # Only keep pairs overlapping by >= min_overlap bp.
    #
    # For each pair, output the ACTUAL shared interval.
    # ========================================================

    # --------------------------------------------------------
    # Rep1 vs Rep2
    # --------------------------------------------------------

    bedtools intersect \
        -a "$conf_dir/${condition}1.bed" \
        -b "$conf_dir/${condition}2.bed" \
        -wo \
        | awk -v m="$min_overlap" '
        BEGIN{OFS="\t"}
        $7 >= m {
            start=($2>$5 ? $2 : $5)
            end=($3<$6 ? $3 : $6)
            print $1,start,end
        }' \
        > "$conf_dir/${condition}_rep1_rep2_overlap.bed"

    # --------------------------------------------------------
    # Rep1 vs Rep3
    # --------------------------------------------------------

    bedtools intersect \
        -a "$conf_dir/${condition}1.bed" \
        -b "$conf_dir/${condition}3.bed" \
        -wo \
        | awk -v m="$min_overlap" '
        BEGIN{OFS="\t"}
        $7 >= m {
            start=($2>$5 ? $2 : $5)
            end=($3<$6 ? $3 : $6)
            print $1,start,end
        }' \
        > "$conf_dir/${condition}_rep1_rep3_overlap.bed"

    # --------------------------------------------------------
    # Rep2 vs Rep3
    # --------------------------------------------------------

    bedtools intersect \
        -a "$conf_dir/${condition}2.bed" \
        -b "$conf_dir/${condition}3.bed" \
        -wo \
        | awk -v m="$min_overlap" '
        BEGIN{OFS="\t"}
        $7 >= m {
            start=($2>$5 ? $2 : $5)
            end=($3<$6 ? $3 : $6)
            print $1,start,end
        }' \
        > "$conf_dir/${condition}_rep2_rep3_overlap.bed"

    # ========================================================
    # Combine pairwise-supported regions
    #
    # A region supported by:
    #
    #   Rep1 + Rep2
    #        OR
    #   Rep1 + Rep3
    #        OR
    #   Rep2 + Rep3
    #
    # is considered supported by >=2/3 replicates.
    #
    # Merge overlapping regions to generate a non-redundant
    # confident peak set.
    # ========================================================

    cat \
        "$conf_dir/${condition}_rep1_rep2_overlap.bed" \
        "$conf_dir/${condition}_rep1_rep3_overlap.bed" \
        "$conf_dir/${condition}_rep2_rep3_overlap.bed" \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > "$conf_dir/LNCD_CHIRP-${condition}_confident_peaks.bed"

    # ========================================================
    # Report peak counts
    # ========================================================

    echo ""
    echo "$condition individual replicate peaks:"
    echo "  Rep1: $(wc -l < "$rep1")"
    echo "  Rep2: $(wc -l < "$rep2")"
    echo "  Rep3: $(wc -l < "$rep3")"

    echo ""
    echo "$condition pairwise overlaps (>=$min_overlap bp):"

    echo "  Rep1 vs Rep2: \
$(wc -l < "$conf_dir/${condition}_rep1_rep2_overlap.bed")"

    echo "  Rep1 vs Rep3: \
$(wc -l < "$conf_dir/${condition}_rep1_rep3_overlap.bed")"

    echo "  Rep2 vs Rep3: \
$(wc -l < "$conf_dir/${condition}_rep2_rep3_overlap.bed")"

    echo ""
    echo "FINAL confident $condition peaks (>=2/3 replicates):"
    echo "  $(wc -l < "$conf_dir/LNCD_CHIRP-${condition}_confident_peaks.bed")"

    echo ""
}

# ============================================================
# Run D condition
# ============================================================

find_confident_peaks D


# ============================================================
# Run V condition
# ============================================================

find_confident_peaks V


# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 5b complete"
echo "========================================"

echo ""
echo "Final D confident peaks:"
echo "$conf_dir/LNCD_CHIRP-D_confident_peaks.bed"

echo ""
echo "Final V confident peaks:"
echo "$conf_dir/LNCD_CHIRP-V_confident_peaks.bed"

echo ""
echo "All Step 5B files are in:"
echo "$conf_dir" 
  


