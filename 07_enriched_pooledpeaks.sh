#!/usr/bin/env bash
# ============================================================
# STEP 07 - Identify D-specific and D-enriched pooled peaks
#
# D = treatment
# V = control
#
# Inputs:
#
#   Step 4 pooled peaks:
#     LNCD_CHIRP-D_FC5_noBL.narrowPeak
#     LNCD_CHIRP-V_FC5_noBL.narrowPeak
#
#   Step 06 CPM-normalized bigWig tracks:
#     CHIRP-D_CPM.bw
#     CHIRP-V_CPM.bw
#
# Definitions:
#
# 1) D-specific peak:
#
#    A D peak with NO genomic overlap with any V peak.
#
#    This is based on peak calling.
#
# 2) D-enriched peak:
#
#    A D peak with:
#
#    mean CHIRP-D CPM / mean CHIRP-V CPM >= 2
#
#    This is based on signal intensity.
#
#    D-enriched peaks can include:
#      - D-specific peaks
#      - D/V-shared peaks
#
#
# IMPORTANT:
#
# D-specific and D-enriched are independent categories.
# A peak can therefore be both D-specific and D-enriched.
#
# These are pooled/descriptive analyses and are NOT
# statistical differential-binding tests.
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=07_enriched_pooled_peaks
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

# ============================================================
# Paths
# ============================================================

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

peak_dir="$out/peaks"
bw_dir="$out/bigwig"

step7_dir="$peak_dir/07_D_specific_and_enriched_pooled"

mkdir -p "$step7_dir"

# Parameter
# ============================================================

enrichment_cutoff=2

# Input peak files from Step 4
# ============================================================

D_np="$peak_dir/LNCD_CHIRP-D_FC5_noBL.narrowPeak"

V_np="$peak_dir/LNCD_CHIRP-V_FC5_noBL.narrowPeak"

# Input bigWig files from Step 06
# ============================================================

D_bw="$bw_dir/CHIRP-D_CPM.bw"

V_bw="$bw_dir/CHIRP-V_CPM.bw"

# Check inputs
# ============================================================

for f in "$D_np" "$V_np" "$D_bw" "$V_bw"; do

    if [[ ! -s "$f" ]]; then

        echo "ERROR: Cannot find or empty file:"
        echo "$f"

        exit 1
    fi

done


# Check software

for cmd in bedtools multiBigwigSummary; do

    if ! command -v "$cmd" >/dev/null 2>&1; then

        echo "ERROR: $cmd not found."
        exit 1

    fi

done


echo ""
echo "========================================"
echo "STEP 07"
echo "D-specific + D-enriched pooled peaks"
echo "========================================"

echo ""
echo "D = treatment"
echo "V = control"

echo ""
echo "D-specific:"
echo "  no overlap with any V peak"

echo ""
echo "D-enriched:"
echo "  D/V CPM >= $enrichment_cutoff"

# A) Convert pooled narrowPeak files to BED4
# ============================================================

echo ""
echo "========================================"
echo "A) Preparing pooled peak BED files"
echo "========================================"

awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$4
}' "$D_np" \
    | sort -k1,1 -k2,2n \
    > "$step7_dir/D_pooled_peaks.bed"


awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$4
}' "$V_np" \
    | sort -k1,1 -k2,2n \
    > "$step7_dir/V_pooled_peaks.bed"


D_bed="$step7_dir/D_pooled_peaks.bed"

V_bed="$step7_dir/V_pooled_peaks.bed"


echo ""
echo "D pooled peaks:"
wc -l "$D_bed"

echo ""
echo "V pooled peaks:"
wc -l "$V_bed"


# ============================================================
# B) Find D peaks overlapping V
#
# Any genomic overlap is considered shared.
#
# -wa : return D peak
# -u  : report each D peak only once
# ============================================================

echo ""
echo "========================================"
echo "B) Finding D peaks shared with V"
echo "========================================"


bedtools intersect \
    -a "$D_bed" \
    -b "$V_bed" \
    -wa -u \
    > "$step7_dir/D_shared_with_V.bed"


echo ""
echo "D peaks overlapping V:"
wc -l "$step7_dir/D_shared_with_V.bed"


# C) Identify D-specific peaks
#
# D-specific =
# D peak with NO overlap with any V peak.
#
# -v returns A peaks that have no overlap with B.
# ============================================================

echo ""
echo "========================================"
echo "C) Identifying D-specific peaks"
echo "========================================"


bedtools intersect \
    -a "$D_bed" \
    -b "$V_bed" \
    -v \
    > "$step7_dir/D_specific_peaks.bed"


# Add new IDs

awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,"D_specific_"NR
}' "$step7_dir/D_specific_peaks.bed" \
    > "$step7_dir/D_specific_peaks_ID.bed"


echo ""
echo "D-specific peaks:"
wc -l "$step7_dir/D_specific_peaks_ID.bed"


# D) Measure D and V CPM over ALL D peaks
#
# IMPORTANT:
#
# Both D and V signal are measured over exactly the SAME
# genomic coordinates: the D pooled peak regions.
#
# This allows direct D/V signal comparison.
# ============================================================

echo ""
echo "========================================"
echo "D) Measuring D and V signal over all D peaks"
echo "========================================"


multiBigwigSummary BED-file \
    --BED "$D_bed" \
    -b \
        "$D_bw" \
        "$V_bw" \
    --labels CHIRP-D CHIRP-V \
    -o "$step7_dir/D_peaks_DV_signal.npz" \
    --outRawCounts "$step7_dir/D_peaks_DV_signal_raw.tsv"


# ============================================================
# E) Create coordinate -> PeakID lookup
#
# multiBigwigSummary keeps coordinates but not our BED4 ID.
# ============================================================

awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$4
}' "$D_bed" \
    > "$step7_dir/D_peak_coordinate_ID.tsv"


# F) Create D/V signal table
#
# Output:
#
# Chr
# Start
# End
# PeakID
# D_CPM
# V_CPM
# D_over_V
# ============================================================

echo ""
echo "========================================"
echo "E) Calculating D/V signal ratios"
echo "========================================"


grep -v '^#' \
    "$step7_dir/D_peaks_DV_signal_raw.tsv" \
    > "$step7_dir/D_peaks_DV_signal_noheader.tmp"


awk '
BEGIN {
    OFS="\t"
}

NR==FNR {

    key=$1 FS $2 FS $3
    id[key]=$4

    next
}

{

    key=$1 FS $2 FS $3

    d=$4+0
    v=$5+0

    if (v > 0) {

        ratio=d/v

    } else if (d > 0) {

        ratio="Inf"

    } else {

        ratio=0

    }

    print $1,$2,$3,id[key],d,v,ratio
}
' \
    "$step7_dir/D_peak_coordinate_ID.tsv" \
    "$step7_dir/D_peaks_DV_signal_noheader.tmp" \
    > "$step7_dir/D_peaks_DV_signal_body.tsv"


{
    printf "Chr\tStart\tEnd\tPeakID\tD_CPM\tV_CPM\tD_over_V\n"

    cat "$step7_dir/D_peaks_DV_signal_body.tsv"

} > "$step7_dir/D_peaks_DV_signal.tsv"


rm -f \
    "$step7_dir/D_peaks_DV_signal_noheader.tmp" \
    "$step7_dir/D_peaks_DV_signal_body.tsv"


# ============================================================
# G) Identify D-enriched peaks
#
# D-enriched =
#
# mean D CPM / mean V CPM >= 2
#
# IMPORTANT:
#
# This is tested across ALL D peaks.
#
# A V peak does NOT need to have been called by MACS2.
#
# Therefore D-enriched peaks can include:
#
#   D-specific peaks
#       +
#   D/V-shared peaks
# ============================================================

echo ""
echo "========================================"
echo "F) Identifying D-enriched peaks"
echo "========================================"


awk -v cutoff="$enrichment_cutoff" '
BEGIN {
    FS=OFS="\t"
}

NR==1 {
    next
}

{

    d=$5+0
    v=$6+0

    if (
        (v > 0 && d/v >= cutoff) ||
        (v == 0 && d > 0)
    ) {

        print $1,$2,$3,$4,$5,$6,$7

    }
}
' "$step7_dir/D_peaks_DV_signal.tsv" \
    > "$step7_dir/D_enriched_peaks_signal.tsv"


# H) Generate BED file for D-enriched peaks
# ============================================================

awk 'BEGIN{OFS="\t"} {

    print $1,$2,$3,"D_enriched_"NR

}' "$step7_dir/D_enriched_peaks_signal.tsv" \
    > "$step7_dir/D_enriched_peaks_ID.bed"


echo ""
echo "D-enriched peaks:"
wc -l "$step7_dir/D_enriched_peaks_ID.bed"


# I) Determine which D-enriched peaks are also D-specific
#
# This is useful because the categories can overlap.
# ============================================================

echo ""
echo "========================================"
echo "G) Comparing D-specific and D-enriched sets"
echo "========================================"


bedtools intersect \
    -a "$step7_dir/D_enriched_peaks_ID.bed" \
    -b "$step7_dir/D_specific_peaks_ID.bed" \
    -wa -u \
    > "$step7_dir/D_specific_AND_enriched.bed"


# ============================================================
# J) D-enriched peaks that overlap V
#
# These are shared peak calls but quantitatively stronger in D.
# ============================================================

bedtools intersect \
    -a "$step7_dir/D_enriched_peaks_ID.bed" \
    -b "$V_bed" \
    -wa -u \
    > "$step7_dir/D_enriched_shared_with_V.bed"


# ============================================================
# K) Create non-redundant treatment-associated peak set
#
# Union of:
#
#   D-specific
#       OR
#   D-enriched
#
# Because some peaks are in BOTH categories, bedtools merge
# removes duplicate/overlapping coordinates.
# ============================================================

cat \
    "$step7_dir/D_specific_peaks_ID.bed" \
    "$step7_dir/D_enriched_peaks_ID.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > "$step7_dir/D_treatment_associated_peaks.bed"


# Add IDs to the union set

awk 'BEGIN{OFS="\t"} {

    print $1,$2,$3,"D_treatment_"NR

}' "$step7_dir/D_treatment_associated_peaks.bed" \
    > "$step7_dir/D_treatment_associated_peaks_ID.bed"

# L) Count categories
# ============================================================

n_D=$(wc -l < "$D_bed")

n_V=$(wc -l < "$V_bed")

n_shared=$(wc -l < "$step7_dir/D_shared_with_V.bed")

n_specific=$(wc -l < "$step7_dir/D_specific_peaks_ID.bed")

n_enriched=$(wc -l < "$step7_dir/D_enriched_peaks_ID.bed")

n_both=$(wc -l < "$step7_dir/D_specific_AND_enriched.bed")

n_enriched_shared=$(wc -l < "$step7_dir/D_enriched_shared_with_V.bed")

n_treatment=$(wc -l < "$step7_dir/D_treatment_associated_peaks_ID.bed")

# M) Print summary
# ============================================================

echo ""
echo "========================================"
echo "STEP 07 SUMMARY"
echo "========================================"

echo ""

echo "D total pooled peaks:             $n_D"

echo "V total pooled peaks:             $n_V"

echo "D peaks overlapping V:            $n_shared"

echo "D-specific peaks:                 $n_specific"

echo "D-enriched peaks (D/V >= 2):      $n_enriched"

echo "Specific AND enriched:            $n_both"

echo "Enriched + shared with V:         $n_enriched_shared"

echo "Treatment-associated union:       $n_treatment"

# N) Save summary
# ============================================================

{
    printf "Category\tCount\n"

    printf "D_total\t%s\n" "$n_D"

    printf "V_total\t%s\n" "$n_V"

    printf "D_shared_with_V\t%s\n" "$n_shared"

    printf "D_specific\t%s\n" "$n_specific"

    printf "D_enriched_DoverV_ge_2\t%s\n" "$n_enriched"

    printf "D_specific_AND_enriched\t%s\n" "$n_both"

    printf "D_enriched_shared_with_V\t%s\n" "$n_enriched_shared"

    printf "D_treatment_associated_union\t%s\n" "$n_treatment"

} > "$step7_dir/D_specific_enriched_summary.tsv"

# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 07 COMPLETE"
echo "========================================"

echo ""
echo "Output directory:"
echo "$step7_dir"

echo ""
echo "1. D-specific peaks:"
echo "$step7_dir/D_specific_peaks_ID.bed"

echo ""
echo "2. D-enriched peaks:"
echo "$step7_dir/D_enriched_peaks_ID.bed"

echo ""
echo "3. D-enriched signal table:"
echo "$step7_dir/D_enriched_peaks_signal.tsv"

echo ""
echo "4. Signal for ALL D peaks:"
echo "$step7_dir/D_peaks_DV_signal.tsv"

echo ""
echo "5. D-specific AND D-enriched:"
echo "$step7_dir/D_specific_AND_enriched.bed"

echo ""
echo "6. D-enriched peaks shared with V:"
echo "$step7_dir/D_enriched_shared_with_V.bed"

echo ""
echo "7. Non-redundant treatment-associated peaks:"
echo "$step7_dir/D_treatment_associated_peaks_ID.bed"

echo ""
echo "8. Summary:"
echo "$step7_dir/D_specific_enriched_summary.tsv"