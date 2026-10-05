#!/usr/bin/env bash
# ============================================================
# STEP 06 - deepTools signal analysis around confident peaks
#
# Experimental design:
#   V = control
#   D = treatment
#
# Input:
#   Merged BAM files from Step 3
#   Confident peaks from Step 5b
#
# This step:
#   1) Generate CPM-normalized bigWig tracks
#   2) Generate heatmap/profile for D confident peaks
#   3) Generate heatmap/profile for V confident peaks
#   4) Create a D/V union peak set
#   5) Compare D, V and their Inputs over the SAME union peaks
#
# NOTE:
#   D/V Pearson correlation is NOT performed because D and V
#   represent treatment and control conditions.
#
# Later:
#   Differential binding analysis can compare D vs V using
#   replicate-level read counts over a common peak set.
# ============================================================
# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_confident_deeptools
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths -------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

bam_dir="$out/bam"

# Step 5b confident peaks (INPUT only)
conf_dir="$out/peaks/confident"

# Step 06 output directories
bw_dir="$out/bigwig/06_confident_signal_deeptools"
fig_dir="$out/figures/06_confident_signal_deeptools"
matrix_dir="$out/matrix/06_confident_signal_deeptools"
peak6_dir="$out/peaks/06_confident_signal_deeptools"
threads=16

mkdir -p "$bw_dir" "$fig_dir" "$matrix_dir" "$peak6_dir"

# ---- Confident peak files ----------------------------------

D_peaks="$conf_dir/LNCD_CHIRP-D_confident_peaks.bed"
V_peaks="$conf_dir/LNCD_CHIRP-V_confident_peaks.bed"

# ============================================================
# Check input files
# ============================================================

for f in \
    "$bam_dir/CHIRP-D_merged.bam" \
    "$bam_dir/Input-D_merged.bam" \
    "$bam_dir/CHIRP-V_merged.bam" \
    "$bam_dir/Input-v_merged.bam" \
    "$D_peaks" \
    "$V_peaks"
do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Cannot find:"
        echo "$f"
        exit 1
    fi
done

# ============================================================
# A) Generate CPM-normalized bigWig tracks
# ============================================================

echo ""
echo "========================================"
echo "A) Generating CPM-normalized bigWig tracks"
echo "========================================"

# Treatment D
bamCoverage \
    -b "$bam_dir/CHIRP-D_merged.bam" \
    -o "$bw_dir/CHIRP-D_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

# Input for treatment D
bamCoverage \
    -b "$bam_dir/Input-D_merged.bam" \
    -o "$bw_dir/Input-D_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

# Control V
bamCoverage \
    -b "$bam_dir/CHIRP-V_merged.bam" \
    -o "$bw_dir/CHIRP-V_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

# Input for control V
bamCoverage \
    -b "$bam_dir/Input-v_merged.bam" \
    -o "$bw_dir/Input-v_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

# ============================================================
# B) Treatment D signal around D confident peaks
# ============================================================

echo ""
echo "========================================"
echo "B) D signal around D confident peaks"
echo "========================================"

computeMatrix reference-point \
    --referencePoint center \
    -S \
        "$bw_dir/CHIRP-D_CPM.bw" \
        "$bw_dir/Input-D_CPM.bw" \
    -R "$D_peaks" \
    -b 3000 \
    -a 3000 \
    --binSize 50 \
    --samplesLabel CHIRP-D Input-D \
    --skipZeros \
    -o "$matrix_dir/D_confident_matrix.gz" \
    -p "$threads"

plotHeatmap \
    -m "$matrix_dir/D_confident_matrix.gz" \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --refPointLabel "peak center" \
    --plotTitle "LNC-D CHIRP-D confident peaks" \
    -out "$fig_dir/D_confident_heatmap.png"

plotProfile \
    -m "$matrix_dir/D_confident_matrix.gz" \
    --refPointLabel "peak center" \
    --plotTitle "CHIRP-D signal around confident peaks" \
    -out "$fig_dir/D_confident_profile.png"

# ============================================================
# C) Control V signal around V confident peaks
# ============================================================

echo ""
echo "========================================"
echo "C) V signal around V confident peaks"
echo "========================================"

computeMatrix reference-point \
    --referencePoint center \
    -S \
        "$bw_dir/CHIRP-V_CPM.bw" \
        "$bw_dir/Input-v_CPM.bw" \
    -R "$V_peaks" \
    -b 3000 \
    -a 3000 \
    --binSize 50 \
    --samplesLabel CHIRP-V Input-v \
    --skipZeros \
    -o "$matrix_dir/V_confident_matrix.gz" \
    -p "$threads"

plotHeatmap \
    -m "$matrix_dir/V_confident_matrix.gz" \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --refPointLabel "peak center" \
    --plotTitle "LNC-D CHIRP-V confident peaks" \
    -out "$fig_dir/V_confident_heatmap.png"

plotProfile \
    -m "$matrix_dir/V_confident_matrix.gz" \
    --refPointLabel "peak center" \
    --plotTitle "CHIRP-V signal around confident peaks" \
    -out "$fig_dir/V_confident_profile.png"

# ============================================================
# D) Create D/V union peak set
#
# Combine treatment and control confident peaks.
#
# This creates one common set of genomic regions so D and V
# can be visualized over exactly the same coordinates.
# ============================================================

echo ""
echo "========================================"
echo "D) Creating D/V union peak set"
echo "========================================"

cat "$D_peaks" "$V_peaks" \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > "$peak6_dir/LNCD_CHIRP-DV_union_peaks.bed"

DV_peaks="$peak6_dir/LNCD_CHIRP-DV_union_peaks.bed"

echo ""
echo "Peak numbers:"
echo "  D confident peaks: $(wc -l < "$D_peaks")"
echo "  V confident peaks: $(wc -l < "$V_peaks")"
echo "  D/V union peaks:   $(wc -l < "$DV_peaks")"

# ============================================================
# E) Compare D/V/Input signal over the SAME union peak set
#
# Four tracks:
#
#   CHIRP-D   = treatment
#   CHIRP-V   = control
#   Input-D   = treatment input
#   Input-v   = control input
#
# +/- 3 kb around the center of each union peak.
# ============================================================

echo ""
echo "========================================"
echo "E) D/V signal over union peaks"
echo "========================================"

computeMatrix reference-point \
    --referencePoint center \
    -S \
        "$bw_dir/CHIRP-D_CPM.bw" \
        "$bw_dir/CHIRP-V_CPM.bw" \
        "$bw_dir/Input-D_CPM.bw" \
        "$bw_dir/Input-v_CPM.bw" \
    -R "$DV_peaks" \
    -b 3000 \
    -a 3000 \
    --binSize 50 \
    --samplesLabel \
        CHIRP-D CHIRP-V Input-D Input-v \
    --skipZeros \
    -o "$matrix_dir/DV_union_matrix.gz" \
    -p "$threads"

# ============================================================
# F) Combined heatmap
# ============================================================

echo ""
echo "========================================"
echo "F) Creating D/V union heatmap"
echo "========================================"

plotHeatmap \
    -m "$matrix_dir/DV_union_matrix.gz" \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --refPointLabel "peak center" \
    --plotTitle "LNC-D ChIRP-seq: treatment vs control" \
    -out "$fig_dir/DV_union_heatmap.png"

# ============================================================
# G) Combined average profile
# ============================================================

echo ""
echo "========================================"
echo "G) Creating D/V union profile"
echo "========================================"

plotProfile \
    -m "$matrix_dir/DV_union_matrix.gz" \
    --refPointLabel "peak center" \
    --plotTitle "LNC-D ChIRP-seq signal around confident peaks" \
    -out "$fig_dir/DV_union_profile.png"

# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 06 COMPLETE"
echo "========================================"

echo ""
echo "Step 06 outputs:"
echo ""

echo "BigWig tracks:"
echo "  $bw_dir"

echo ""
echo "Matrices:"
echo "  $matrix_dir"

echo ""
echo "Figures:"
echo "  $fig_dir"

echo ""
echo "Step 06 peak sets:"
echo "  $peak6_dir"

echo ""
echo "D/V union peaks:"
echo "  $DV_peaks"