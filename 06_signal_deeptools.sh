#!/usr/bin/env bash
# =============================================================================
# STEP 6 - Make signal tracks and plots with deepTools (the final figures).
#
# This step:
#   1) generates CPM-normalized bigWig tracks
#   2) summarizes signal over the Step-4 peaks
#   3) generates heatmaps and profiles
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 06_signal_deeptools.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_deeptools
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

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"
bam_dir="$out/bam"
peak_dir="$out/peaks"
bw_dir="$out/bigwig"       # output: signal tracks
fig_dir="$out/figures"     # output: plots
threads=16

mkdir -p "$bw_dir" "$fig_dir"

# ============================================================
# A) Generate CPM-normalized bigWig tracks
# ============================================================

echo "=== Generating CPM-normalized bigWig tracks ==="

# CHIRP-D
bamCoverage \
    -b "$bam_dir/CHIRP-D_merged.bam" \
    -o "$bw_dir/CHIRP-D_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

# Input-D
bamCoverage \
    -b "$bam_dir/Input-D_merged.bam" \
    -o "$bw_dir/Input-D_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

# CHIRP-V
bamCoverage \
    -b "$bam_dir/CHIRP-V_merged.bam" \
    -o "$bw_dir/CHIRP-V_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

# Input-v
bamCoverage \
    -b "$bam_dir/Input-v_merged.bam" \
    -o "$bw_dir/Input-v_CPM.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    -p "$threads"

echo "BigWig generation complete."
# ============================================================
# B) Define Step-4 peak files
# ============================================================

D_peaks="$peak_dir/LNCD_CHIRP-D_FC5_noBL.narrowPeak"
V_peaks="$peak_dir/LNCD_CHIRP-V_FC5_noBL.narrowPeak"

# Check files exist

for f in "$D_peaks" "$V_peaks"; do

    if [[ ! -f "$f" ]]; then
        echo "ERROR: Cannot find peak file:"
        echo "$f"
        exit 1
    fi

done


echo ""
echo "Step-4 peak counts:"
echo "D peaks: $(wc -l < "$D_peaks")"
echo "V peaks: $(wc -l < "$V_peaks")"

# ============================================================
# C) Summarize D and Input-D signal over D peaks
# ============================================================
echo ""
echo "=== Summarizing signal over D peaks ==="

multiBigwigSummary BED-file \
    --BED "$D_peaks" \
    -b \
        "$bw_dir/CHIRP-D_CPM.bw" \
        "$bw_dir/Input-D_CPM.bw" \
    --labels CHIRP-D Input-D \
    -o "$peak_dir/D_peak_signal.npz" \
    --outRawCounts "$peak_dir/D_peak_signal.tab" \
    -p "$threads"

# ============================================================
# D) Summarize V and Input-v signal over V peaks
# ============================================================

echo ""
echo "=== Summarizing signal over V peaks ==="

multiBigwigSummary BED-file \
    --BED "$V_peaks" \
    -b \
        "$bw_dir/CHIRP-V_CPM.bw" \
        "$bw_dir/Input-v_CPM.bw" \
    --labels CHIRP-V Input-v \
    -o "$peak_dir/V_peak_signal.npz" \
    --outRawCounts "$peak_dir/V_peak_signal.tab" \
    -p "$threads"

# ============================================================
# E) Heatmap around D peaks
#
# Show:
#   CHIRP-D
#   Input-D
#
# +/- 3 kb around peak center
# ============================================================

echo ""
echo "=== Generating D peak heatmap ==="

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
    -o "$fig_dir/D_peak_matrix.gz" \
    -p "$threads"

plotHeatmap \
    -m "$fig_dir/D_peak_matrix.gz" \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --refPointLabel "peak center" \
    --heatmapHeight 12 \
    -out "$fig_dir/D_peak_heatmap.png"

# ============================================================
# F) Average profile around D peaks
# ============================================================

plotProfile \
    -m "$fig_dir/D_peak_matrix.gz" \
    --refPointLabel "peak center" \
    --plotTitle "Signal around D peaks" \
    -out "$fig_dir/D_peak_profile.png"

# ============================================================
# G) Heatmap around V peaks
#
# Show:
#   CHIRP-V
#   Input-v
# ============================================================

echo ""
echo "=== Generating V peak heatmap ==="

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
    -o "$fig_dir/V_peak_matrix.gz" \
    -p "$threads"

plotHeatmap \
    -m "$fig_dir/V_peak_matrix.gz" \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --refPointLabel "peak center" \
    --heatmapHeight 12 \
    -out "$fig_dir/V_peak_heatmap.png"

# ============================================================
# H) Average profile around V peaks
# ============================================================

plotProfile \
    -m "$fig_dir/V_peak_matrix.gz" \
    --refPointLabel "peak center" \
    --plotTitle "Signal around V peaks" \
    -out "$fig_dir/V_peak_profile.png"

# ============================================================
# Finished
# ============================================================

echo ""
echo "========================================"
echo "STEP 5 complete"
echo "========================================"

echo ""
echo "BigWig tracks:"
echo "$bw_dir"

echo ""
echo "D peak set:"
echo "$D_peaks"

echo ""
echo "V peak set:"
echo "$V_peaks"

echo ""
echo "Figures:"
echo "$fig_dir"\
