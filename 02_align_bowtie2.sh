#!/usr/bin/env bash
# =============================================================================
# STEP 2 - Align the cleaned reads to the human genome (GRCh38) with bowtie2.
#
# "Alignment" means figuring out where in the genome each short read came from.
# The output is a BAM file (a compressed table of aligned reads), sorted by
# genome position so later tools can use it quickly.
#
# HOW TO RUN (from a TSCC login node):
#     sbatch 02_align_bowtie2.sh
# =============================================================================

# ---- Slurm settings ---------------------------------------------------------
#SBATCH --job-name=LNCD_align
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24         # hotel nodes allow up to 28 cores per node
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type END
#SBATCH --mail-user jiw169@ucsd.edu

# ---- Turn on our software ---------------------------------------------------
source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths ------------------------------------------------------------------
out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"
fastq_dir="$out/fastq_trimmed"    # input: cleaned reads from step 01b
bam_dir="$out/bam"                # output: aligned BAM files
ref_dir="$out/ref"                # reference genome/index directory
index="$out/ref/GRCh38_noalt_as"  # the genome index 
bt2_threads=16                    # cores for bowtie2 (the aligning)
sort_threads=6                    # cores for samtools sort (the sorting)
sort_mem="4G"                     # memory per sort thread (6 x 4G = 24G)

mkdir -p "$bam_dir"
mkdir -p "$ref_dir"

# ---- Download GRCh38 Bowtie2 index and hg38 blacklist -----------------------

# Download the GRCh38 Bowtie2 index only if it is missing.
if [ ! -e "$ref_dir/GRCh38.1.bt2" ] && [ ! -e "$ref_dir/GRCh38.1.bt2l" ]; then

    echo "=== downloading GRCh38 Bowtie2 index ==="
    cd "$ref_dir"
    wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
    unzip -o GRCh38_noalt_as.zip
    rm -f GRCh38_noalt_as.zip

    cd -
fi

# Download the hg38 ENCODE blacklist only if it is missing.
if [ ! -e "$ref_dir/hg38-blacklist.v2.bed" ]; then

    echo "=== downloading hg38 ENCODE blacklist ==="

    wget -O "$ref_dir/hg38-blacklist.v2.bed.gz" \
        https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz

    gunzip -f "$ref_dir/hg38-blacklist.v2.bed.gz"
fi


echo "Reference genome is ready in: $ref_dir"

# ---- Bowtie2 alignment -------------------------------------------------------

for sample in CHIRP-D1_S10_L002 \
    CHIRP-D2_S11_L002 \
    CHIRP-D3_S12_L002 \
    CHIRP-V1_S7_L002 \
    CHIRP-V2_S8_L002 \
    CHIRP-V3_S9_L002 \
    Input-D1_S4_L002 \
    Input-D2_S5_L002 \
    Input-D3_S6_L002 \
    Input-v1_S1_L002 \
    Input-v2_S2_L002 \
    Input-v3_S3_L002; do
    echo "=== aligning $sample ==="
    # bowtie2 aligns the reads and prints results; we pipe ("|") those results
    # straight into "samtools sort" so we never write the huge unsorted file.
    #   -x   the genome index
    #   -1/-2  the two paired read files
    #   -X 1000  allow read pairs up to 1000 bp apart
    #   --no-unal  drop reads that didn't align anywhere
    #   --mm  share the genome index in memory across samples (faster)
    bowtie2 \
        -x "$index" \
        -1 "$fastq_dir/${sample}_R1.fastq.gz" \
        -2 "$fastq_dir/${sample}_R2.fastq.gz" \
        -X 1000 \
        --no-unal \
        --mm \
        -p "$bt2_threads" \
        2> "$bam_dir/${sample}_bowtie2.log" \
      | samtools sort -@ "$sort_threads" -m "$sort_mem" \
            -T "$bam_dir/${sample}_sorttmp" \
            -o "$bam_dir/${sample}_sorted.bam" -

    # Build an index for the BAM (a .bai file) so tools can jump around it fast.
    samtools index -@ "$sort_threads" "$bam_dir/${sample}_sorted.bam"
    echo "Done: $sample"
done

echo "Aligned BAM files are in: $bam_dir"
