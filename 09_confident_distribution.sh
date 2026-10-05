#!/usr/bin/env bash
# ============================================================
# STEP 09 - Genomic distribution of confident D peaks
#
# Input:
#   HOMER annotation from Step 08
#
# Output:
#   Genomic category for each peak
#   Counts and percentages for each category
# ============================================================
#SBATCH --job-name=09_confident_distribution
#SBATCH --account=htl145
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/tscc/lustre/ddn/scratch/jiw169/chirp_seq_analysis/LNC-D/%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiw169@ucsd.edu
# Activate software
# ============================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate CHIRP-SEQ

# ---- Paths -----------------------------------------------------------------

out="/tscc/lustre/ddn/scratch/$USER/chirp_seq_analysis/LNC-D"

input="$out/confident_peak_annotation/08_confident_D/confident_D_HOMER_annotation.tsv"

step9_dir="$out/genomic_distribution/09_confident_D"

mkdir -p "$step9_dir"


# ============================================================
# Check input
# ============================================================

if [[ ! -s "$input" ]]; then

    echo "ERROR: Cannot find Step 08 annotation:"
    echo "$input"

    exit 1

fi

# ============================================================
# Find HOMER Annotation column
# ============================================================

annotation_col=$(awk -F'\t' '
NR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="Annotation") {
            print i
            exit
        }
    }
}' "$input")


if [[ -z "$annotation_col" ]]; then

    echo "ERROR: Cannot find Annotation column."

    exit 1

fi


echo ""
echo "Annotation column:"
echo "$annotation_col"


# ============================================================
# Classify genomic regions
#
# Output:
# peak_ID    HOMER_annotation    genomic_category
# ============================================================

awk -F'\t' \
    -v OFS='\t' \
    -v col="$annotation_col" '

NR==1 {

    print "Peak_ID","HOMER_Annotation","Genomic_Category"

    next
}

{

    annotation=$col

    lower=tolower(annotation)


    if(lower ~ /promoter/) {

        category="Promoter"

    } else if(lower ~ /exon/) {

        category="Exon"

    } else if(lower ~ /intron/) {

        category="Intron"

    } else if(lower ~ /intergenic/) {

        category="Intergenic"

    } else {

        category="Other"

    }


    print $1,annotation,category

}' "$input" \
> "$step9_dir/confident_D_genomic_categories.tsv"


# ============================================================
# Calculate counts
# ============================================================

awk -F'\t' \
    -v OFS='\t' '

NR>1 {

    count[$3]++

    total++

}

END {

    print "Category","Count","Percentage"

    categories[1]="Promoter"
    categories[2]="Exon"
    categories[3]="Intron"
    categories[4]="Intergenic"
    categories[5]="Other"


    for(i=1;i<=5;i++) {

        c=categories[i]

        n=count[c]+0

        if(total > 0) {
            pct=n/total*100
        } else {
            pct=0
        }

        printf "%s\t%d\t%.2f\n",c,n,pct
    }

}' \
"$step9_dir/confident_D_genomic_categories.tsv" \
> "$step9_dir/confident_D_genomic_distribution.tsv"


# ============================================================
# Print results
# ============================================================

echo ""
echo "========================================"
echo "STEP 09 - Genomic distribution"
echo "========================================"

echo ""

column -t \
"$step9_dir/confident_D_genomic_distribution.tsv"


echo ""
echo "Total confident D peaks:"

awk -F'\t' '
NR>1 {
    total += $2
}
END {
    print total
}' \
"$step9_dir/confident_D_genomic_distribution.tsv"


echo ""
echo "========================================"
echo "STEP 09 COMPLETE"
echo "========================================"

echo ""
echo "Peak-level categories:"
echo "$step9_dir/confident_D_genomic_categories.tsv"

echo ""
echo "Distribution summary:"
echo "$step9_dir/confident_D_genomic_distribution.tsv"