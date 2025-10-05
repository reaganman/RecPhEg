#!/bin/bash
RESULTS_DIR=$1
GENONMES=$2
CLUSTER_THRESHOLD=$3


# Initial lovis4u run to generate base feature annotation table
lovis4u -gff $RESULTS_DIR/pharokka_results/single_gffs/ \
    -o $RESULTS_DIR/lovis4u/ \
    -scc -hl -sxa -c A4L

# Create separate annotation tables for each possible combo group
python Scripts/update_annotations_table.py \
    --overlaps $RESULTS_DIR/shared_kmers.csv \
    --annotations $RESULTS_DIR/lovis4u/feature_annotation_table.tsv \
    --out $RESULTS_DIR/lovis4u/feature_annotation_table_updated.tsv \
    --genomes $GENONMES \
    --cluster_threshold $CLUSTER_THRESHOLD

# Path to combo groups tracking file
COMBO_GROUPS=$RESULTS_DIR/lovis4u/feature_annotation_table_updated_combo_groups.csv

# Loop over each generated combo annotation table and run lovis4u separately
while IFS=, read -r group accessions; do
    # skip header line
    if [ "$group" == "group" ]; then
        continue
    fi

    # Clean accessions field → space-separated list
    CLEAN_ACCESSIONS=$(echo "$accessions" | sed 's/\[//g; s/\]//g; s/"//g; s/'\''//g; s/,/ /g')

    # Make tmp GFF directory for this group
    GROUP_GFF_DIR="$RESULTS_DIR/pharokka_results/${group}_gffs"
    mkdir -p "$GROUP_GFF_DIR"

    # Copy relevant GFFs
    for acc in $CLEAN_ACCESSIONS; do
        GFF_FILE="$RESULTS_DIR/pharokka_results/single_gffs/${acc}.gff"
        if [ -f "$GFF_FILE" ]; then
            cp "$GFF_FILE" "$GROUP_GFF_DIR/"
        else
            echo "⚠️ Warning: GFF file not found for $acc ($GFF_FILE)"
        fi
    done


    # Run lovis4u with group-specific annotation + GFF subset

    mkdir -p $RESULTS_DIR/lovis4u_overlaps
    ANNO_FILE="$RESULTS_DIR/lovis4u/feature_annotation_table_updated_${group}.tsv"
    OUT_DIR="$RESULTS_DIR/lovis4u_overlaps/${group}"

    echo "▶ Running lovis4u for $group with $ANNO_FILE"
    lovis4u -gff "$GROUP_GFF_DIR" \
        -o "$OUT_DIR" \
        -scc -sxa -hl -c A4L \
        --feature-annotation-file "$ANNO_FILE"

done < "$COMBO_GROUPS"