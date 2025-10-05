#!/bin/bash

usage() {
    echo "Usage: $0 -i <INPUT_FASTA> -d <PHAROKKA_DBS> -o <OVERLAP_SIZE> [-O <OUTPUT_DIR>] [-t <CLUSTER_THRESHOLD>]"
    echo "Example: bash find_assembly_overlaps.sh -i /home/remc1874/phage_genomics/fastas/bas14-18.fasta -d /scratch/alpine/remc1874/pharokka -o 30 -t 5"
    exit 1
}

# Parse arguments
while getopts ":i:d:o:O:t:" opt; do
  case $opt in
    i) INPUT_FASTA="$OPTARG" ;;
    d) PHAROKKA_DBS="$OPTARG" ;;
    o) OVERLAP_SIZE="$OPTARG" ;;
    O) OUTPUT_DIR="$OPTARG" ;;
    t) CLUSTER_THRESHOLD="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check required arguments
if [ -z "$INPUT_FASTA" ] || [ -z "$PHAROKKA_DBS" ] || [ -z "$OVERLAP_SIZE" ]; then
    echo "Error: -i, -d, and -o arguments are required."
    usage
fi

# If OUTPUT_DIR not provided, default to ./results_<basename of input fasta>
if [ -z "$OUTPUT_DIR" ]; then
    BASENAME=$(basename "$INPUT_FASTA" .fasta)
    BASENAME=$(basename "$BASENAME" .fa)  # handle .fa as well
    OUTPUT_DIR="./results_${BASENAME}"
fi

echo "Running with:"
echo "  INPUT_FASTA       = $INPUT_FASTA"
echo "  PHAROKKA_DBS      = $PHAROKKA_DBS"
echo "  OVERLAP_SIZE      = $OVERLAP_SIZE"
echo "  OUTPUT_DIR        = $OUTPUT_DIR"
echo "  CLUSTER_THRESHOLD = ${CLUSTER_THRESHOLD:-not set}"

# Create main output directory
mkdir -p "$OUTPUT_DIR"

# Activate pharokka environment
module load anaconda
conda activate pharokka_env

# Run pharokka
export PYTHONWARNINGS="ignore"
pharokka.py -i "$INPUT_FASTA" -o "$OUTPUT_DIR/pharokka_results" -d "$PHAROKKA_DBS" -t 8 -f --meta --split --skip_mash

# Run mmseqs2
conda activate mmseqs2

FAA_DIR="$OUTPUT_DIR/pharokka_results/single_faas"

# Make output directories
mkdir -p "$OUTPUT_DIR/mmseqs/DB" "$OUTPUT_DIR/mmseqs/Results" "$OUTPUT_DIR/mmseqs/tmp"

# Combine all proteomes into one FASTA
cat "$FAA_DIR"/*.faa > "$OUTPUT_DIR/proteomes.faa"

# Create MMseqs2 database
mmseqs createdb "$OUTPUT_DIR/proteomes.faa" "$OUTPUT_DIR/mmseqs/DB/proteomesDB"

# Cluster proteins
mmseqs cluster "$OUTPUT_DIR/mmseqs/DB/proteomesDB" "$OUTPUT_DIR/mmseqs/Results/clusters" "$OUTPUT_DIR/mmseqs/tmp" \
    --cluster-mode 0 --cov-mode 0 -c 0.8 -s 7 --min-seq-id 0.25

# Convert results to TSV
mmseqs createtsv "$OUTPUT_DIR/mmseqs/DB/proteomesDB" "$OUTPUT_DIR/mmseqs/DB/proteomesDB" \
    "$OUTPUT_DIR/mmseqs/Results/clusters" "$OUTPUT_DIR/mmseqs/clusters.tsv"

# Make cluster fastas
conda activate find_assembly_overlaps_env

CLUSTER_THRESHOLD_ARG=""
if [ -n "$CLUSTER_THRESHOLD" ]; then
    CLUSTER_THRESHOLD_ARG="--cluster_threshold $CLUSTER_THRESHOLD"
fi

python Scripts/create_cluster_fastas.py \
    --genomes "$INPUT_FASTA" \
    --cds "$OUTPUT_DIR/pharokka_results/pharokka_cds_final_merged_output.tsv" \
    --clusters "$OUTPUT_DIR/mmseqs/clusters.tsv" \
    --outdir "$OUTPUT_DIR/cluster_fastas" \
    $CLUSTER_THRESHOLD_ARG

# Make cluster context region fastas
python Scripts/create_cluster_context_regions.py \
    --genomes "$INPUT_FASTA" \
    --cds "$OUTPUT_DIR/pharokka_results/pharokka_cds_final_merged_output.tsv" \
    --clusters "$OUTPUT_DIR/mmseqs/clusters.tsv" \
    --outdir "$OUTPUT_DIR/cluster_context_regions" \
    --flank "$OVERLAP_SIZE" \
    $CLUSTER_THRESHOLD_ARG

# Look for overlaps in cluster context regions
rm -f "$OUTPUT_DIR/shared_kmers.csv"
python Scripts/find_shared_kmers_clusters.py "$OUTPUT_DIR/cluster_context_regions/" -k "$OVERLAP_SIZE" -o $OUTPUT_DIR/shared_kmers.csv -n $CLUSTER_THRESHOLD

# Visualize assembly overlaps with lovis4u
conda activate ilovis4u_env
bash Scripts/lovis_overlaps.sh $OUTPUT_DIR $INPUT_FASTA $CLUSTER_THRESHOLD
