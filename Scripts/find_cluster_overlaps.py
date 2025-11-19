import argparse
import pandas as pd
from Bio import SeqIO, Seq
from collections import defaultdict
import re


def get_genome_key(gene):
    match = re.match(r"([A-Z]{2}_\d+\.\d+)", gene)
    return match.group(1) if match else gene.split("_")[0]

def shared_overlap(cluster_df, genomes, ref_point, direction, sense):
    """
    Find the overlap length across all members of a cluster starting at
    ref_point (CDS start or stop) and going in direction.

    
    Parameters:
        cluster_df (pd.DataFrame): cluster annotations
        genomes (list): SeqRecords of genomes
        ref_point (str): "start" or "stop"
        direction (str): "+" or "-"
        sense (str): "up" or "down"
    
    Returns:
        int: overlap length in nucleotides
    """
    genome_dict = {rec.id: rec.seq for rec in genomes}

    mem_genome_dict = {}
    for _, mem in cluster_df.iterrows():
        gene = mem["gene"]
        genome_key = get_genome_key(gene)

        mem_genome_dict[gene] = genome_dict.get(genome_key)

    base_mem = cluster_df.iloc[0]
    base_genome_key = get_genome_key(base_mem["gene"])
    base_genome = genome_dict[base_genome_key]

    ov_len = 0

    while True:
        all_overlap = True
        for _, mem in cluster_df.iterrows():
            seq = mem_genome_dict[mem["gene"]]
            coord = int(mem[ref_point]) - 1

            # Determine position index based on sense and strand
            if direction == "+":
                if sense == "up":
                    pos = coord - ov_len
                else:  # down
                    pos = coord + ov_len
            else:  # reverse strand
                if sense == "up":
                    pos = coord + ov_len
                else:  # down
                    pos = coord - ov_len

            if pos < 0 or pos >= len(seq):
                all_overlap = False
                break

            base_coord = int(base_mem[ref_point]) - 1
            if direction == "+":
                base_pos = base_coord - ov_len if sense == "up" else base_coord + ov_len
            else:
                base_pos = base_coord + ov_len if sense == "up" else base_coord - ov_len

            if base_pos < 0 or base_pos >= len(base_genome):
                all_overlap = False
                break

            if base_genome[base_pos] != seq[pos]:
                all_overlap = False
                break

        if not all_overlap:
            break
        ov_len += 1

    return ov_len


def find_overlaps(cluster_df, genomes):
    """
    For a cluster:
      - Detect shared overlaps upstream and downstream of start and stop codons
      - Merge up/down overlaps into single start and stop records
      - Extract overlap sequences and per-member coordinates

    Returns:
      pd.DataFrame with one row per region (start, stop)
    """

    genome_dict = {rec.id: rec.seq for rec in genomes}
    cluster_id = cluster_df["cluster"].iloc[0]
    direction = cluster_df["frame"].astype(str).iloc[0]

    # --- Step 1: Compute all four overlap lengths ---
    ov_dict = {}
    for region in ["start", "stop"]:
        for sense in ["up", "down"]:
            key = f"{region}_{sense}"
            ov_dict[key] = shared_overlap(cluster_df, genomes, region, direction, sense)

    base_mem = cluster_df.iloc[0]
    base_genome = genome_dict[get_genome_key(base_mem["gene"])]

    merged_rows = []

    # --- Step 2: Merge upstream and downstream overlaps for each region ---
    for region in ["start", "stop"]:
        up_len = ov_dict.get(f"{region}_up", 0)
        down_len = ov_dict.get(f"{region}_down", 0)
        total_len = up_len + down_len

        if total_len == 0:
            continue

        coord = int(base_mem[region]) - 1

        # Compute overall start/stop positions relative to base member
        if direction == "+":
            ov_start = coord - up_len + 1
            ov_stop = coord + down_len - 1
        else:  # reverse strand
            ov_start = coord - down_len + 1
            ov_stop = coord + up_len - 1

        ov_start = max(0, ov_start)
        ov_stop = min(len(base_genome) - 1, ov_stop)
        ov_seq = base_genome[ov_start:ov_stop + 1]

        row_data = {
            "cluster": cluster_id,
            "region": region,
            "ov_seq": str(ov_seq),
        }

        # Record per-member coordinates for the merged overlap
        for _, mem in cluster_df.iterrows():
            seq = genome_dict[get_genome_key(mem["gene"])]
            coord_mem = int(mem[region]) - 1

            if direction == "+":
                mem_ov_start = coord_mem - up_len + 1
                mem_ov_stop = coord_mem + down_len - 1
            else:
                mem_ov_start = coord_mem - down_len + 1
                mem_ov_stop = coord_mem + up_len - 1

            mem_ov_start = max(0, mem_ov_start)
            mem_ov_stop = min(len(seq) - 1, mem_ov_stop)

            genome_key = get_genome_key(mem["gene"])
            row_data[f"{genome_key}_ov_start"] = mem_ov_start + 1
            row_data[f"{genome_key}_ov_stop"] = mem_ov_stop + 1

        merged_rows.append(row_data)

    return pd.DataFrame(merged_rows)


def get_min_cds_num(members):
    # Extract the smallest CDS number (e.g. 0001 from "MZ501107.1_CDS_0001")
    cds_nums = []
    for mem in members:
        try:
            num = int(mem.split("_CDS_")[-1])
            cds_nums.append(num)
        except ValueError:
            continue
    return min(cds_nums) if cds_nums else 999999

def main():
    parser = argparse.ArgumentParser(description="Identify overlaps in the start or stop codon regions of CDS clusters using pharokka and mmseqs output")
    parser.add_argument("--genomes", required=True, help="Combined genome FASTA file used for pharokka input")
    parser.add_argument("--cds", required=True, help="cds_final_merged.tsv from pharokka output")
    parser.add_argument("--clusters", required=True, help="MMseqs2 clusters.tsv file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--cluster_threshold", type=int, required=False, help="Minimum number of members for a cluster to be valid. Default is the number of input genomes.")
    args = parser.parse_args()

    # --- Load genomes ---
    genomes = [rec for rec in SeqIO.parse(args.genomes, "fasta")]

    # --- Set cluster threshold n ---
    n = args.cluster_threshold if args.cluster_threshold is not None else len(genomes)

    # --- Load CDS annotations from pharroka ---
    cds_df = pd.read_csv(args.cds, sep="\t", dtype=str)

    # --- Load clusters from mmseqs ---
    clusters = defaultdict(set)
    with open(args.clusters) as f:
        for line in f:
            rep, mem = line.strip().split("\t")
            clusters[rep].add(mem)
    
    # --- Remove clusters with fewer than n members ---
    filtered_clusters = {rep: list(mems) for rep, mems in clusters.items() if len(mems) >= n}


    # Sort clusters by their smallest CDS number
    sorted_clusters = sorted(filtered_clusters.values(), key=get_min_cds_num)

    # Reassign cluster names in order
    clusters = {f"cluster{i+1}": mems for i, mems in enumerate(sorted_clusters)}

    # --- Make dfs for clusters ---
    cluster_dfs = []
    for cluster, mems in clusters.items():
        cluster_df = cds_df[cds_df["gene"].isin(mems)].copy()
        cluster_df["cluster"] = cluster
        cluster_dfs.append(cluster_df)


    
    # --- Find and record overlaps for each cluster ---
    all_overlaps = []

    for cluster_df in cluster_dfs:
        overlap_df = find_overlaps(cluster_df, genomes)
        all_overlaps.append(overlap_df)

    # --- Save a file describing which CDSs are in each cluster ---
    cluster_membership = []
    for cluster_df in cluster_dfs:
        cluster_id = cluster_df["cluster"].iloc[0]
        for gene in cluster_df["gene"]:
            cluster_membership.append({"cluster": cluster_id, "gene": gene})

    cluster_membership_df = pd.DataFrame(cluster_membership)
    cluster_membership_df.to_csv(f"{args.outdir}/cluster_membership.tsv", sep="\t", index=False)
    print(f"[✓] Saved cluster membership file to {args.outdir}/cluster_membership.tsv")


    final_overlap_df = pd.concat(all_overlaps, ignore_index=True)
    final_overlap_df.to_csv(f"{args.outdir}/cluster_overlap_summary.tsv", sep="\t", index=False)
        


if __name__ == "__main__":
    main()







#python Scripts/find_cluster_overlaps.py --genomes ../../fastas/bas14-15-16-18.fasta --cds results_bas14-15-16-18/pharokka_results/pharokka_cds_final_merged_output.tsv --clusters results_bas14-15-16-18/mmseqs/clusters.tsv --outdir results_bas14-15-16-18/
