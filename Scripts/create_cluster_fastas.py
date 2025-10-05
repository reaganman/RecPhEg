#!/usr/bin/env python3
import argparse
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import warnings


def main():
    parser = argparse.ArgumentParser(description="Generate FASTAs of nucleotide sequences for protein clusters from pharokka and mmseqs output.")
    parser.add_argument("--genomes", required=True, help="Combined genome FASTA file used for pharokka input")
    parser.add_argument("--cds", required=True, help="cds_final_merged.tsv from pharokka output")
    parser.add_argument("--clusters", required=True, help="MMseqs2 clusters.tsv file")
    parser.add_argument("--outdir", required=True, help="Output directory for cluster FASTAs")
    parser.add_argument("--cluster_threshold", type=int, required=False, help="Minimum number of members for a cluster to be valid. Default is the number of input genomes.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Load genome sequences ---
    genome_seqs = {rec.id: rec.seq for rec in SeqIO.parse(args.genomes, "fasta")}

    # --- Set cluster threshold n ---
    n = args.cluster_threshold if args.cluster_threshold is not None else len(genome_seqs)

    # --- Load CDS annotations ---
    df = pd.read_csv(args.cds, sep="\t", dtype=str)

    required = {"gene", "start", "stop", "contig", "frame"}
    if not required.issubset(df.columns):
        raise ValueError(f"CDS file must contain columns: {', '.join(required)}\nCheck pharroka output")

    cds_dict = {}
    header_dict = {}
    for _, row in df.iterrows():
        gene = row["gene"]
        start, stop = int(row["start"]), int(row["stop"])
        contig = row["contig"]

        # ensure start < stop
        start, stop = sorted([start, stop])

        seq = genome_seqs[contig][start-1:stop]  # 1-based to 0-based

        cds_dict[gene] = seq
        frame = row["frame"]
        header_dict[gene] = f"{contig}|{row['start']}-{row['stop']}|{frame}|{gene}"

    # --- Load clusters from mmseqs ---
    clusters = defaultdict(set)
    with open(args.clusters) as f:
        for line in f:
            rep, mem = line.strip().split("\t")
            clusters[rep].add(mem)

    # --- Check if representatives are missing ---
    for rep, mems in clusters.items():
        if rep not in mems:
            warnings.warn(f"Representative {rep} not listed as a member in its own cluster. Adding it manually.")
            mems.add(rep)

    # --- Remove clusters with fewer than n members ---
    filtered_clusters = {rep: list(mems) for rep, mems in clusters.items() if len(mems) >= n}

    # --- Write FASTAs per cluster ---
    for i, (rep, members) in enumerate(filtered_clusters.items(), start=1):
        fname = outdir / f"cluster_{rep}.fasta"
        with open(fname, "w") as out:
            for mem in members:
                if mem in cds_dict:
                    seq = cds_dict[mem]
                    header = header_dict[mem]
                    out.write(f">{header}\n{seq}\n")


if __name__ == "__main__":
    main()
