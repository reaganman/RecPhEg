#!/usr/bin/env python3
import argparse
from collections import defaultdict
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import warnings


def main():
    parser = argparse.ArgumentParser(description="Generate FASTAs of nucleotide start/stop regions for protein clusters from pharokka and mmseqs output.")
    parser.add_argument("--genomes", required=True, help="Combined genome FASTA file used for pharokka input")
    parser.add_argument("--cds", required=True, help="cds_final_merged.tsv from pharokka output")
    parser.add_argument("--clusters", required=True, help="MMseqs2 clusters.tsv file")
    parser.add_argument("--outdir", required=True, help="Output directory for cluster FASTAs")
    parser.add_argument("--flank", type=int, required=True, help="Number of bp to include on each side of start/stop codon")
    parser.add_argument("--cluster_threshold", type=int, required=False,
                        help="Minimum number of members for a cluster to be valid. Default is the number of input genomes.")
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
        raise ValueError(f"CDS file must contain columns: {', '.join(required)}")

    cds_info = {}
    for _, row in df.iterrows():
        gene = row["gene"]
        start, stop = int(row["start"]), int(row["stop"])
        contig = row["contig"]
        frame = row["frame"]
        cds_info[gene] = (contig, start, stop, frame)

    # --- Load clusters ---
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

    # --- Filter clusters by threshold ---
    filtered_clusters = {rep: list(mems) for rep, mems in clusters.items() if len(mems) >= n}

    # --- Extract regions ---
    flank = args.flank - 1
    for i, (rep, members) in enumerate(filtered_clusters.items(), start=1):
        start_fasta = outdir / f"cluster_{rep}_start.fasta"
        stop_fasta = outdir / f"cluster_{rep}_stop.fasta"

        with open(start_fasta, "w") as fout_start, open(stop_fasta, "w") as fout_stop:
            for mem in members:
                if mem not in cds_info:
                    continue
                contig, start, stop, frame = cds_info[mem]
                seq = genome_seqs[contig]

                # --- start codon region ---
                if frame == "+":
                    left = max(1, start + 2 - flank)  # include start codon
                    right = start + flank
                elif frame == "-":
                    left = start - flank
                    right = min(len(seq), start - 2 + flank)
                else:
                    raise ValueError("Could not determine context region boundaries")

                start_seq = seq[left-1:right]  # convert to 0-based
                fout_start.write(f">{contig}|{left}-{right}|{frame}|{mem}|start_region\n{start_seq}\n")

                # --- stop codon region ---
                if frame == "+":
                    left = stop - flank
                    right = min(len(seq), stop - 2 + flank)
                elif frame == "-":
                    left = max(1, stop + 2 - flank)
                    right = stop + flank
                else:
                    raise ValueError("Could not determine context region boundaries")

                stop_seq = seq[left-1:right]  # convert to 0-based
                fout_stop.write(f">{contig}|{left}-{right}|{frame}|{mem}|stop_region\n{stop_seq}\n")


if __name__ == "__main__":
    main()
