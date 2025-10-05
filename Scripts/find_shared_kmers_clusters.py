#!/usr/bin/env python3
import argparse
import csv
import os
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_cluster_kmers(seqs, k):
    """
    Generate all subsequences of length k for each seq in a cluster
    seqs = list of cluster sequences [(ID1, seq1), (ID2, seq2)...]
    Returns dictionary with {kmer_seq:(CDS1, CDS2,...)} where 'CDS* are the sequences in the cluster that have the kmer
    """

    kmers_dict = {}

    for id, seq in seqs:
        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        cds_id = id.split('|')[3]

        for kmer in kmers:
            if kmer not in kmers_dict:
                kmers_dict[kmer] = [cds_id]
            else:
                kmers_dict[kmer].append(cds_id)

    return kmers_dict



def main():
    parser = argparse.ArgumentParser(description="Find shared subsequences between sequences in all fasta files found in input directory")
    parser.add_argument("cluster_dir", help="directory containing fasta files of sequence clusters")
    parser.add_argument("-k", type=int, required=True, help="Subsequence length (k-mer size)")
    parser.add_argument("-n", type=int, required=False, help="Minimum number of cluster members with shared kmer to report overlap. Default is only kmers shared by all members")
    parser.add_argument("-o", default="shared_kmers.csv", help="Output CSV file (default: shared_kmers.csv)")
    args = parser.parse_args()

    with open(args.o, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Cluster", "n", "accessions", "CDSs", "Sequence"])

        # Loop through each cluster fasta
        for filename in os.listdir(args.cluster_dir):
            if filename.endswith(".fasta") and os.path.isfile(os.path.join(args.cluster_dir, filename)):
                cluster_name = os.path.splitext(filename)[0]
                filepath = os.path.join(args.cluster_dir, filename)

                # Load cluster sequences
                sequences = [(record.id, str(record.seq)) for record in SeqIO.parse(filepath, "fasta")]

                # Default n = all sequences in cluster
                n = args.n if args.n is not None else len(sequences)

                #print(f"Searching for shared {args.k}-mers among {len(sequences)} sequences in {filename}")
                #print(f"kmers found in {n} or more sequences in {filename} will be saved to {args.o}\n")

                # Make dictionary with kmers
                cluster_kmers = get_cluster_kmers(sequences, args.k)

                # Count kmer occurrences and write those >= n to output
                for kmer_seq, CDSs in cluster_kmers.items():
                    #print(f"Cluster members: {CDSs}\nShared K-mer: {kmer_seq}")
                    CDS_accessions = [CDS.split('_')[0] for CDS in CDSs]
                    if len(CDSs) >= n:
                        result = [cluster_name, len(CDSs), CDS_accessions, ";".join(CDSs), kmer_seq]
                        print(f"Writing {result} to output")
                        writer.writerow(result)


if __name__ == "__main__":
    main()
