import argparse
import pandas as pd
import itertools
from Bio import SeqIO, pairwise2
import os



def percent_id(seq1, seq2, out_path=None, seq1_name="seq1", seq2_name="seq2"):
    """
    Align seq1 and seq2 and return the percent identity.
    If out_path is provided, save the alignment in FASTA format.
    """
    # Perform global alignment (Needleman–Wunsch)
    aligned_seqs = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    seq1_aln, seq2_aln, score, start, end = aligned_seqs[0]

    # Sanity check
    if len(seq1_aln) != len(seq2_aln):
        raise ValueError("Aligned sequences are different lengths!")

    # Calculate percent identity
    matches = sum(a == b for a, b in zip(seq1_aln, seq2_aln))
    pid = matches / len(seq1_aln) * 100

    # Optionally write alignment in FASTA format
    if out_path:
        with open(out_path, "w") as f:
            f.write(f">{seq1_name}\n")
            f.write(f"{seq1_aln}\n")
            f.write(f">{seq2_name}\n")
            f.write(f"{seq2_aln}\n")
            f.write(f"# Percent Identity: {pid:.2f}%\n")
            f.write(f"# Alignment Score: {score}\n")

    return pid
    




def main():
    parser = argparse.ArgumentParser(description="Make pairwise alignments and calculate percent identity for sequences in input fastas")
    parser.add_argument("--fasta_dir", required=True, help="Directory with fragment fastas")
    parser.add_argument("--csv_outdir", required=True, help="Path to save csv with alignment metrics")
    parser.add_argument("--aln_outdir", required=False, help="Directory to save pairwise fragment alignments")

    args = parser.parse_args()


    # Create output directories if needed
    os.makedirs(args.csv_outdir, exist_ok=True)
    if args.aln_outdir:
        os.makedirs(args.aln_outdir, exist_ok=True)

    # Loop over all FASTA files in the input directory
    for fasta_file in os.listdir(args.fasta_dir):
        if not fasta_file.endswith(".fasta") and not fasta_file.endswith(".fa"):
            continue

        fasta_path = os.path.join(args.fasta_dir, fasta_file)
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            print(f"[!] Skipping empty FASTA: {fasta_file}")
            continue

        names = [rec.id for rec in records]
        n = len(records)
        pid_matrix = pd.DataFrame(0.0, index=names, columns=names)

        # Compute all pairwise alignments
        for rec1, rec2 in itertools.combinations(records, 2):
            seq1_acc, seq2_acc = rec1.id, rec2.id
            seq1, seq2 = str(rec1.seq), str(rec2.seq)

            # Set output path for alignment if required
            out_path = None
            if args.aln_outdir:
                # Example: alignment/<fasta_name>_seq1_vs_seq2.fasta
                fasta_base = os.path.splitext(fasta_file)[0]
                aln_subdir = os.path.join(args.aln_outdir, fasta_base)
                os.makedirs(aln_subdir, exist_ok=True)
                out_path = os.path.join(aln_subdir, f"{seq1_acc}_vs_{seq2_acc}.fasta")

            pid = percent_id(seq1, seq2, out_path, seq1_acc, seq2_acc)

            # Fill matrix symmetrically
            pid_matrix.loc[seq1_acc, seq2_acc] = pid
            pid_matrix.loc[seq2_acc, seq1_acc] = pid

        # Diagonal = 100% self-identity
        for name in names:
            pid_matrix.loc[name, name] = 100.0

        # Save CSV matrix for this FASTA
        csv_outpath = os.path.join(args.csv_outdir, f"{os.path.splitext(fasta_file)[0]}_pid_matrix.csv")
        pid_matrix.to_csv(csv_outpath)
        print(f"[✓] Saved PID matrix for {fasta_file} → {csv_outpath}")

if __name__ == "__main__":
    main()
    





