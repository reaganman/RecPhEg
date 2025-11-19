import argparse
import pandas as pd
import re
from Bio.Seq import Seq
import os



def check_overlap_gc(ov_seq, min_size=20):
    """
    Check that the GC content of ov_seq is between 40–60%. If yes, return (sequence, GC%).
    If not, search for subsequences >= min_size within ov_seq that meet the GC range
    and return the longest one and its GC%. If none found, return (None, None).
    """
    if not ov_seq:
        return None, None

    seq = ov_seq.upper()

    def gc_content(s):
        return sum(base in {"G", "C"} for base in s) / len(s)

    # Check full sequence
    gc = gc_content(seq)
    if 0.4 <= gc <= 0.6:
        return seq, round(gc * 100, 2)

    # Otherwise, search for subsequences
    best_subseq = None
    best_gc = 0

    for i in range(len(seq)):
        for j in range(i + min_size, len(seq) + 1):
            subseq = seq[i:j]
            gc_sub = gc_content(subseq)
            if 0.4 <= gc_sub <= 0.6:
                # prefer longer subsequences, then higher GC
                if best_subseq is None or len(subseq) > len(best_subseq) or (
                    len(subseq) == len(best_subseq) and gc_sub > best_gc
                ):
                    best_subseq = subseq
                    best_gc = gc_sub

    if best_subseq:
        return best_subseq, round(best_gc * 100, 2)
    else:
        return None, None

def find_homopolymers(ov_seq, max_homo=4):
    """
    Return a list of tuples [(start, stop, sequence), ...] for each homopolymer >= max_homo.
    Coordinates are 0-based, and 'stop' is exclusive (Python-style).
    """
    ov_seq = ov_seq.upper()
    homos = []
    i = 0
    n = len(ov_seq)

    while i < n:
        j = i + 1
        while j < n and ov_seq[j] == ov_seq[i]:
            j += 1
        # If this run is long enough, record it
        if j - i >= max_homo:
            homos.append((i, j, ov_seq[i:j]))
        i = j  # move to next segment

    return homos

        
def recommend_no_homo(ov_seq, homos, max_homo_size=3, min_len=1):
    """
    Produce fragments for every continuous subsequence that does NOT contain
    a homopolymer longer than max_homo_size. Fragments *keep* up to
    `max_homo_size` bases of each trimmed homopolymer on both sides, so
    adjacent fragments will overlap by `max_homo_size`.

    Parameters
    ----------
    ov_seq : str
        Original overlap sequence.
    homos : list of (start, stop, seq)
        Output from find_homopolymers(); start is 0-based, stop is exclusive.
    max_homo_size : int
        Maximum allowed homopolymer length to keep on each side.
    min_len : int
        Minimum fragment length to include in the result.

    Returns
    -------
    List[str]
        Fragments (may overlap) each guaranteed to contain no homopolymer >
        max_homo_size and to be at least min_len long.
    """
    ov_seq = ov_seq.upper()
    if not homos:
        return [ov_seq] if len(ov_seq) >= min_len else []

    fragments = []
    last_cut = 0

    # assume homos are sorted by start; if not, sort them
    homos_sorted = sorted(homos, key=lambda h: h[0])

    for start, stop, seq in homos_sorted:
        run_len = stop - start
        if run_len > max_homo_size:
            # left fragment: from last_cut up to (and including) max_homo_size bases
            # from the start of this homopolymer
            left_frag = ov_seq[last_cut : start + max_homo_size]
            if len(left_frag) >= min_len:
                fragments.append(left_frag)

            # next fragment should begin max_homo_size bases before the end of the run,
            # i.e. we keep max_homo_size bases of the homopolymer at the start of the next fragment
            last_cut = max(0, stop - max_homo_size)

    # final trailing fragment
    if len(ov_seq) - last_cut >= min_len:
        fragments.append(ov_seq[last_cut:])

    return fragments


def find_repeats(seq, motif_len=2, min_repeats=3):
    """
    Find tandem repeats of motifs of length motif_len repeated >= min_repeats times.
    Returns [(start, stop, motif, repeat_count, sequence), ...]
    """
    seq = seq.upper()
    pattern = re.compile(rf'((\w{{{motif_len}}})\2{{{min_repeats-1},}})')
    results = []
    for match in pattern.finditer(seq):
        start, stop = match.span()
        full = match.group(1)
        motif = match.group(2)
        repeat_count = len(full) // motif_len
        results.append((start, stop, motif, repeat_count, full))
    return results


def find_palindromes(seq, min_len=4, max_len=10, min_size=20):
    """
    Find palindromic subsequences within a DNA sequence, filter out nested ones,
    and return both:
      1) a list of (start, stop, seq) palindromes
      2) the longest non-palindromic subsequence (>= min_size), if any

    Parameters
    ----------
    seq : str
        DNA sequence
    min_len : int
        Minimum palindrome length to detect
    max_len : int
        Maximum palindrome length to detect
    min_size : int
        Minimum length of the returned non-palindromic region

    Returns
    -------
    filtered_palindromes : list of (start, stop, str)
    longest_non_pal_seq : str or None
    """
    seq = seq.upper()
    palindromes = []

    # --- find all palindromic subsequences ---
    for i in range(len(seq)):
        for l in range(min_len, max_len + 1):
            frag = seq[i:i + l]
            if len(frag) == l and frag == str(Seq(frag).reverse_complement()):
                palindromes.append((i, i + l, frag))

    if not palindromes:
        # no palindromes at all
        if len(seq) >= min_size:
            return [], seq
        else:
            return [], None

    # --- sort and filter nested ones ---
    palindromes.sort(key=lambda x: (x[0], -(x[1] - x[0])))
    filtered = []
    current_end = -1

    for start, stop, frag in palindromes:
        if start >= current_end:
            filtered.append((start, stop, frag))
            current_end = stop

    # --- find non-palindromic fragments ---
    non_pal_regions = []
    prev_end = 0
    for start, stop, frag in filtered:
        # region before this palindrome
        if start > prev_end:
            non_pal_regions.append(seq[prev_end:start])
        prev_end = stop
    # tail after last palindrome
    if prev_end < len(seq):
        non_pal_regions.append(seq[prev_end:])

    # --- select the longest valid fragment ---
    longest_non_pal = None
    if non_pal_regions:
        longest_non_pal = max(
            (frag for frag in non_pal_regions if len(frag) >= min_size),
            key=len,
            default=None
        )

    return filtered, longest_non_pal

def optimize_overlap(ov_seq, min_size=20, max_homo_size=4):
    """
    Returns the best overlapping subsequnce without mono- di- or inverted repeats, GC content 40-60%
    and size > min_size 
    """

    best_ov = ov_seq

    # Check for homopolymers
    homos = find_homopolymers(ov_seq)
    if homos:
        print(f"{len(homos)} homopolymers found: {homos}")
        # Check if sequences <= min_size exists without homopolymers
        no_homo = recommend_no_homo(ov_seq, homos, min_len=min_size)
        if no_homo:
            # Recommend largest ov_seq without homopolymers
            max_no_homo = max(no_homo, key=len)
            best_ov = max_no_homo
            print(f"Recommended overlap without homopolymers >= {max_homo_size}:\n{max_no_homo}  Length: {len(max_no_homo)}")
        else:
            print(f"No valid overlap <= {min_size} without homopolymers >= {max_homo_size} :(")
            return None

    # Check for dinucleotide repeats
    di_repeats = find_repeats(best_ov, 2, 3)
    if di_repeats:
        print(f"Dinucleotide repeat region(s) found in {best_ov}\n{di_repeats}")
        return None

    # Check for inverted repeats 
    iv_repeats = find_palindromes(best_ov, 6, 10, min_size)
    if iv_repeats[0]:
        print(f"Inverted repeat(s) found in {best_ov}\n{iv_repeats}")
        if iv_repeats[1]:
            print(f"Recomended overlap witout inverted repeat(s) {iv_repeats[1]}\t Length: {len(iv_repeats[1])}")
        else:
            print(f"No valid overlap without inverted repeaats :(")
            return None
    best_ov = iv_repeats[1]


    # Check GC%
    best_ov = check_overlap_gc(best_ov, min_size)
    if best_ov[0]:
        print(f"Optimized overlap for {ov_seq}\n{best_ov[0]}\tGC%: {best_ov[1]}")
    else:
        print(f"Could not find optimized overlap (40-60% GC) for {ov_seq} :(\n")
        return None

    if len(best_ov[0]) < min_size:
        print(f"Overlap below minimum size ({len(best_ov[0])} < {min_size}) — skipping.")
        return None
    return best_ov[0]





def main():
    parser = argparse.ArgumentParser(description="Select optimal overlap sequences from those identified by find_cluster_overlaps.py")
    parser.add_argument("--overlaps", required=True, help="cluster_overlaps.tsv from find_cluster_overlaps.py")
    parser.add_argument("--outdir", required=False, help="directory to save filetered overlaps. Default is same dir as overlaps input.")
    parser.add_argument("--min_size", required=False, help="Minimum size for an overlap. Default is 20 bp.")
    
    args = parser.parse_args()

    # Determine output directory
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = os.path.dirname(os.path.abspath(args.overlaps)) or "."

    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "optimized_overlaps.tsv")

    # Determine min overlap size
    if args.min_size:
        min_size = int(args.min_size)
    else:
        min_size = 20 

    # Load original overlaps
    ov_df = pd.read_csv(args.overlaps, sep="\t")

    # Filter out those smaller that 20 bp
    ov_df = ov_df[ov_df["ov_seq"].str.len() >= min_size]

    
    optimized_overlaps_rows = []
    # Remove overlaps with repeats
    for i, row in ov_df.iterrows():
        ov = row['ov_seq']
        print(f"Analyzing: {ov}")
        optimized_ov = optimize_overlap(ov, min_size)
        # Make updated row for optimized overlap
        if optimized_ov:
            delta_start = ov.find(optimized_ov)
            delta_stop = delta_start + len(optimized_ov) - len(ov)
            print(f"Optimized overlap delta_start:\t{delta_start} delta stop:\t{delta_stop}\tSize:\t{len(optimized_ov)}\n")
            
            opt_ov_row = row.copy()
            opt_ov_row['ov_seq'] = optimized_ov

            # adjust all *_start and *_stop columns
            for col in ov_df.columns:
                if col.endswith("_start") and pd.notnull(row[col]):
                    opt_ov_row[col] = row[col] + delta_start
                elif col.endswith("_stop") and pd.notnull(row[col]):
                    opt_ov_row[col] = row[col] + delta_stop
            
            optimized_overlaps_rows.append(opt_ov_row)

    else:
        print(f"Could not optimize {ov}\n")

    # create new DataFrame and save
    optimized_ov_df = pd.DataFrame(optimized_overlaps_rows)
    optimized_ov_df.to_csv(outfile, sep="\t", index=False)
    print(f"\n[✓] Saved {len(optimized_ov_df)} optimized overlaps to: {outfile}")



    
        






if __name__ == "__main__":
    main()





    

