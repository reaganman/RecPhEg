#!/usr/bin/env python3
import argparse
import pandas as pd
from itertools import combinations
from Bio import SeqIO
import ast
import os




def group_CDSs(CDDs, n):
    """
    Return all combinations of CDSs elements (genome accessions) 
    with size >= n as a list of lists.
    """
    CDS_combos = []
    for r in range(n, len(CDDs) + 1):
        CDS_combos.extend([list(c) for c in combinations(CDDs, r)])
    return CDS_combos



def main():
    parser = argparse.ArgumentParser(
        description="Update lovis4u annotation table to color code CDSs where overlaps were found"
    )
    parser.add_argument("--genomes", required=True, help="original fasta input with all genomes")
    parser.add_argument("--cluster_threshold", type=int, required=False,
                        help="Minimum number of members for a cluster to be valid. Default is the number of input genomes.")
    parser.add_argument("--overlaps", required=True, help="shared_kmers.csv file with start/stop overlaps")
    parser.add_argument("--annotations", required=True, help="feature_annotations_table.tsv from lovis4u output")
    parser.add_argument("--out", required=False, help="output prefix (default: overwrite annotations file)")
    args = parser.parse_args()


    def assign_color(feature_id):
            if feature_id not in CDS_overlap_regions:
                return "#cccccc"  # grey
            regions = CDS_overlap_regions[feature_id]
            if "start" in regions and "stop" in regions:
                return "#d367eb"  # purple
            elif "start" in regions:
                return "#15ed2b"  # green
            elif "stop" in regions:
                return "#ed1515"  # red
            return "#cccccc"  # fallback grey

    def assign_cluster(feature_id):
        name = annotations_df.loc[
            annotations_df["feature_id"] == feature_id, "name"
        ].iloc[0]
        if name == "hypothetical protein":
            name = "HP"
        if feature_id in CDS_clusters:
            name = "; ".join([name, CDS_clusters[feature_id]])
        return name


    def filter_annotations(group_accessions, annotations_df):
        # make a new df with each row where annotations_df[locus_id] is in the group_accessions list
        return annotations_df[annotations_df["locus_id"].isin(group_accessions)].copy()

    # make different combos of input genomes to visualize overlaps for
    n = args.cluster_threshold or len(list(SeqIO.parse(args.genomes, "fasta")))
    genome_accessions = [rec.id for rec in SeqIO.parse(args.genomes, "fasta")]
    overlaps_df = pd.read_csv(args.overlaps)

    combo_groups = {}           # dict to track accessions in different combo groups

    # ensure accessions column is a list (in case it was saved as a string)
    if overlaps_df["accessions"].dtype == "object":
        try:
            overlaps_df["accessions"] = overlaps_df["accessions"].apply(ast.literal_eval)
        except Exception:
            pass  # already list-like


    print(f"combo groups: {group_CDSs(genome_accessions, n)}")
    for i, combo in enumerate(group_CDSs(genome_accessions, n)):
        combo_overlaps_df = overlaps_df[overlaps_df["accessions"].apply(lambda x: set(combo).issubset(set(x)))]   

        combo_group = f"group{i}"
        combo_groups[combo_group] = [combo]

        CDS_overlap_regions = {}
        CDS_clusters = {}
        overlap_accessions = []  # initialize empty list here ✅

        for _, row in combo_overlaps_df.iterrows():
            cluster = "_".join(row["Cluster"].split("_")[0:4])
            region = row["Cluster"].split("_")[-1]
            CDSs = row["CDSs"].split(";")
            overlap_accessions.extend([CDS.split("_")[0] for CDS in CDSs])

            for CDS in CDSs:
                if CDS in CDS_clusters and CDS_clusters[CDS] != cluster:
                    raise ValueError(
                        f"CDS: {CDS} already in cluster: {CDS_clusters[CDS]}\nCannot add to {cluster}"
                    )
                CDS_clusters[CDS] = cluster

                CDS_overlap_regions.setdefault(CDS, set()).add(region)

        # === Only filter if we actually had overlap accessions ===
        annotations_df = pd.read_csv(args.annotations, sep="\t", dtype=str)
        if overlap_accessions:
            print(f"filtering annotation table for group: {combo}")
            group_annotations_df = filter_annotations(combo, annotations_df)
            print(f"filtered annotation table:\n{group_annotations_df}")

            group_annotations_df["fill_colour"] = group_annotations_df["feature_id"].apply(assign_color)
            group_annotations_df["name"] = group_annotations_df["feature_id"].apply(assign_cluster)

            out_prefix = args.out if args.out else args.annotations
            out_prefix = os.path.splitext(out_prefix)[0]
            out_file = f"{out_prefix}_{combo_group}.tsv"

            group_annotations_df.to_csv(out_file, sep="\t", index=False)
            print(f"Saved updated annotations for combo {combo} → {out_file}")
        else:
            print(f"⚠️ Skipping {combo} (no overlaps found)")

    combo_group_df = pd.DataFrame(list(combo_groups.items()), columns=["group", "accessions"])
    combo_group_df.to_csv(f"{out_prefix}_combo_groups.csv", index=False)



if __name__ == "__main__":
    main()


# ex usage: 
# python update_annotations_table.py --overlaps results_bas14-18/shared_kmers.csv --annotations results_bas14-18/ilovis4u/feature_annotation_table.tsv --out results_bas14-18/ilovis4u/feature_annotation_table_updated.tsv
# python Scripts/update_annotations_table.py --overlaps results_bas14-18/shared_kmers.csv --annotations results_bas14-18/lovis4u/feature_annotation_table.tsv --out results_bas14-18/lovis4u/feature_annotation_table.tsv --genomes ../../fastas/bas14-18.fasta --cluster_threshold 3