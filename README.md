# 🧬 RecPhEg: Recombinatorial Phage Engineering

**RecPhEg** is an **interactive pipeline** for designing **Re**combinatorial **Ph**age **E**n**g**ineering assemblies.

**NOTE: THIS PIPELINE IS IN DEVELOPMENT AND HAS NOT BEEN TESTED IN-VITRO**

---

## 🚀 Overview

RecPhEg can be used for designing interchangeable genome fragments from related phage genomes.  

The pipeline currently requires **interactive input** and operates on a **combined FASTA file** containing the complete genome sequences of related phages.

---

## 🧩 Workflow

The RecPhEg workflow is as follows:

1. **Call Coding Sequences (CDSs)**  
   Identify coding regions across all input genomes.

2. **Cluster Homologous CDSs**  
   Group homologous genes to identify shared and unique genomic elements.

3. **Identify Candidate Overlaps**  
   Detect potential overlaps within CDS clusters.  
   Overlaps must include a **start** or **stop codon**.

4. **Filter Candidate Overlaps**  
   Apply sequence and positional filters to remove unsuitable overlaps.

5. **Determine Optimal Overlap Combinations**  
   Compute the most suitable overlap set for a specified number of fragments (**N**) and maximum fragment size (**max_frag_size**).

6. **Build and Analyze Fragments**  
   Construct fragment assemblies based on optimal overlap sets.  
   Design and evaluate **PCR primers** for fragment synthesis and assembly.

---

## 📥 Input Requirements

- A **combined FASTA file** containing complete genome sequences of related phages.

---

## ⚙️ Output

- Optimal fragment assemblies  
- Sequence and structural analysis of designed fragments  
- PCR primer designs for downstream validation

---

## 🧠 Notes

- The pipeline is **interactive** (automation features are under development).  
- Intended for **comparative phage genomics** and **rational recombination design**.

---

## 📅 Future Work

- Full automation of the interactive steps  
- Integration of codon optimization and assembly simulation modules  
- Expansion to support metagenomic phage datasets
