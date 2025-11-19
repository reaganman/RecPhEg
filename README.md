# 🧬 RecPhEg: Recombinatorial Phage Engineering

**RecPhEg** is an **interactive pipeline** for designing **Re**combinatorial **Ph**age **E**n**g**ineering assemblies.

**NOTE: THIS PIPELINE IS IN DEVELOPMENT AND HAS NOT BEEN TESTED IN-VITRO**

---

## 🚀 Overview

RecPhEg can be used for designing interchangeable genome fragments wigth shared overlaps from related phage genomes.  

The pipeline currently requires a **combined FASTA file** containing the complete genome sequences of related phages, and a minimum overlap length.

---

## 🧩 Workflow

The RecPhEg workflow is as follows:

1. **Call Coding Sequences (CDSs)**  
   Identify putative coding regions across all input genomes using pharokka and assign each to PHROG categrory with mmseqs2 where possible.

2. **Cluster Homologous CDSs**  
   Group homologous CDSs using mmseqs2 with **--cluster-mode 0 --cov-mode 1 -c 0.8 -s 7.5 --min-seq-id 0.25**

3. **Identify Candidate Overlaps**  
   Detect potential overlaps within CDS clusters.  
   Overlaps must include a **start** or **stop codon**.

4. **Filter/Optimize Candidate Overlaps**  
   Remove repetatve and low GC sequences and those smaller that the size threshold.

5. **Make output diagram**
   Color code CDSs with suitable overlaps containing start or stop codons

---

## 📥 Input 

- **Combined FASTA file** containing complete genome sequences of related phages.
- **Minimum Overlap size**

---

## ⚙️ Output

- Optimized overlap sequences  
- Diagrams showing overlaps and pharokka annotations
---

## 🧠 Extra Scripts

- **fetch_accessions.py** Create a FASTA file for each value in the accession column of input csv


---

## 📅 Future Work
- Test overlap selection stategy in-vitro
- Use blastn to identify interchangeable fragments across NCBI nt databases based on overlap selection
- Automate overlap selection? 

