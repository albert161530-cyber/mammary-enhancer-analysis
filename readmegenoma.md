# Comparative Genomics Analysis Suite: Primate Regulatory Evolution

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![C++](https://img.shields.io/badge/C++-Pybind11-orange.svg)](https://github.com/pybind/pybind11)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance computational pipeline designed to test hypotheses regarding primate regulatory evolution, specifically focusing on lineage-specific sequence changes, transcription factor (TF) binding site turnover (PWM/Consensus scoring), and epistasis in candidate enhancers across key loci (`AR`, `ESR1`, `ESR2`, `CYP19A1`).

---

## 🧬 Overview

This repository contains the complete program used to evaluate how non-coding genomic differences between humans and ancestral/extant great apes (chimpanzees, gorillas, orangutans, and bonobos) potentially impact hormone receptor binding and regulatory dynamics (e.g., Androgen Receptor [AR], Estrogen Receptor Alpha/Beta [ESR1/ESR2], and shared Estrogen Response Elements [ERE]).

To handle intensive sliding-window Position Weight Matrix (PWM) scans, reverse-complements, and combinatorial motif interaction checks across large genomic windows, the engine integrates a **native C++ kernel module** compiled dynamically via `pybind11`.

---

## 🚀 Key Features

* **Hybrid C++ / Python Architecture:** Computationally heavy inner loops (PWM scoring, centered motif searches, reverse complements) are executed in optimized C++ for maximum throughput.
* **Multi-Primate Alignment Parsing:** Automatically processes Multiple Alignment Format (`.maf`) blocks to infer ancestral states and identify candidate single-nucleotide variants (SNPs) and structural insertions/deletions.
* **Enhancer Annotation Integration:** Intersects variants with SCREEN Registry V4 cCREs (ELS) and EpiMap ChromHMM cell state tracks (focusing on mammary epithelial/myoepithelial cells and general adipose tissues).
* **TF Binding & Consensus Scoring:** Evaluates Position Weight Matrices from JASPAR (2026 CORE vertebrate non-redundant set) alongside structured consensus motifs (ARE, ERE).
* **Non-Additive Interaction Scanning:** Scans for potential combinatorial effects (proxy epistasis) between closely spaced variant pairs within regulatory windows.
* **Automated Data Management:** Verifies, downloads, and caches required annotation files and sequences from UCSC, ENCODE, and EpiMap repositories.

---

## 🛠️ Installation & Prerequisites

### 1. Dependencies
Ensure you have Python 3.8+ installed along with the following packages:
```bash
pip install requests certifi
```

### 2. C++ Native Extension (`genoma_kernels`)
The script automatically attempts to compile the inline C++ code via `pybind11`. Ensure you have a C++ compiler (like `g++` or `clang`) and `pybind11` installed if you wish to build the acceleration module:
```bash
pip install pybind11 setuptools
```

---

## 📂 Repository Structure

* `Genoma1.py`: Main execution script containing the alignment parser, annotation loader, C++ kernel wrapper, motif scorer, and summary statistics generator.
* `README.md`: Project documentation and usage guide.
* *Generated Outputs (Created upon running)*:
  * `variantes_{gene}_anotadas.csv`: Annotated SNPs with regulatory overlaps and motif effect scores.
  * `efectos_motivo_{gene}.csv`: Detailed transcription factor motif disruption/gain predictions.
  * `deleciones_{gene}_anotadas.csv`: Candidate indel gaps and reconstructed ancestral sequences.
  * `interacciones_motivos_{gene}.csv`: Non-additive paired variant evaluations.
  * `estado_anotaciones_{gene}.csv`: Download and coverage status for genomic tracks.

---

## ⚙️ Usage

Run the main analysis script directly from your terminal:

```bash
python Genoma1.py
```

Upon execution, the program will:
1. Check for local `.maf` alignment files (and download them if configured/available).
2. Download and index required regulatory BED files (SCREEN V4, EpiMap).
3. Query the UCSC Genome Browser API (`hg38`) for target gene regions.
4. Execute high-performance C++ kernel scans for TF binding affinity shifts.
5. Export structured CSV reports for downstream statistical analysis and visualization.

---

## 🔍 Disclaimer & Context

This software is an **exploratory computational research tool**. Sequence-based PWM scoring, inferred ancestral consensus states, and proximity overlaps serve as molecular proxies and *in silico* hypotheses. They do not substitute for functional validation assays (such as ChIP-seq, luciferase reporter assays, or CRISPR editing) to confirm *in vivo* transcription factor binding or transcriptional activation.

---

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.