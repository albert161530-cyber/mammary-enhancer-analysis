# Comparative Genomics Analysis Suite: Primate Regulatory Evolution

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![C++](https://img.shields.io/badge/C++-Pybind11-orange.svg)](https://github.com/pybind/pybind11)
[![AlphaGenome](https://img.shields.io/badge/Deep%20Learning-AlphaGenome-green.svg)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A high-performance, multi-stage computational pipeline designed to test hypotheses regarding primate cis-regulatory evolution—specifically evaluating whether the evolution of permanently enlarged human breasts was driven primarily by the loss of repressive Androgen Receptor (`AR_loss`) "brakes" or the gain of Estrogen Receptor Alpha (`ESR1_gain`) "accelerators" across key loci (`AR`, `ESR1`, `ESR2`, `CYP19A1`).

---

## 🧬 Overview

This repository contains the complete multi-class framework used to evaluate how non-coding genomic differences (SNVs, short indels, deletions, and human-specific insertions) between humans and ancestral/extant great apes (chimpanzees, gorillas, orangutans, and bonobos) impact hormone receptor binding, chromatin architecture, and transcriptional dynamics.

To handle intensive sliding-window Position Weight Matrix (PWM) scans, reverse-complements, and combinatorial motif interaction checks across 1 Mb genomic windows, the engine integrates a **native C++ kernel module** compiled dynamically via `pybind11`, coupled with **AlphaGenome** deep-learning sequence predictions across 279 cell and tissue tracks.

---

## 🚀 Key Features

* **Hybrid C++ / Python Architecture:** Computationally heavy inner loops (PWM scoring, centered motif searches, reverse complements) are executed in optimized C++ (`genoma_kernels`) for maximum throughput.
* **Multi-Primate Alignment Parsing:** Processes Multiple Alignment Format (`.maf`) blocks across 30 mammals to reconstruct ancestral primate states and identify candidate SNVs, deletions, and human-specific insertions (HSIs).
* **Cross-Species Epigenetic & Enhancer Integration:** Intersects variants with SCREEN Registry V4 cCREs (ELS) and EpiMap ChromHMM cell state tracks (focusing on mammary luminal/myoepithelial cells, breast fibroblasts, and adipose tissues), while scanning unannotated regions to capture ancestral enhancers deactivated in humans.
* **AlphaGenome Deep-Learning Evaluation:** Predicts single-variant and joint-haplotype impacts (`human` vs. `reconstructed ancestral`) across ATAC-seq, DNase-seq, Histone ChIP-seq, TF ChIP-seq, CAGE, PRO-cap, and RNA-seq tracks, isolating non-additive epistasis.
* **Endocrine & Physiological Modeling:** Isolates the net regulatory proportion of androgenic versus estrogenic activity—holding basal hormone levels and affinities equal—across physiological states (Tanner stages, menstrual cycle, pregnancy, lactation, gynecomastia, PCOS, and CAIS).

---

## 🛠️ Installation & Prerequisites

### 1. Python Dependencies
Ensure you have Python 3.8+ installed along with the required scientific, statistical, and AlphaGenome client packages:
```bash
pip install numpy pandas scipy matplotlib requests certifi alphagenome
```

### 2. C++ Native Extension (`genoma_kernels`)
The pipeline automatically compiles the inline C++ acceleration module via `pybind11`. Ensure you have a C++ compiler (`g++`, `clang`, or MSVC) installed:
```bash
pip install pybind11 setuptools
```

### 3. AlphaGenome API Key Configuration
For scripts that query the AlphaGenome deep-learning model (`automatizar_bloques.py` and `auditoria_inserciones.py`), provide your API key using either an environment variable or a local text file in the project root:
```bash
export ALPHAGENOME_API_KEY="your_api_key_here"
# Or save the key directly into a file named api_key.txt in the root directory
```

---

## 📂 Repository Structure

### Core Pipeline Scripts
1. `Genoma1.py`: **Stage 1 — Comparative MAF Parser & Motif Scanner.** Extracts human-specific substitutions and deletions from `.maf` alignments, intersects with SCREEN/EpiMap mammary/adipose beds, runs C++ PWM/consensus scans (JASPAR 2026), and exports consolidated candidate variants.
2. `automatizar_bloques.py`: **Stage 2 — AlphaGenome Joint Haplotype & Epistasis Engine (SNVs/Deletions).** Evaluates individual and joint haplotype effects (`AR_loss` vs. `ESR1_gain`) across 1 Mb windows and 279 tracks, quantifying non-additive regulatory interactions.
3. `auditoria_inserciones.py`: **Stage 3 — Human-Specific Insertions (HSIs) Auditor.** Detects and evaluates structural insertions across the 4 candidate loci using AlphaGenome across mammary and adipose tracks (supports smart auto-resume).
4. `analisis_posthoc_mamario.py`: **Stage 4 — Post-Hoc Statistical & Tissue-Specific Auditor.** Performs exact binomial and paired Wilcoxon tests, filters mammary/adipose compartments (`n=178`), ranks top driver variants, analyzes MCF-7 direct factor cross-talk, and generates publication-grade figures.
5. `generar_reporte_multiclase.py`: **Stage 5 — Multi-Class Evolutionary Consolidator.** Integrates SNV/indel and structural insertion outcomes into a unified evolutionary dominance table.
6. `modelo_modular_potenciadores.py`: **Stage 6 — Physiological & Hormonal Enhancer Simulator.** Simulates net estrogenic/androgenic regulatory balances across clinical and developmental states comparing wild-type human, wild-type chimpanzee, and "humanized-enhancer" chimpanzee models.

---

## ⚙️ Execution Order & Commands

To reproduce the entire analysis from raw alignments to final consolidated tables and figures, execute the scripts in the following sequential order:

### Step 1: Parse MAF Alignments, Annotate Enhancers & Score PWM Motifs
Runs the initial comparative genomic screening across `AR`, `ESR1`, `ESR2`, and `CYP19A1`, generating `variantes_todos_loci.csv` and baseline hormonal templates.
```bash
python Genoma1.py
```
* **Key Outputs:** `variantes_todos_loci.csv`, `variantes_{gene}_anotadas.csv`, `deleciones_todos_loci.csv`, `comparacion_modelos_3_escenarios.csv`, `resultados_modelo_hormonal.csv`.

---

### Step 2: Run AlphaGenome Predictions for Substitutions & Short Indels
Takes `variantes_todos_loci.csv` and queries AlphaGenome for single-variant and joint-haplotype effects (`AR_loss` vs. `ESR1_gain`) across 1 Mb windows.
```bash
python automatizar_bloques.py
```
* **Key Outputs:** `resultados_impacto_conjunto_comparacion.csv`, `resultados_impacto_conjunto_por_track.csv`, `resultados_impacto_individual_por_track.csv`, `resultados_impacto_conjunto_regiones.csv`.

---

### Step 3: Run AlphaGenome Predictions for Human-Specific Insertions (HSIs)
Evaluates lineage-specific insertions across the same 1 Mb windows and mammary/adipose tracks. *(Note: Features smart auto-resume if interrupted).*
```bash
python auditoria_inserciones.py
```
* **Key Outputs:** `variantes_inserciones_loci.csv`, `tabla_5_impacto_inserciones_mamario.csv`, `tabla_5_resumen_inserciones_por_modalidad.csv`.

---

### Step 4: Execute Post-Hoc Statistical Auditing & Generate Publication Figures
Processes the AlphaGenome outputs from Step 2 to generate tissue-stratified statistical tests (Binomial & Wilcoxon), locus-level epistasis tables, top driver rankings, and high-resolution plots.
```bash
python analisis_posthoc_mamario.py
```
* **Key Outputs:** 
  * `tabla_1_resumen_global_y_mamario_por_modalidad.csv`
  * `tabla_2_epistasis_no_aditividad_por_locus.csv`
  * `tabla_3_top_variantes_drivers_mamario_adiposo.csv`
  * `tabla_4_cross_talk_factores_directos_MCF7.csv`
  * `resumen_validacion_alphagenome_mamario.png` (300 DPI)

---

### Step 5: Consolidate Multi-Class Evolutionary Report (SNVs + Insertions)
Combines the results from Step 2 (`resultados_impacto_conjunto_comparacion.csv`) and Step 3 (`tabla_5_impacto_inserciones_mamario.csv`) into the final multi-class summary table.
```bash
python generar_reporte_multiclase.py
```
* **Key Output:** `tabla_6_consolidado_multiclase_evolutivo.csv`.

---

### Step 6: Run the Modular Physiological & Endocrine Simulation (Standalone/Complementary)
Evaluates how the derived enhancer activity weights (`AR` brake vs. `ER` accelerator) interact with systemic estradiol ($E_2$) and testosterone ($T$) levels across physiological scenarios (Tanner stages, pregnancy, lactation, gynecomastia).
```bash
python modelo_modular_potenciadores.py
```
* **Key Outputs:** `entradas_modelo_hormonal.csv` and console/CSV comparative readouts of Human vs. Chimpanzee vs. Humanized-Chimpanzee net tissue signaling.

---

### ⚡ Quick One-Liner (Full Pipeline Execution)
If your `.maf` files and `api_key.txt` are already in place, you can run the entire sequential pipeline in a single command:
```bash
python Genoma1.py && python automatizar_bloques.py && python auditoria_inserciones.py && python analisis_posthoc_mamario.py && python generar_reporte_multiclase.py && python modelo_modular_potenciadores.py
```

---

## 🔍 Disclaimer & Context

This software is an **exploratory computational research suite**. Sequence-based PWM scoring, inferred ancestral consensus states, proximity overlaps, and deep-learning chromatin/transcription predictions (`AlphaGenome`) serve as high-resolution *in silico* hypotheses and molecular proxies. They are designed to prioritize candidate loci and variants for downstream wet-lab functional validation (such as Massively Parallel Reporter Assays [MPRA], ChIP-seq, or CRISPRi/a editing in 3D mammary organoids).

Since the total execution of the repository takes approximately 20 hours, here is a folder with all the data and results after execution: https://drive.google.com/drive/u/0/folders/1mF9IdvVCiWCk4xWaQ_8Kg0_93LUjDDF-
<img width="4500" height="3300" alt="resumen_validacion_alphagenome_mamario" src="https://github.com/user-attachments/assets/42de1eed-c065-4d66-907a-5af6d07c0670" />


---

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.
