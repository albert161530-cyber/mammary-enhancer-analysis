# Comparative Genomics Analysis Suite: Primate Regulatory Evolution

[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![C++](https://img.shields.io/badge/C++-pybind11-orange.svg)](https://github.com/pybind/pybind11)
[![AlphaGenome](https://img.shields.io/badge/Deep%20Learning-AlphaGenome-green.svg)](https://github.com/google-deepmind/alphagenome)
[![Comparative Genomics](https://img.shields.io/badge/Comparative-Genomics-purple.svg)](#comparative-genomics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A high-performance, multi-stage computational genomics pipeline investigating the regulatory evolution of permanently enlarged human breasts.**

This project evaluates two competing evolutionary hypotheses: whether the emergence of permanently enlarged human breasts was driven more by the loss of activity in androgen-associated regulatory brakes, particularly androgen receptor (*AR*) and estrogen receptor beta (*ESR2*) pathways, or by the gain of activity in estrogen-associated regulatory accelerators, particularly estrogen receptor alpha (*ESR1*) pathways.

The framework combines comparative genomics across primates, a 30-mammal Multiple Alignment Format (MAF) alignment, cross-species epigenomic evidence, native C++ sequence-analysis kernels, Python-based statistical and physiological modeling, and AlphaGenome deep-learning predictions.

The primary objective is to quantify and compare the predicted regulatory contributions of these two evolutionary mechanisms, identify candidate causal variants, test whether independent computational approaches converge on the same mechanism, and prioritize regulatory regions for experimental investigation.

## Table of Contents

- [Overview](#overview)
- [Research Questions](#research-questions)
- [Comparative Genomics](#comparative-genomics)
- [Regulatory Biology and the Brake-versus-Accelerator Model](#regulatory-biology-and-the-brake-versus-accelerator-model)
- [Architecture and Computational Methods](#architecture-and-computational-methods)
- [Epigenomic and Tissue-Specific Analysis](#epigenomic-and-tissue-specific-analysis)
- [AlphaGenome Prediction Framework](#alphagenome-prediction-framework)
- [Endocrine and Physiological Modeling](#endocrine-and-physiological-modeling)
- [External Tools and Data Resources](#external-tools-and-data-resources)
- [Repository Structure](#repository-structure)
- [Installation and Environment Setup](#installation-and-environment-setup)
- [AlphaGenome API Key Configuration](#alphagenome-api-key-configuration)
- [Input Data and Configuration](#input-data-and-configuration)
- [Execution Workflow](#execution-workflow)
- [Statistical Analysis and Outputs](#statistical-analysis-and-outputs)
- [Translational Research Directions](#translational-research-directions)
- [Computational Requirements and Runtime](#computational-requirements-and-runtime)
- [Reproducibility and Interpretation](#reproducibility-and-interpretation)
- [Data and Results](#data-and-results)
- [References and External Documentation](#references-and-external-documentation)
- [License](#license)

## Overview

Human breasts exhibit permanent enlargement beyond the developmental periods associated with pregnancy and lactation. This project investigates whether evolutionary changes in non-coding regulatory DNA contributed to that phenotype by altering the balance between androgen-associated and estrogen-associated transcriptional regulation in mammary tissue.

Instead of limiting the analysis to protein-coding mutations or previously annotated enhancers, the pipeline examines sequence changes across candidate regulatory regions and evaluates their genomic, epigenomic, transcriptional, and modeled physiological consequences.

The analysis focuses on four genes:

| Gene | Primary role in the analytical framework |
|---|---|
| `AR` | Androgen receptor; candidate androgen-associated regulatory brake |
| `ESR1` | Estrogen receptor alpha; candidate estrogen-associated regulatory accelerator |
| `ESR2` | Estrogen receptor beta; candidate regulatory component of the estrogen signaling network |
| `CYP19A1` | Aromatase; candidate regulator of local estrogen biosynthesis |

These categories represent the working mechanistic hypothesis. The pipeline tests their regulatory contributions rather than assuming that every enhancer associated with a particular receptor necessarily has the same functional effect.

### Main analytical principle

The project separates three related questions:

1. **What changed during evolution?** Identify human-lineage substitutions, insertions, deletions, and other candidate regulatory differences by comparing primate sequences against a multispecies mammalian alignment.
2. **What is the predicted regulatory effect?** Evaluate sequence motifs, tissue-specific regulatory activity, and AlphaGenome-predicted changes in molecular outputs, individually and jointly.
3. **Which mechanism contributes more to the modeled phenotype?** Integrate the regulatory results into comparative statistical summaries and an endocrine model that isolates enhancer-derived contributions to androgen-associated versus estrogen-associated activity.

The central comparison is between `AR_loss` and `ESR1_gain`, while incorporating `ESR2`, `CYP19A1`, and additional interacting regulatory regions.

## Research Questions

The pipeline is designed to investigate the following questions.

### Evolutionary origin

- Which candidate regulatory changes distinguish the human lineage from chimpanzees, bonobos, gorillas, and orangutans?
- Which substitutions, deletions, and insertions occur within or near conserved regulatory regions?
- Which human-specific substitutions coincide with predicted changes in transcription-factor motif strength?
- Are candidate regulatory sequences associated with activity in human mammary tissue but different activity in the corresponding chimpanzee tissue?
- Do variants affecting repressive regulatory mechanisms contribute more to the modeled change than variants affecting activating mechanisms?

### Functional impact

- What is the predicted effect of each candidate variant on relevant molecular outputs?
- Do variants behave differently when analyzed individually versus jointly?
- Which loci exhibit non-additive interactions between variants?
- Do sequence-based motif scores and deep-learning predictions produce concordant rankings?
- Which candidates show the strongest predicted effects in mammary epithelial, stromal, or adipose compartments?

### Physiological interpretation

- How do the inferred regulatory differences alter the modeled balance between androgen-associated and estrogen-associated activity?
- How does the model behave across developmental and reproductive states?
- Which regulatory candidates merit further investigation in breast cancer, gynecomastia, or mammary reconstruction research?

## Comparative Genomics

### Multispecies alignment

The initial evolutionary scan uses a MAF alignment containing 30 mammalian species. The analysis compares the human sequence with four closely related non-human primates:

- *Pan troglodytes* — chimpanzee
- *Pan paniscus* — bonobo
- *Gorilla gorilla* — gorilla
- *Pongo abelii* — Sumatran orangutan

The additional mammalian sequences provide broader evolutionary context for distinguishing lineage-specific changes from ancestral or more widely conserved sequence states.

The analysis concentrates on regulatory regions associated with `AR`, `ESR1`, `ESR2`, and `CYP19A1`, including annotated regulatory elements and candidate regions that merit scanning even when they lack existing enhancer annotations.

### Variant classes

The pipeline considers several classes of sequence variation:

| Variant class | Analytical purpose |
|---|---|
| Single-nucleotide variants (SNVs) | Identify substitutions that may alter transcription-factor recognition or regulatory sequence activity |
| Short insertions and deletions (indels) | Evaluate sequence-length changes and their effects on motif organization |
| Larger deletions | Identify candidate regulatory losses and missing ancestral sequence segments |
| Human-specific insertions | Evaluate newly introduced regulatory sequence and possible motif gains |
| Combined variants | Test the joint consequences of multiple changes within the same modeled haplotype or genomic interval |

The distinction between a sequence change and a regulatory change is fundamental. A single nucleotide substitution can alter transcription-factor recognition without deleting the enhancer itself. Conversely, a sequence can remain physically present while its activity changes because its motif architecture or chromatin context has changed.

For this reason, the analysis does not classify every human-specific variant as an enhancer gain or loss solely on the basis of its sequence difference.

### Ancestral and lineage-specific comparisons

Where alignment coverage and the implemented reconstruction procedure permit, the pipeline uses multispecies conservation to contextualize candidate ancestral and derived states.

A human-specific state is interpreted relative to the aligned non-human sequences, their coverage, and the ancestral-state assumptions used in the analysis. Missing alignment data and ambiguous states must not be automatically treated as evidence of a genuine deletion.

## Regulatory Biology and the Brake-versus-Accelerator Model

The principal mechanistic comparison contrasts candidate changes that reduce androgen-associated regulatory effects with those that increase estrogen-associated regulatory effects.

### Model A: Loss of regulatory brakes

`AR_loss` represents the candidate evolutionary reduction of androgen-associated regulatory repression.

The analysis evaluates whether sequence changes in relevant regulatory regions reduce predicted regulatory activity, alter transcription-factor motif recognition, or change downstream molecular predictions associated with the androgen receptor pathway.

`ESR2` is included because estrogen receptor beta may contribute to the broader receptor-regulatory balance. Its contribution is evaluated separately rather than automatically equated with the effect of `AR`.

### Model B: Gain of regulatory accelerators

`ESR1_gain` represents candidate evolutionary increases in estrogen-associated regulatory activity.

The pipeline examines variants that may increase predicted regulatory output or alter the activity of enhancers associated with estrogen receptor alpha.

`CYP19A1` is included because aromatase regulates estrogen biosynthesis and can provide an additional mechanistic connection between local endocrine regulation and mammary tissue behavior.

### Comparative dominance

The central result is not simply the number of candidate variants assigned to each category. The framework compares predicted effect sizes, regulatory directions, tissue specificity, and the contribution of interacting variants.

A larger number of variants in one category does not automatically imply a larger cumulative effect. The consolidated analyses are designed to evaluate the magnitude and consistency of the predicted contributions.

The primary outcome is the relative regulatory contribution of the two candidate mechanisms under the model's stated assumptions.

## Architecture and Computational Methods

The software combines Python orchestration with a native C++ acceleration layer.

### Python orchestration

Python manages:

- MAF parsing and variant annotation.
- Regulatory element and epigenomic track processing.
- Input and output management.
- AlphaGenome API requests.
- Variant batching and checkpointing.
- Statistical tests and aggregation.
- Physiological model calculations.
- CSV exports and publication-quality figures.

### Native C++ kernels

The `genoma_kernels` extension is compiled through `pybind11` and provides native implementations for computationally intensive sequence operations.

These include:

- Position Weight Matrix (PWM) scoring.
- Centered motif searches.
- Reverse-complement sequence analysis.
- Repeated sequence-scanning operations.
- Motif-level calculations across candidate windows and variant configurations.

The purpose of the native extension is to reduce the computational overhead of repeated sequence operations while allowing Python to retain control of the overall analysis.

### Matrix-based evaluation

The sequence-analysis and modeling stages represent candidate regions, motif scores, regulatory effects, and comparative scenarios as structured arrays or matrices where appropriate.

This enables batch processing and aggregation across multiple variants and candidate loci.

The pipeline can then pass the resulting candidate sets to AlphaGenome for individual and combined predictions, instead of treating every downstream comparison as an isolated manual operation.

## Epigenomic and Tissue-Specific Analysis

A central component of the framework is the integration of genomic sequence differences with regulatory activity information.

A sequence motif change alone does not establish that an enhancer is active in mammary tissue. Regulatory annotations and epigenomic signals provide additional context for evaluating whether a candidate region is relevant to the tissue under investigation.

### Human and chimpanzee epigenomes

Where suitable datasets are available, the analysis compares the regulatory context of corresponding human and chimpanzee regions.

The key question is whether a sequence difference coincides with a difference in regulatory activity between the two species.

For example, a candidate region may retain a recognizable enhancer-like sequence in both species while containing a human-specific substitution predicted to alter a relevant transcription-factor motif. Comparing the corresponding epigenomic evidence helps determine whether the proposed regulatory difference is consistent with the tissue-specific activity data.

The interpretation depends on the availability of compatible biological samples, tissue annotations, genome assemblies, and assay types.

### Mammary and adipose compartments

The analysis prioritizes the following tissue or cell categories where data are available:

- Mammary luminal epithelial cells.
- Mammary myoepithelial cells.
- Breast fibroblasts and other stromal compartments.
- Mammary-associated adipose tissue.
- Additional relevant epithelial, endocrine, and adipose datasets when needed for context.

This tissue stratification is important because regulatory activity can differ substantially between mammary epithelium, stromal tissue, and adipose cells.

### Regulatory annotations

The current framework includes, or is designed to work with, the following annotation categories:

- SCREEN Registry V4 candidate cis-regulatory elements (cCREs), including enhancer-like signatures.
- EpiMap ChromHMM regulatory-state annotations.
- ATAC-seq and DNase-seq accessibility tracks.
- Histone modification ChIP-seq tracks.
- Transcription-factor binding tracks.
- RNA-seq and related transcriptional measurements.
- BigWig-formatted signal tracks accessed through `pyBigWig`.

Annotated regulatory elements are combined with sequence-level scanning so that the investigation can also consider candidate regions outside previously catalogued enhancers.

## AlphaGenome Prediction Framework

AlphaGenome is used as a deep-learning component for evaluating the predicted molecular consequences of selected sequence changes.

Official project: https://github.com/google-deepmind/alphagenome

The pipeline sends prioritized candidate variants or sequence configurations to AlphaGenome through its API, then consolidates the returned outputs for downstream comparisons.

### Individual-variant predictions

Single-variant analysis evaluates each candidate separately.

This stage is used to identify variants with substantial predicted effects, compare the direction of regulatory changes, and rank candidate variants across the target loci and relevant tracks.

### Joint-haplotype predictions

Combined analysis evaluates the joint effect of multiple variants within the same modeled sequence context.

This is important because the joint effect of a group of variants is not necessarily equal to the sum of their isolated effects.

For each appropriate comparison, the analysis examines whether the combined prediction agrees with the sum of the individual predictions or exhibits a non-additive effect.

### Human versus reference or reconstructed sequence configurations

The scripts compare model-supported alternative sequence configurations, including relevant reference, derived, or reconstructed sequence states where implemented.

The corresponding sequence coordinates and allele definitions must match the selected reference assembly and the inputs expected by the AlphaGenome client.

### Molecular outputs

Depending on the requested prediction types and available model tracks, the analysis considers outputs such as:

- ATAC-seq and DNase-seq.
- Histone ChIP-seq.
- Transcription-factor ChIP-seq.
- CAGE and PRO-cap.
- RNA-seq.
- Additional supported regulatory outputs where configured.

The project configuration describes a collection of 279 cell and tissue tracks for its prediction and aggregation workflow. The available results depend on the requested prediction types, returned model outputs, and the track-selection logic implemented in the scripts.

### Three complementary evaluation modes

The framework compares three analytical layers:

**1. Sequence and motif analysis.** C++-accelerated PWM scoring and related motif calculations estimate how candidate substitutions or indels alter sequence-based regulatory evidence.

**2. Individual-variant prediction.** AlphaGenome estimates the molecular effect of each selected variant independently.

**3. Combined-variant prediction.** AlphaGenome evaluates multiple variants jointly, allowing direct comparison between individual and combined predictions.

Agreement across these layers supports a more consistent candidate ranking. Disagreement is retained as an analytical result because it can identify non-additivity, context dependence, or differences between motif-based and deep-learning estimates.

The downstream consolidation compares `AR_loss` and `ESR1_gain` using the outputs generated by the configured analyses.

### AlphaGenome processing and recovery

The scripts use batch processing to organize requests across candidate genomic regions. Where implemented, checkpointing and automatic resume allow an interrupted insertion or variant-prediction run to continue without unnecessarily repeating completed work.

Because this stage depends on an external API, total execution time also depends on request volume, permitted query rates, successful responses, and the configuration of the individual scripts.

## Endocrine and Physiological Modeling

The project includes a modular physiological model that translates regulatory activity estimates into a comparative androgen-associated versus estrogen-associated signaling framework.

The purpose is to examine how the predicted enhancer contributions affect the modeled balance between the two regulatory sides.

### Controlled baseline

The primary comparative model holds basal hormone levels and receptor affinities equal when isolating the contribution of enhancer-derived regulation.

This separates the intended regulatory comparison from differences caused by starting hormone concentrations or receptor-binding affinity.

The resulting values describe the model's controlled comparison. They must not be interpreted as direct measurements of circulating hormone concentrations, receptor occupancy, or total physiological signaling in an individual.

### Physiological scenarios

The simulator includes, or is configured to compare, relevant developmental and physiological states:

- Fetal developmental stages.
- Male and female developmental profiles.
- Tanner stages.
- Menstrual-cycle stages.
- Normal pregnancy.
- Lactation.
- Gynecomastia.
- Polycystic ovary syndrome (PCOS).
- Wild-type human regulatory configurations.
- Wild-type chimpanzee regulatory configurations.
- Humanized-chimpanzee configurations incorporating selected human regulatory changes.

These scenarios provide a common modeling framework for examining how the candidate regulatory differences behave under distinct physiological inputs or state-specific configurations.

Where the model uses the same basal concentrations across scenarios, the outputs isolate the regulatory component under that controlled assumption rather than reproducing the actual endocrine concentrations of each physiological state.

### Comparative configurations

The simulator compares three main genetic-regulatory configurations:

1. Wild-type human.
2. Wild-type chimpanzee.
3. Humanized chimpanzee, incorporating the selected human regulatory changes.

This comparison is designed to evaluate how much of the modeled difference follows from the introduction of the selected regulatory states into an otherwise comparative configuration.

The input assumptions, normalization rules, and exact implementation of each scenario should be retained with the corresponding output files.

## External Tools and Data Resources

The main pipeline integrates comparative sequence analysis and AlphaGenome prediction with a set of external tools and resources for regulatory annotation, expression context, and network-level investigation.

The specific capabilities available in a given run depend on the scripts, data files, endpoints, and packages configured for that run.

| Tool or resource | Intended analytical role |
|---|---|
| ENCODE REST API | Retrieve or cross-reference genomic experiments, biosamples, and regulatory datasets |
| `pyBigWig` | Read and query signal tracks stored in BigWig format |
| GSEApy | Gene-set enrichment and pathway-level interpretation of appropriate gene lists or ranked results |
| CellOracle | Regulatory-network modeling and in silico perturbation analysis when compatible network inputs are available |
| SCENIC | Gene-regulatory network and transcription-factor activity analysis using compatible expression data |
| ARCHS4 | Expression and co-expression context using its available transcriptomic resources |
| GEPIA2 / EpiA2 | Expression-analysis cross-checks using the configured expression resource and supported interface |
| SCREEN | Candidate cis-regulatory element annotations |
| EpiMap | Epigenomic regulatory-state annotations |
| JASPAR | Position Weight Matrices for transcription-factor motif analysis |
| AlphaGenome | Sequence-based deep-learning predictions of regulatory and transcriptional effects |

These resources serve complementary purposes. They do not all measure the same biological quantity, and their outputs should not be combined as if they were interchangeable scores.

In particular, expression-based network inference, epigenomic signal, motif strength, and sequence-based molecular predictions represent different types of evidence.

## Repository Structure

### Core pipeline scripts

| Script | Stage | Main responsibility |
|---|---|---|
| `Genoma1.py` | 1 | Comparative MAF parsing, motif scanning, regulatory annotation, and candidate-variant generation |
| `automatizar_bloques.py` | 2 | Individual-variant and joint-haplotype AlphaGenome predictions for substitutions and short indels |
| `auditoria_inserciones.py` | 3 | Human-specific insertion analysis and AlphaGenome evaluation |
| `analisis_posthoc_mamario.py` | 4 | Tissue-specific statistical analyses, variant ranking, epistasis evaluation, and figure generation |
| `generar_reporte_multiclase.py` | 5 | Consolidation of sequence-variant and insertion results into a common evolutionary summary |
| `modelo_modular_potenciadores.py` | 6 | Comparative endocrine and physiological simulation |

### Main input files

The exact paths are defined by the current script configuration. The project uses or expects the following input categories:

- A 30-mammal MAF alignment.
- Human and corresponding primate genomic sequences or alignable regions.
- SCREEN, EpiMap, or other configured regulatory annotations.
- BigWig epigenomic tracks when available.
- Candidate-variant tables generated during the initial stage.
- AlphaGenome API credentials.
- State-specific parameters and input tables for the physiological simulator.

### Main output files

The core workflow produces intermediate and consolidated tables, including:

- `variantes_todos_loci.csv`
- `variantes_{gene}_anotadas.csv`
- `deleciones_todos_loci.csv`
- `comparacion_modelos_3_escenarios.csv`
- `resultados_modelo_hormonal.csv`
- `resultados_impacto_conjunto_comparacion.csv`
- `resultados_impacto_conjunto_por_track.csv`
- `resultados_impacto_individual_por_track.csv`
- `resultados_impacto_conjunto_regiones.csv`
- `variantes_inserciones_loci.csv`
- `tabla_5_impacto_inserciones_mamario.csv`
- `tabla_5_resumen_inserciones_por_modalidad.csv`

Later stages aggregate these intermediate results into statistical summaries, driver rankings, and an overall multi-class evolutionary comparison.

## Installation and Environment Setup

The project uses a Python environment with a native C++ toolchain. The tested workflow also required installing `pyBigWig` in a Conda environment running inside a Docker container with Ubuntu.

This separates the genomic processing environment from the Windows host and provides an appropriate Linux environment for compiling and running the relevant dependencies.

### 1. Prepare Ubuntu inside Docker

Open a shell inside the Ubuntu container and install the native build tools:

```bash
apt-get update
apt-get install -y build-essential git ca-certificates
```

If the container already includes these packages, this step does not need to be repeated.

### 2. Create the Conda environment

Using Mamba:

```bash
mamba create -n primate-regulatory python=3.11 pip -y
mamba activate primate-regulatory
```

If only Conda is installed, substitute `conda` for `mamba`.

### 3. Install `pyBigWig` through Conda

Install the package in the active Linux environment:

```bash
mamba install -c conda-forge -c bioconda pybigwig -y
```

This step is particularly relevant when the pipeline needs to access epigenomic signal stored in BigWig files.

The installation must take place inside the same Conda environment used to execute the scripts.

### 4. Install the Python packages

```bash
python -m pip install --upgrade pip

python -m pip install \
    numpy \
    pandas \
    scipy \
    matplotlib \
    requests \
    certifi \
    biopython \
    alphagenome \
    pybind11 \
    setuptools
```

Additional packages are required when running optional analysis modules that use CellOracle, GSEApy, SCENIC, or other external interfaces. Install and configure those packages according to the corresponding module's actual requirements.

### 5. Verify `pyBigWig`

Run the following command inside the activated environment:

```bash
python -c "import pyBigWig; print('pyBigWig import successful')"
```

If the import succeeds, Python can load the installed module from that environment.

### 6. Verify the native C++ toolchain

```bash
g++ --version
python -c "import pybind11; print(pybind11.get_include())"
```

The scripts that build the native `genoma_kernels` extension must run in an environment where the compiler, Python development environment, and `pybind11` are available.

### 7. Verify the AlphaGenome client

```bash
python -c "import alphagenome; print('AlphaGenome client import successful')"
```

Successful import verifies that the client package is available. A successful live prediction also requires a valid API key and access to the API.

## AlphaGenome API Key Configuration

The AlphaGenome API key is stored in a local file named `api_key.txt` in the project root.

The expected structure is:

```text
project_root/
├── api_key.txt
├── Genoma1.py
├── automatizar_bloques.py
├── auditoria_inserciones.py
├── analisis_posthoc_mamario.py
├── generar_reporte_multiclase.py
└── modelo_modular_potenciadores.py
```

The file should contain the API key in the format expected by the current scripts.

Example of the file structure:

```text
YOUR_ALPHAGENOME_API_KEY
```

Replace the placeholder with your own key. Do not include quotation marks unless the code specifically expects them.

The prediction scripts that use this configuration include:

- `automatizar_bloques.py`
- `auditoria_inserciones.py`

Keep `api_key.txt` out of version control. Add the following entry to `.gitignore`:

```gitignore
api_key.txt
```

The API key should not appear in source code, committed configuration files, screenshots, logs, or publicly shared outputs.

Official AlphaGenome project and installation documentation:

https://github.com/google-deepmind/alphagenome

## Input Data and Configuration

Before running the complete pipeline, verify that the required input files and configuration parameters are available.

### Genome and alignment consistency

The MAF input must contain the relevant human and non-human primate alignments. The genome assembly and coordinate system used for variant annotation must match the reference expected by the downstream scripts.

Coordinate mismatches can produce incorrect overlaps between variants, enhancer annotations, and epigenomic tracks.

### Epigenomic tracks

BigWig files and other epigenomic inputs must correspond to the relevant genome assembly. When comparing human and chimpanzee data, confirm the species, tissue, biosample, and assay represented by each track.

The availability of matched samples determines which cross-species comparisons can be performed directly.

### Candidate tables

The output of the first stage supplies candidate variants to the AlphaGenome stages. Intermediate tables must have the columns and variant representations expected by the downstream scripts.

The insertion auditor operates on its own insertion-candidate workflow and produces outputs used by the multi-class consolidation stage.

### API access

Ensure that `api_key.txt` is located where the scripts expect it and that the active Python environment can access the AlphaGenome client.

For large runs, retain the script checkpoints and existing results where the resume functionality depends on them.

## Execution Workflow

The full analysis is organized into six sequential stages.

### Step 1 — Parse MAF alignments, annotate regulatory regions, and scan motifs

```bash
python Genoma1.py
```

This stage processes the multispecies alignment, identifies candidate sequence differences, performs motif-level analysis, and combines the findings with configured regulatory annotations.

Main outputs:

- `variantes_todos_loci.csv`
- `variantes_{gene}_anotadas.csv`
- `deleciones_todos_loci.csv`
- `comparacion_modelos_3_escenarios.csv`
- `resultados_modelo_hormonal.csv`

The specific outputs depend on the available inputs and the current script configuration.

### Step 2 — Evaluate substitutions and short indels using AlphaGenome

```bash
python automatizar_bloques.py
```

This stage processes the selected variants and generates individual-variant and combined-variant predictions.

The analysis compares the predicted molecular consequences associated with `AR_loss` and `ESR1_gain` and stores outputs for downstream statistical evaluation.

Main outputs:

- `resultados_impacto_conjunto_comparacion.csv`
- `resultados_impacto_conjunto_por_track.csv`
- `resultados_impacto_individual_por_track.csv`
- `resultados_impacto_conjunto_regiones.csv`

### Step 3 — Audit human-specific insertions

```bash
python auditoria_inserciones.py
```

This stage evaluates candidate human-specific insertions and their predicted regulatory consequences within the target loci and available mammary or adipose tracks.

Where supported by the current script configuration, the process uses checkpointing and automatic resume to preserve completed work.

Main outputs:

- `variantes_inserciones_loci.csv`
- `tabla_5_impacto_inserciones_mamario.csv`
- `tabla_5_resumen_inserciones_por_modalidad.csv`

### Step 4 — Perform post-hoc statistical analysis

```bash
python analisis_posthoc_mamario.py
```

This stage aggregates prediction results, applies tissue-specific filters, performs the configured statistical tests, ranks candidate variants, and generates figures.

The analysis includes exact binomial tests, paired Wilcoxon tests, locus-level non-additivity analysis, and candidate ranking across mammary and adipose compartments.

The current project configuration identifies a mammary/adipose subset of `n = 178`. The sample count and the inclusion criteria should be verified against the actual generated tables when interpreting a particular run.

Main outputs:

- `tabla_1_resumen_global_y_mamario_por_modalidad.csv`
- `tabla_2_epistasis_no_aditividad_por_locus.csv`
- `tabla_3_top_variantes_drivers_mamario_adiposo.csv`
- `tabla_4_cross_talk_factores_directos_MCF7.csv`
- `resumen_validacion_alphagenome_mamario.png`

### Step 5 — Consolidate the evolutionary comparison

```bash
python generar_reporte_multiclase.py
```

This stage integrates the relevant variant-prediction results and the human-specific insertion analysis into a common multi-class summary.

Main output:

- `tabla_6_consolidado_multiclase_evolutivo.csv`

The consolidated report is intended to make the relative contributions of the candidate mechanisms easier to compare across the evaluated variant classes.

### Step 6 — Run the physiological model

```bash
python modelo_modular_potenciadores.py
```

This stage runs the modular enhancer and endocrine simulation using the configured input parameters.

It compares the modeled regulatory behavior of wild-type human, wild-type chimpanzee, and humanized-chimpanzee configurations.

Main outputs:

- `entradas_modelo_hormonal.csv`
- Comparative console or CSV results generated by the current implementation.

Review the output fields and model assumptions before interpreting the final numerical values.

### Full pipeline command

If the alignment files, annotations, model inputs, dependencies, and API key are ready, execute the stages in sequence:

```bash
python Genoma1.py && \
python automatizar_bloques.py && \
python auditoria_inserciones.py && \
python analisis_posthoc_mamario.py && \
python generar_reporte_multiclase.py && \
python modelo_modular_potenciadores.py
```

The `&&` operators stop the sequence when a command fails. This prevents later stages from silently running against missing or incomplete outputs.

If only one stage fails, correct the error and resume from the appropriate stage, checking whether the relevant intermediate files and checkpoint state are valid.

## Statistical Analysis and Outputs

The post-hoc analysis evaluates the consistency and relative magnitude of the predicted regulatory effects across the selected tracks, variants, and candidate loci.

### Global and tissue-specific summaries

`tabla_1_resumen_global_y_mamario_por_modalidad.csv` summarizes the comparative results by analysis category and tissue compartment.

### Non-additive interactions

`tabla_2_epistasis_no_aditividad_por_locus.csv` evaluates the relationship between combined predictions and the corresponding individual-variant effects.

This analysis distinguishes additive behavior from non-additive outcomes under the joint-prediction procedure used by the scripts.

### Candidate drivers

`tabla_3_top_variantes_drivers_mamario_adiposo.csv` ranks candidate variants by the criteria implemented in the post-hoc analysis.

These rankings can be used to prioritize variants for more focused molecular testing.

### Direct-factor cross-talk

`tabla_4_cross_talk_factores_directos_MCF7.csv` contains the configured cross-talk analysis associated with direct-factor data in the MCF-7 context.

Its interpretation should follow the variables and evidence represented in the actual table, including the distinction between measured data and modeled predictions.

### Insertion-specific effects

`tabla_5_impacto_inserciones_mamario.csv` and `tabla_5_resumen_inserciones_por_modalidad.csv` summarize the insertion analysis.

### Integrated comparison

`tabla_6_consolidado_multiclase_evolutivo.csv` combines the relevant outputs for the final multi-class evolutionary comparison.

Together, these files provide the main outputs for comparing candidate regulatory brakes and accelerators.

The direction, magnitude, statistical support, and cross-method agreement of the generated results determine the interpretation of the analysis. The pipeline does not establish a winner in advance.

## Translational Research Directions

The regulatory framework is also intended to prioritize mechanistic hypotheses for future experimental investigation.

### Breast cancer

The `AR`–`ESR1` regulatory axis motivates investigation of whether particular regulatory elements influence hormone-associated signaling in mammary cells and breast cancer models.

The candidate ranking can help identify regulatory regions for focused functional assays and subsequent investigation of their relevance to tumor biology.

### Gynecomastia

The model provides a framework for investigating whether particular regulatory configurations alter the modeled balance between androgen-associated and estrogen-associated activity in male mammary tissue.

The output can support the selection of regulatory candidates for future experimental study.

### Mammary construction and reconstruction

The framework explores whether regulatory principles associated with mammary development, epithelial regulation, and local adipose biology can inform future investigations of mammary tissue engineering and reconstructive biology.

These directions are hypothesis-generating research applications. The computational outputs do not constitute validated cancer treatments, gynecomastia therapies, or established protocols for breast construction or reconstruction.

Any therapeutic or tissue-engineering strategy derived from these findings requires independent functional testing and appropriate subsequent validation.

## Computational Requirements and Runtime

The pipeline combines comparatively inexpensive local operations with computationally intensive deep-learning API requests.

Actual execution time depends on:

- The number of candidate variants and genomic intervals.
- The size and structure of the alignment input.
- The number and configuration of prediction tracks.
- The number of individual and joint AlphaGenome requests.
- The availability and use of checkpoints.
- CPU performance and local data processing.
- API availability and permitted request rates.

The complete workflow has an estimated runtime of approximately 20 hours in the current configuration. This is a project-specific runtime estimate, not a guaranteed execution time on every system.

Running the native kernels inside an appropriately configured Linux/Conda environment and preserving valid intermediate outputs helps make long-running analyses easier to manage.

## Reproducibility and Interpretation

For a reproducible analysis, retain:

- The genome assembly and MAF input version.
- The species represented in the alignment.
- The annotations and epigenomic track versions.
- The configuration parameters used by each script.
- The AlphaGenome client version and prediction settings.
- The input variant tables and resulting intermediate outputs.
- The statistical procedures and physiological model assumptions.

Several distinctions are essential when interpreting the results.

**Sequence differences versus regulatory effects:** A human-specific substitution is evidence of sequence divergence. Its predicted effect on enhancer activity is a separate analytical result.

**Predicted versus measured activity:** PWM scores and AlphaGenome outputs are computational estimates, not direct experimental measurements of enhancer activity.

**Individual versus combined effects:** Single-variant effects do not necessarily predict the effect of a joint sequence configuration. The joint-prediction and post-hoc stages address this distinction.

**Regulatory balance versus complete physiology:** The endocrine simulator isolates contributions under its declared assumptions. It does not independently reconstruct all endocrine, developmental, cellular, and systemic processes responsible for mammary phenotypes.

**Model comparison versus experimental confirmation:** Cross-method agreement strengthens the internal consistency of a computational result. Functional experiments are required to establish the biological effect of specific candidate regulatory sequences.

Potential follow-up experiments include massively parallel reporter assays (MPRA), targeted ChIP-seq or accessibility assays, CRISPRi/CRISPRa perturbation, and experiments in appropriate mammary cellular models or three-dimensional organoids.

## Data and Results

The folder below contains the project's collected data and results after execution:

[Open the project data and results on Google Drive](https://drive.google.com/drive/u/0/folders/1mF9IdvVCiWCk4xWaQ_8Kg0_93LUjDDF-)

### AlphaGenome mammary validation summary

![AlphaGenome mammary validation summary](https://github.com/user-attachments/assets/42de1eed-c065-4d66-907a-5af6d07c0670)

The figure corresponds to the project's mammary AlphaGenome validation summary. The numerical interpretation should be based on the associated CSV outputs, statistical procedures, and prediction settings.

## References and External Documentation

### Core technologies

- Python: https://www.python.org/
- pybind11: https://github.com/pybind/pybind11
- AlphaGenome: https://github.com/google-deepmind/alphagenome
- pyBigWig: https://github.com/deeptools/pyBigWig
- JASPAR: https://jaspar.genereg.net/

### Regulatory and epigenomic resources

- ENCODE Project: https://www.encodeproject.org/
- ENCODE SCREEN Registry: https://screen.encodeproject.org/
- EpiMap: https://epigenomegateway.wustl.edu/
- ARCHS4: https://maayanlab.cloud/archs4/
- GEPIA2: http://gepia2.cancer-pku.cn/

### Regulatory-network and enrichment tools

- GSEApy: https://gseapy.readthedocs.io/
- CellOracle: https://github.com/morris-lab/CellOracle
- SCENIC: https://github.com/aertslab/SCENIC

Each resource should be cited according to the analysis for which it was used. Cite the underlying datasets and software versions in any resulting manuscript or report.

## License

Distributed under the MIT License. See `LICENSE` for more information.
