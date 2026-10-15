# main_traslacional.py
"""
Programa Principal de Validación Traslacional In Silico.
Orquesta las pruebas de:
1. Supresión Clínica (Ginecomastia y Cáncer de Mama ER+ en MCF-7)
2. Reconstrucción Estética e Ingeniería Tisular Adiposa
3. Auditoría de simetría y generación de figuras en inglés para pitch académico.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from config_dianas import TOP_AR_LOSS_DRIVERS, TISSUE_MODALITY_RATIOS
from modulo_supresion_clinica import ClinicalSuppressionSimulator
from modulo_reconstruccion_estetica import AestheticReconstructionSimulator

def ejecutar_pruebas_unitarias():
    """Verifica que los módulos y datos cumplan con la integridad matemática antes de simular."""
    print("[1/4] Ejecutando pruebas de verificación e integridad de módulos...")
    
    # Prueba 1: Verificar que las 10 variantes driver estén completas y ordenadas
    assert len(TOP_AR_LOSS_DRIVERS) == 10, "Error: Deben existir exactamente 10 variantes driver."
    scores = [d.mean_l1_impact for d in TOP_AR_LOSS_DRIVERS]
    assert scores == sorted(scores, reverse=True), "Error: Las variantes no están ordenadas por impacto L1."
    
    # Prueba 2: Verificar que la respuesta monotónica del simulador clínico sea coherente
    sim_clinico = ClinicalSuppressionSimulator()
    df_clin = sim_clinico.run_all_clinical_panels()
    rescate_1 = df_clin.loc[df_clin["Num_Guias_sgRNA"] == 1, "Rescate_Epistatico_Haplotype_Pct"].values[0]
    rescate_4 = df_clin.loc[df_clin["Num_Guias_sgRNA"] == 4, "Rescate_Epistatico_Haplotype_Pct"].values[0]
    assert rescate_4 > rescate_1, "Error lógico: El panel de 4 guías debe rescatar más señal que 1 sola guía."

    # Prueba 3: Verificar selectividad adiposa en DNase y ChIP-TF
    sim_estetico = AestheticReconstructionSimulator()
    df_sel = sim_estetico.evaluate_adipose_selectivity()
    dnase_ratio = df_sel.loc[df_sel["Modalidad_Epigenetica"] == "DNASE", "Indice_Selectividad_Adiposa"].values[0]
    assert dnase_ratio > 1.5, "Error: El índice de selectividad adiposa en DNase no refleja los tracks n=25."

    print("      ✔ Todas las pruebas de integridad pasaron exitosamente.\n")

def generar_grafica_traslacional_en(df_clinico: pd.DataFrame, df_selectividad: pd.DataFrame, output_path: str):
    """Genera una figura de publicación en inglés con los resultados de las pruebas traslacionales."""
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.5))
    fig.subplots_adjust(wspace=0.28)

    # Panel A: Rescate In Silico por Panel Terapéutico (Ginecomastia & Cáncer ER+)
    panels = [
        "Single Top-1\n(ESR1)",
        "Gynecomastia Dual\n(ESR1 + CYP19A1)",
        "Oncology Triad\n(Top-3 ESR1)",
        "Severe Gynecomastia\n(4-Guide Panel)",
        "Full Top-10\nAncestral Rescue"
    ]
    linear_vals = df_clinico["Rescate_Lineal_Pct"].values
    epistatic_vals = df_clinico["Rescate_Epistatico_Haplotype_Pct"].values
    rnaseq_vals = df_clinico["Supresion_Transcripcional_RNAseq_Pct"].values

    x = np.arange(len(panels))
    width = 0.26

    axs[0].bar(x - width, linear_vals, width, color="#7f7f7f", edgecolor="black", label="Additive L1 Rescue (%)")
    axs[0].bar(x, epistatic_vals, width, color="#1f77b4", edgecolor="black", label="Epistatic Chromatin Rescue (%)")
    axs[0].bar(x + width, rnaseq_vals, width, color="#d62728", edgecolor="black", label="Predicted RNA-seq Suppression (%)")

    for i in range(len(x)):
        axs[0].text(x[i], epistatic_vals[i] + 1.5, f"{epistatic_vals[i]:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#1f77b4")

    axs[0].set_title("A. In Silico CRISPR Base-Editing Rescue (Gynecomastia & ER+ Cancer)", fontweight="bold", fontsize=11)
    axs[0].set_ylabel("% of Hyperactive Phenotype Reversed (n=149 Mammary Tracks)")
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(panels, fontsize=9)
    axs[0].set_ylim(0, 110)
    axs[0].axhline(y=70.0, color="black", linestyle="--", alpha=0.5, label="Clinical Efficacy Threshold (70%)")
    axs[0].legend(loc="upper left", fontsize=8.5)

    # Panel B: Ventaja de Remodelación en Tejido Adiposo vs Epitelio (Reconstrucción Estética)
    modalities = df_selectividad["Modalidad_Epigenetica"].values
    adip_scores = df_selectividad["Ratio_AR_vs_ESR1_Adiposo_n25"].values
    mamm_scores = df_selectividad["Ratio_AR_vs_ESR1_Epitelio_n149"].values

    x_b = np.arange(len(modalities))
    width_b = 0.35

    axs[1].bar(x_b - width_b/2, mamm_scores, width_b, color="#306089", edgecolor="black", label="Mammary Epithelium/Stroma (n=149)")
    axs[1].bar(x_b + width_b/2, adip_scores, width_b, color="#388E3C", edgecolor="black", label="Adipose Tissue Tracks (n=25) [Aesthetic Target]")

    for i in range(len(x_b)):
        pct_gain = df_selectividad["Ventaja_Volumetrica_Adiposa_Pct"].values[i]
        if pct_gain > 0:
            axs[1].text(x_b[i] + width_b/2, adip_scores[i] + 0.04, f"+{pct_gain:.0f}%", ha="center", va="bottom", fontsize=9, fontweight="bold", color="#2e7d32")

    axs[1].axhline(y=1.0, color="red", linestyle="--", label="Neutral Equilibrium (Ratio = 1.0)")
    axs[1].set_title("B. Adipose vs. Epithelial Responsiveness for Tissue Engineering", fontweight="bold", fontsize=11)
    axs[1].set_ylabel("Median Impact Ratio (AR_loss / ESR1_gain)")
    axs[1].set_xticks(x_b)
    axs[1].set_xticklabels(modalities, rotation=15, fontsize=9)
    axs[1].set_ylim(0, 2.5)
    axs[1].legend(loc="upper right", fontsize=8.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"[4/4] Gráfica en inglés exportada exitosamente en: '{output_path}'")

def main():
    print("=" * 75)
    print("  SUITE DE VALIDACIÓN TRASLACIONAL IN SILICO - ALPHAGENOME (AR vs ESR1)")
    print("=" * 75)

    # 1. Ejecutar auditoría de pruebas
    ejecutar_pruebas_unitarias()

    # 2. Correr módulo de Supresión Clínica (Ginecomastia y Cáncer de Mama ER+)
    print("[2/4] Ejecutando simulaciones de rescate CRISPR en Ginecomastia y Oncología...")
    sim_clinico = ClinicalSuppressionSimulator()
    df_clinico = sim_clinico.run_all_clinical_panels()
    df_clinico.to_csv("resultados_simulacion_ginecomastia_oncologia.csv", index=False)
    print(df_clinico[[
        "Panel_Terapeutico", "Num_Guias_sgRNA", 
        "Rescate_Lineal_Pct", "Rescate_Epistatico_Haplotype_Pct", 
        "Supresion_Transcripcional_RNAseq_Pct"
    ]].to_string(index=False))
    print("\n      -> Guardado en 'resultados_simulacion_ginecomastia_oncologia.csv'\n")

    # 3. Correr módulo de Reconstrucción Estética e Ingeniería Tisular Adiposa
    print("[3/4] Evaluando selectividad adiposa y protocolos de reconstrucción estética...")
    sim_estetico = AestheticReconstructionSimulator()
    df_selectividad = sim_estetico.evaluate_adipose_selectivity()
    df_protocolos_esteticos = sim_estetico.simulate_regenerative_protocols()
    
    df_selectividad.to_csv("resultados_selectividad_adiposa_estetica.csv", index=False)
    df_protocolos_esteticos.to_csv("protocolos_bioingenieria_reconstruccion.csv", index=False)
    
    print(df_protocolos_esteticos[[
        "Protocolo_Bioingenieria", "Coordenadas_Diana", 
        "Apertura_Cromatina_Adiposa_DNASE", "Indice_Seguridad_Epitelial"
    ]].to_string(index=False))
    print("\n      -> Guardados en 'resultados_selectividad_adiposa_estetica.csv' y 'protocolos_bioingenieria_reconstruccion.csv'\n")

    # 4. Generar figura en inglés lista para el pitch
    generar_grafica_traslacional_en(
        df_clinico, 
        df_selectividad, 
        output_path="translational_in_silico_trials_EN.png"
    )
    print("=" * 75)

if __name__ == "__main__":
    main()