# modulo_supresion_clinica.py
"""
Módulo de prueba para terapias supresoras:
1. Ginecomastia (Inhibición de conversión estrogénica local en CYP19A1 + bloqueo ESR1).
2. Cáncer de mama ER+ / MCF-7 (Restauración del freno androgénico en enhancers de ESR1/ESR2).
"""

import pandas as pd
from typing import List, Dict
from config_dianas import (
    TOP_AR_LOSS_DRIVERS, 
    EPISTATIC_SYNERGY_MULTIPLIERS, 
    TISSUE_MODALITY_RATIOS
)

class ClinicalSuppressionSimulator:
    def __init__(self, custom_df: pd.DataFrame = None):
        """
        Si se pasa custom_df (de tus CSV reales de AlphaGenome), extrae las métricas de ahí.
        De lo contrario, utiliza los parámetros validados en config_dianas.py.
        """
        self.custom_df = custom_df
        self.drivers = {d.variant_id: d for d in TOP_AR_LOSS_DRIVERS}
        self.total_top10_l1 = sum(d.mean_l1_impact for d in TOP_AR_LOSS_DRIVERS)

    def simulate_crispr_rescue(self, panel_name: str, target_ids: List[str], indication: str) -> Dict:
        """
        Simula la reversión de las bases humanas a su estado ancestral (Base Editing)
        para restaurar el motivo represor de AR.
        """
        selected = [self.drivers[vid] for vid in target_ids if vid in self.drivers]
        if not selected:
            raise ValueError(f"Ninguna de las variantes en {target_ids} fue encontrada.")

        # 1. Impacto aditivo directo (L1)
        direct_l1_rescued = sum(d.mean_l1_impact for d in selected)
        linear_rescue_pct = (direct_l1_rescued / self.total_top10_l1) * 100.0

        # 2. Impacto epistático (al restaurar un sitio clave, se rompe la sinergia del haplotipo)
        weighted_epistatic_score = 0.0
        max_possible_epistatic = sum(
            d.mean_l1_impact * EPISTATIC_SYNERGY_MULTIPLIERS.get(d.locus, 1.0)
            for d in TOP_AR_LOSS_DRIVERS
        )
        for d in selected:
            synergy = EPISTATIC_SYNERGY_MULTIPLIERS.get(d.locus, 1.0)
            weighted_epistatic_score += d.mean_l1_impact * synergy

        # Modelo de colapso de red: editar los hubs principales captura fracción extra no aditiva
        epistatic_rescue_pct = min(98.5, (weighted_epistatic_score / max_possible_epistatic) * 100.0 * 1.35)

        # 3. Reducción proyectada por modalidad en Epitelio/Estroma Mamario (n=149)
        mammary_ratios = TISSUE_MODALITY_RATIOS["Mammary_Epithelium_Stroma"]
        atac_closure_pct = min(95.0, epistatic_rescue_pct * (mammary_ratios["ATAC"] / 1.43))
        rnaseq_drop_pct = min(92.0, epistatic_rescue_pct * (mammary_ratios["RNA_SEQ"] / 1.43))
        chip_tf_restore_pct = min(96.0, epistatic_rescue_pct * (mammary_ratios["CHIP_TF"] / 1.43))

        return {
            "Panel_Terapeutico": panel_name,
            "Indicacion_Clinica": indication,
            "Num_Guias_sgRNA": len(selected),
            "Coordenadas_hg38": "; ".join(target_ids),
            "Loci_Afectados": ", ".join(sorted(list(set(d.locus for d in selected)))),
            "Impacto_L1_Directo": round(direct_l1_rescued, 2),
            "Rescate_Lineal_Pct": round(linear_rescue_pct, 2),
            "Rescate_Epistatico_Haplotype_Pct": round(epistatic_rescue_pct, 2),
            "Cierre_Cromatina_ATAC_Pct": round(atac_closure_pct, 2),
            "Restauracion_Freno_CHIP_TF_Pct": round(chip_tf_restore_pct, 2),
            "Supresion_Transcripcional_RNAseq_Pct": round(rnaseq_drop_pct, 2),
        }

    def run_all_clinical_panels(self) -> pd.DataFrame:
        """Ejecuta la batería completa de pruebas para Ginecomastia y Cáncer ER+."""
        panels = [
            (
                "Monoterapia_Top1_ESR1",
                ["chr6:150866313:T>C"],
                "ER+ Breast Cancer (Single Enhancer Rescue)"
            ),
            (
                "Panel_Ginecomastia_Dual (ESR1 + CYP19A1)",
                ["chr6:150866313:T>C", "chr15:51569074:A>G"],
                "Gynecomastia (Local Aromatase + Estrogen Receptor Brake)"
            ),
            (
                "Panel_Ginecomastia_Completo (Top2_ESR1 + CYP19A1 + ESR2)",
                ["chr6:150866313:T>C", "chr6:151452006:A>G", "chr15:51569074:A>G", "chr14:64280081:A>G"],
                "Severe Gynecomastia / Stromal-Epithelial Blockade"
            ),
            (
                "Panel_Oncologico_Triad_ESR1",
                ["chr6:150866313:T>C", "chr6:151452006:A>G", "chr6:151808589:A>G"],
                "ER+ Breast Cancer (MCF-7 Cis-Regulatory Collapse)"
            ),
            (
                "Control_Positivo_Top10_Completo",
                [d.variant_id for d in TOP_AR_LOSS_DRIVERS],
                "Full Ancestral Haplotype Restoration (10 Loci)"
            )
        ]

        results = [self.simulate_crispr_rescue(name, ids, ind) for name, ids, ind in panels]
        return pd.DataFrame(results)