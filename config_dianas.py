# config_dianas.py
"""
Configuración central de coordenadas genómicas (hg38), impactos L1 base y 
parámetros de tejidos obtenidos del pipeline AlphaGenome.
"""

from dataclasses import dataclass
from typing import List, Dict

@dataclass(frozen=True)
class DriverVariant:
    variant_id: str
    chrom: str
    pos: int
    ref_human: str      # Estado actual en hg38 (freno roto)
    alt_ancestral: str  # Estado ancestral (freno AR funcional)
    locus: str
    mean_l1_impact: float

# Top 10 variantes 'Driver' de AR_loss (Panel D)
TOP_AR_LOSS_DRIVERS: List[DriverVariant] = [
    DriverVariant("chr6:150866313:T>C", "chr6", 150866313, "C", "T", "ESR1", 419.5),
    DriverVariant("chr6:151452006:A>G", "chr6", 151452006, "G", "A", "ESR1", 338.0),
    DriverVariant("chr6:151808589:A>G", "chr6", 151808589, "G", "A", "ESR1", 253.7),
    DriverVariant("chr15:51569074:A>G", "chr15", 51569074, "G", "A", "CYP19A1", 240.6),
    DriverVariant("chr6:150866751:T>C", "chr6", 150866751, "C", "T", "ESR1", 235.8),
    DriverVariant("chr14:64280081:A>G", "chr14", 64280081, "G", "A", "ESR2", 225.9),
    DriverVariant("chr6:151299349:A>G", "chr6", 151299349, "G", "A", "ESR1", 220.6),
    DriverVariant("chr6:151014628:A>G", "chr6", 151014628, "G", "A", "ESR1", 215.3),
    DriverVariant("chr6:151322799:T>C", "chr6", 151322799, "C", "T", "ESR1", 210.1),
    DriverVariant("chr6:151763572:G>C", "chr6", 151763572, "C", "G", "ESR1", 206.6),
]

# Factores de sinergia epistática por locus (Panel C: Joint vs Ind. Sum)
EPISTATIC_SYNERGY_MULTIPLIERS: Dict[str, float] = {
    "ESR1": 4.0,     # 4034 vs 1020
    "ESR2": 5.5,     # 3143 vs 571
    "CYP19A1": 5.1,  # 1651 vs 322
    "AR": 3.1        # 556 vs 178
}

# Ratios de magnitud de efecto por tejido y modalidad (Panel B)
TISSUE_MODALITY_RATIOS: Dict[str, Dict[str, float]] = {
    "Mammary_Epithelium_Stroma": { # n=149 tracks
        "ATAC": 1.43, "DNASE": 1.15, "CHIP_HISTONE": 1.71, 
        "CHIP_TF": 1.48, "CAGE": 1.44, "RNA_SEQ": 1.11
    },
    "Adipose_Tissue": {            # n=25 tracks
        "ATAC": 1.81, "DNASE": 2.03, "CHIP_HISTONE": 1.70, 
        "CHIP_TF": 2.20, "CAGE": 1.52, "RNA_SEQ": 1.50
    },
    "Other_Endocrine": {           # n=101 tracks
        "ATAC": 1.96, "DNASE": 1.98, "CHIP_HISTONE": 2.07, 
        "CHIP_TF": 2.10, "CAGE": 1.51, "RNA_SEQ": 1.31
    }
}

# Impacto total acumulado de referencia por locus en tejido mamario/adiposo
LOCUS_BASELINE_BURDEN: Dict[str, float] = {
    "ESR1": 4034.0,
    "ESR2": 3143.0,
    "CYP19A1": 1651.0,
    "AR": 556.0
}