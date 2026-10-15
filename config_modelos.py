# config_modelos.py
"""
Configuración central de coordenadas hg38, firmas génicas y endpoints de APIs
oficiales (ENCODE, ARCHS4, Enrichr/GSEApy y GEPIA2).
"""

from dataclasses import dataclass
from typing import List, Dict

@dataclass(frozen=True)
class TargetLocus:
    variant_id: str
    chrom: str
    pos: int
    ref_human: str
    alt_ancestral: str
    gene: str
    l1_score: float

# Tus 10 variantes Driver principales de AR_loss
DRIVER_VARIANTS: List[TargetLocus] = [
    TargetLocus("chr6:150866313:T>C", "chr6", 150866313, "C", "T", "ESR1", 419.5),
    TargetLocus("chr6:151452006:A>G", "chr6", 151452006, "G", "A", "ESR1", 338.0),
    TargetLocus("chr6:151808589:A>G", "chr6", 151808589, "G", "A", "ESR1", 253.7),
    TargetLocus("chr15:51569074:A>G", "chr15", 51569074, "G", "A", "CYP19A1", 240.6),
    TargetLocus("chr6:150866751:T>C", "chr6", 150866751, "C", "T", "ESR1", 235.8),
    TargetLocus("chr14:64280081:A>G", "chr14", 64280081, "G", "A", "ESR2", 225.9),
    TargetLocus("chr6:151299349:A>G", "chr6", 151299349, "G", "A", "ESR1", 220.6),
    TargetLocus("chr6:151014628:A>G", "chr6", 151014628, "G", "A", "ESR1", 215.3),
    TargetLocus("chr6:151322799:T>C", "chr6", 151322799, "C", "T", "ESR1", 210.1),
    TargetLocus("chr6:151763572:G>C", "chr6", 151763572, "C", "G", "ESR1", 206.6),
]

# Red regulatoria para GSEApy, ARCHS4 y GEPIA2
GENE_NETWORK_SIGNATURE: Dict[str, List[str]] = {
    "CORE_TARGETS": ["ESR1", "ESR2", "CYP19A1", "AR"],
    "PIONEER_COFACTORS": ["FOXA1", "GATA3", "EP300", "NCOA1", "NCOA3", "PGR"],
    "UP_IN_BROKEN_BRAKE": ["ESR1", "ESR2", "CYP19A1", "GREB1", "TFF1", "CCND1", "MYC", "PGR"],
    "DOWN_IN_BROKEN_BRAKE": ["AR", "KLK3", "FKBP5", "TMPRSS2", "ZBTB16", "PTEN"]
}

# Endpoints oficiales de ARCHS4
ARCHS4_API_URLS: Dict[str, str] = {
    "correlation": "https://maayanlab.cloud/sigpy/data/correlation",
    "knn_signature": "https://maayanlab.cloud/sigpy/data/knn/signature",
    "diffexp": "https://maayanlab.cloud/sigpy/data/diffexp"
}

ENCODE_SEARCH_URL: str = "https://www.encodeproject.org/search/"

# Términos ontológicos oficiales en ENCODE para capturar todas las pistas disponibles
ENCODE_BIOSAMPLES: Dict[str, List[str]] = {
    "Mammary_Oncology_Epithelium": [
        "MCF-7", 
        "T47D", 
        "mammary epithelial cell", 
        "breast epithelium"
    ],
    "Adipose_Stromal_Reconstruction": [
        "adipose tissue", 
        "adipocyte", 
        "subcutaneous adipose tissue"
    ]
}