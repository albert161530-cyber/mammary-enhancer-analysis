# prueba_scenic_celloracle.py
"""
Módulo 3 (Windows): Adaptador oficial para SCENIC+ (eRegulon) y CellOracle (base_GRN).
"""

import pandas as pd
from config_modelos import DRIVER_VARIANTS


class SingleCellGRNAdapter:
    def __init__(self, enhancer_half_width: int = 250):
        self.half_width = enhancer_half_width

    def export_scenicplus_eregulon_metadata(self) -> pd.DataFrame:
        rows = []
        for d in DRIVER_VARIANTS:
            start = d.pos - self.half_width
            end = d.pos + self.half_width
            scenic_region_name = f"{d.chrom}:{start}-{end}"

            rows.append({
                "Region_Name_SCENIC": scenic_region_name,
                "Chromosome": d.chrom,
                "Start": start,
                "End": end,
                "Variant_ID": d.variant_id,
                "TF_Regulator": "AR",
                "Target_Gene": d.gene,
                "TF_to_Gene_Relationship": "negative (-)",
                "Region_to_Gene_Relationship": "positive (+)",
                "Evolutionary_State_hg38": "AR_loss_derepressed",
                "AlphaGenome_L1_Weight": d.l1_score
            })
        return pd.DataFrame(rows)

    def export_celloracle_base_grn_bed(self) -> pd.DataFrame:
        rows = []
        for d in DRIVER_VARIANTS:
            start = d.pos - self.half_width
            end = d.pos + self.half_width
            peak_id = f"{d.chrom}_{start}_{end}"
            rows.append({
                "peak_id": peak_id,
                "gene_short_name": d.gene,
                "AR_motif_ancestral": 1,
                "AR_motif_human_hg38": 0,
                "AR_motif_CRISPR_rescue": 1,
                "ESR1_motif": 1 if d.gene in ["ESR1", "CYP19A1"] else 0,
                "FOXA1_pioneer_motif": 1
            })
        return pd.DataFrame(rows)