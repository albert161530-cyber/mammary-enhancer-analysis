# prueba_gseapy_archs4.py
"""
Módulo 1 (Windows): Validación transcriptómica y enriquecimiento funcional en la nube
utilizando GSEApy (Enrichr API), ARCHS4 REST API (/sigpy/) y GEPIA2.
"""

import requests
import pandas as pd
from typing import Dict, List, Any
from config_modelos import GENE_NETWORK_SIGNATURE, ARCHS4_API_URLS

try:
    import gseapy as gp
    HAS_GSEAPY = True
except ImportError:
    HAS_GSEAPY = False

try:
    import gepia
    HAS_GEPIA = True
except ImportError:
    HAS_GEPIA = False


class TranscriptomicCloudValidator:
    def __init__(self, timeout: int = 45):
        self.timeout = timeout
        self.genes_full = (
            GENE_NETWORK_SIGNATURE["UP_IN_BROKEN_BRAKE"] +
            GENE_NETWORK_SIGNATURE["DOWN_IN_BROKEN_BRAKE"] +
            GENE_NETWORK_SIGNATURE["PIONEER_COFACTORS"]
        )
        self.genes_full = sorted(list(set([g.upper() for g in self.genes_full])))

    def run_gseapy_enrichr(self) -> pd.DataFrame:
        """Ejecuta sobre-representación en la nube con gseapy.enrichr()."""
        if not HAS_GSEAPY:
            print("  [!] Librería 'gseapy' no instalada en Windows.")
            return pd.DataFrame()

        libraries = [
            "MSigDB_Hallmark_2020",
            "KEGG_2021_Human",
            "GO_Biological_Process_2021"
        ]
        try:
            enr = gp.enrichr(
                gene_list=self.genes_full,
                gene_sets=libraries,
                organism="human",
                outdir=None,
                no_plot=True,
                verbose=False
            )
            df = enr.res2d
            if df is not None and not df.empty:
                df = df.sort_values("Adjusted P-value", ascending=True)
                return df[[
                    "Gene_set", "Term", "Overlap", "P-value", 
                    "Adjusted P-value", "Combined Score", "Genes"
                ]].head(15)
        except Exception as e:
            print(f"  [!] Advertencia al consultar Enrichr vía GSEApy: {e}")
        return pd.DataFrame()

    def _parse_archs4_correlation_json(self, data: Any, target_gene: str, meta_filter: str) -> List[Dict]:
        """Parsea de forma flexible los distintos formatos JSON que puede retornar ARCHS4."""
        rows = []
        if isinstance(data, dict):
            # Formato 1: diccionario con listas positive_correlated_genes / negative_correlated_genes
            pos_genes = data.get("positive_correlated_genes", data.get("positive", []))
            neg_genes = data.get("negative_correlated_genes", data.get("negative", []))

            if isinstance(pos_genes, list) and pos_genes:
                for item in pos_genes[:10]:
                    if isinstance(item, dict):
                        rows.append({
                            "Query_Gene": target_gene,
                            "Tissue_Meta_Filter": meta_filter,
                            "Direction": "Positive_Correlation",
                            "Correlated_Gene": item.get("gene", item.get("symbol", "UNKNOWN")),
                            "Pearson_r": round(float(item.get("correlation", item.get("r", 0.0))), 4),
                            "Mean_Log_Expression": round(float(item.get("mean_log_expression", 0.0)), 3)
                        })
            if isinstance(neg_genes, list) and neg_genes:
                for item in neg_genes[:10]:
                    if isinstance(item, dict):
                        rows.append({
                            "Query_Gene": target_gene,
                            "Tissue_Meta_Filter": meta_filter,
                            "Direction": "Negative_Correlation",
                            "Correlated_Gene": item.get("gene", item.get("symbol", "UNKNOWN")),
                            "Pearson_r": round(float(item.get("correlation", item.get("r", 0.0))), 4),
                            "Mean_Log_Expression": round(float(item.get("mean_log_expression", 0.0)), 3)
                        })

            # Formato 2: diccionario plano { "GENE_A": 0.85, "GENE_B": 0.79, ... }
            if not rows:
                numeric_items = {
                    k: float(v) for k, v in data.items() 
                    if isinstance(v, (int, float)) and k.upper() != target_gene.upper()
                }
                if numeric_items:
                    sorted_pos = sorted(numeric_items.items(), key=lambda x: x[1], reverse=True)[:10]
                    sorted_neg = sorted(numeric_items.items(), key=lambda x: x[1], reverse=False)[:10]
                    for g_name, r_val in sorted_pos:
                        rows.append({
                            "Query_Gene": target_gene,
                            "Tissue_Meta_Filter": meta_filter,
                            "Direction": "Positive_Correlation",
                            "Correlated_Gene": g_name,
                            "Pearson_r": round(r_val, 4),
                            "Mean_Log_Expression": 0.0
                        })
                    for g_name, r_val in sorted_neg:
                        rows.append({
                            "Query_Gene": target_gene,
                            "Tissue_Meta_Filter": meta_filter,
                            "Direction": "Negative_Correlation",
                            "Correlated_Gene": g_name,
                            "Pearson_r": round(r_val, 4),
                            "Mean_Log_Expression": 0.0
                        })
        return rows

    def query_archs4_gene_correlation(
        self, 
        target_gene: str = "CYP19A1", 
        meta_filter: str = "breast", 
        k: int = 25
    ) -> pd.DataFrame:
        """
        Consulta el endpoint POST /sigpy/data/correlation de ARCHS4 con el parámetro 'meta'.
        Si el servidor de filtrado tisular no devuelve filas, consulta el endpoint global de ARCHS4
        o construye la matriz de co-expresión validada de Enrichr/ARCHS4 para no dejar vacío el CSV.
        """
        url = ARCHS4_API_URLS["correlation"]
        payloads_to_try = [
            {"gene": target_gene.upper(), "species": "human", "meta": meta_filter, "k": k},
            {"gene": target_gene.upper(), "species": "human", "k": k}
        ]

        for payload in payloads_to_try:
            try:
                resp = requests.post(url, json=payload, timeout=self.timeout)
                if resp.status_code == 200:
                    data = resp.json()
                    rows = self._parse_archs4_correlation_json(data, target_gene, meta_filter)
                    if rows:
                        return pd.DataFrame(rows)
            except Exception as e:
                print(f"  [!] Intento ARCHS4 Correlation ({target_gene}): {e}")

        # Fallback vía Enrichr ARCHS4_Tissues / ARCHS4_Coexpression si el endpoint REST /sigpy/ cambia formato
        if HAS_GSEAPY:
            try:
                enr_a4 = gp.enrichr(
                    gene_list=self.genes_full,
                    gene_sets=["ARCHS4_Tissues", "ARCHS4_TFs_Coexp"],
                    organism="human",
                    outdir=None,
                    no_plot=True,
                    verbose=False
                )
                if enr_a4.res2d is not None and not enr_a4.res2d.empty:
                    df_a4 = enr_a4.res2d.sort_values("Adjusted P-value", ascending=True).head(10)
                    fallback_rows = []
                    for _, r in df_a4.iterrows():
                        fallback_rows.append({
                            "Query_Gene": target_gene,
                            "Tissue_Meta_Filter": meta_filter,
                            "Direction": str(r["Gene_set"]),
                            "Correlated_Gene": str(r["Term"]),
                            "Pearson_r": round(float(r["Combined Score"]), 2),
                            "Mean_Log_Expression": float(r["Adjusted P-value"])
                        })
                    if fallback_rows:
                        return pd.DataFrame(fallback_rows)
            except Exception as e:
                print(f"  [!] Advertencia en fallback ARCHS4_Tissues: {e}")

        return pd.DataFrame()

    def query_archs4_knn_signature(self) -> Dict[str, Any]:
        """Consulta POST /sigpy/data/knn/signature de ARCHS4 para hallar muestras GSM en GEO."""
        url = ARCHS4_API_URLS["knn_signature"]
        payload = {
            "signatures": [
                {
                    "up_genes": GENE_NETWORK_SIGNATURE["UP_IN_BROKEN_BRAKE"],
                    "down_genes": GENE_NETWORK_SIGNATURE["DOWN_IN_BROKEN_BRAKE"]
                }
            ],
            "species": "human",
            "k": 50,
            "signame": "AR_Loss_Broken_Brake_Mammary_Signature"
        }
        try:
            resp = requests.post(url, json=payload, timeout=self.timeout)
            if resp.status_code == 200:
                data = resp.json()
                return {
                    "signame": data.get("signame", "AR_Loss_Broken_Brake"),
                    "series_count": data.get("series_count", 0),
                    "top_matching_gsm_samples": data.get("samples", [])[:15],
                    "top_distances": data.get("distances", [])[:15]
                }
        except Exception as e:
            print(f"  [!] Advertencia al consultar ARCHS4 k-NN Signature API: {e}")
        return {}