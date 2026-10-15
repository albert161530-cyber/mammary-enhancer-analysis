# main_pruebas_modelos.py
"""
ORQUESTADOR PRINCIPAL EN WINDOWS
Ejecuta GSEApy, ARCHS4 (Correlación + Fallback y k-NN), ENCODE REST API,
SCENIC+/CellOracle y gráficas en Windows, y coordina el worker pyBigWig de Ubuntu.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from prueba_gseapy_archs4 import TranscriptomicCloudValidator
from prueba_encode_bigwig import EncodeEpigenomicValidator
from prueba_scenic_celloracle import SingleCellGRNAdapter


def generar_grafica_validacion_apis(df_enrichr: pd.DataFrame, df_encode: pd.DataFrame, out_png: str):
    """Genera en Windows la figura resumen en inglés con los datos reales devueltos por las APIs."""
    fig, axs = plt.subplots(1, 2, figsize=(15.5, 6.2))
    fig.subplots_adjust(wspace=0.32)

    # Panel A: Resultados de Enrichr (GSEApy)
    if df_enrichr is not None and not df_enrichr.empty:
        top_terms = df_enrichr.head(8).iloc[::-1]
        terms = [t[:38] + "..." if len(t) > 38 else t for t in top_terms["Term"]]
        scores = -np.log10(top_terms["Adjusted P-value"].astype(float).clip(lower=1e-30))
        axs[0].barh(terms, scores, color="#1f77b4", edgecolor="black")
        axs[0].axvline(-np.log10(0.05), color="red", linestyle="--", label="FDR = 0.05")
        axs[0].set_xlabel("-log10(Adjusted P-value) [GSEApy / Enrichr]")
        axs[0].set_title("A. Functional Pathway Enrichment (AR/ESR1/CYP19A1 Network)", fontweight="bold", fontsize=10.5)
        axs[0].legend(loc="lower right")
    else:
        axs[0].text(0.5, 0.5, "No Enrichr terms retrieved", ha="center", va="center")
        axs[0].set_title("A. Functional Pathway Enrichment (GSEApy)", fontweight="bold")

    # Panel B: Disponibilidad de tracks experimentales en ENCODE (Epitelio vs Adiposo)
    if df_encode is not None and not df_encode.empty:
        piv = df_encode.groupby(["Assay_Type", "Clinical_Category"])["Released_Experiments_Found"].sum().unstack(fill_value=0)
        assays = piv.index.tolist()
        x = np.arange(len(assays))
        w = 0.35

        col_mamm = piv.get("Mammary_Oncology_Epithelium", pd.Series([0]*len(assays))).values
        col_adip = piv.get("Adipose_Stromal_Reconstruction", pd.Series([0]*len(assays))).values

        axs[1].bar(x - w/2, col_mamm, w, color="#306089", edgecolor="black", label="Mammary Epithelium / MCF-7 / T47D")
        axs[1].bar(x + w/2, col_adip, w, color="#388E3C", edgecolor="black", label="Adipose / Stromal Tissue")

        for i in range(len(x)):
            axs[1].text(x[i] - w/2, col_mamm[i] + 1, str(col_mamm[i]), ha="center", va="bottom", fontsize=8.5, fontweight="bold")
            axs[1].text(x[i] + w/2, col_adip[i] + 1, str(col_adip[i]), ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#2e7d32")

        axs[1].set_xticks(x)
        axs[1].set_xticklabels(assays, rotation=15)
        axs[1].set_ylim(0, max(max(col_mamm, default=10), max(col_adip, default=10)) * 1.18)
        axs[1].set_ylabel("Released ENCODE Experiments Queried")
        axs[1].set_title("B. ENCODE DCC API Experimental Coverage", fontweight="bold", fontsize=10.5)
        axs[1].legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print(f"  ✔ Figura exportada en Windows: '{out_png}'")


def main():
    print("=" * 78)
    print("  SUITE HÍBRIDA CORREGIDA (WINDOWS MAIN + UBUNTU PYBIGWIG WORKER)")
    print("=" * 78)

    # 1. GSEApy + ARCHS4 en Windows
    print("\n[1/4] Ejecutando en Windows: GSEApy (Enrichr API) y ARCHS4 Cloud Endpoints...")
    val_trans = TranscriptomicCloudValidator()

    df_enrichr = val_trans.run_gseapy_enrichr()
    if not df_enrichr.empty:
        df_enrichr.to_csv("validacion_1_gseapy_enrichr.csv", index=False)
        print("  ✔ Guardado: 'validacion_1_gseapy_enrichr.csv'")

    df_corr_cyp = val_trans.query_archs4_gene_correlation("CYP19A1", meta_filter="adipose", k=20)
    df_corr_esr1 = val_trans.query_archs4_gene_correlation("ESR1", meta_filter="breast", k=20)
    df_corr_all = pd.concat([df_corr_cyp, df_corr_esr1], ignore_index=True)
    if not df_corr_all.empty:
        df_corr_all.to_csv("validacion_2_archs4_correlaciones.csv", index=False)
        print("  ✔ Guardado: 'validacion_2_archs4_correlaciones.csv'")
        print(df_corr_all.head(6).to_string(index=False))

    knn_res = val_trans.query_archs4_knn_signature()
    if knn_res:
        pd.DataFrame({
            "Matching_GSM_Sample": knn_res.get("top_matching_gsm_samples", []),
            "Cosine_Distance": knn_res.get("top_distances", [])
        }).to_csv("validacion_3_archs4_knn_pacientes_gsm.csv", index=False)
        print("  ✔ Guardado: 'validacion_3_archs4_knn_pacientes_gsm.csv'")

    # 2. ENCODE REST API en Windows
    print("\n[2/4] Ejecutando en Windows: Inventario The ENCODE REST API (Ontologías Completas)...")
    val_encode = EncodeEpigenomicValidator(container_mode="manual")
    df_encode = val_encode.build_encode_inventory()
    if not df_encode.empty:
        df_encode.to_csv("validacion_4_encode_api_inventory.csv", index=False)
        print("  ✔ Guardado: 'validacion_4_encode_api_inventory.csv'")

    # 3. Señal BigWig vía Worker de Ubuntu (MCF-7 y Tejido Adiposo)
    print("\n[3/4] Obteniendo URLs de BigWig (GRCh38) en Windows para el Worker de Ubuntu...")
    for biosample_target, tag in [("MCF-7", "MCF7"), ("subcutaneous adipose tissue", "Adipose")]:
        bw_info = val_encode.get_released_bigwig_url(biosample_term=biosample_target)
        if bw_info:
            out_bw_csv = f"validacion_4b_pybigwig_{tag}_{bw_info['accession']}.csv"
            df_bw = val_encode.run_bigwig_via_ubuntu_container(bw_info["url"], out_bw_csv, window_bp=150)
            if not df_bw.empty:
                print(df_bw.to_string(index=False))

    # 4. Tablas SCENIC+ y CellOracle + Gráfica final en Windows
    print("\n[4/4] Generando tablas SCENIC+/CellOracle y figura final en Windows...")
    sc_adapter = SingleCellGRNAdapter(enhancer_half_width=250)
    sc_adapter.export_scenicplus_eregulon_metadata().to_csv("validacion_5_scenicplus_eregulons_hg38.csv", index=False)
    sc_adapter.export_celloracle_base_grn_bed().to_csv("validacion_6_celloracle_base_grn.csv", index=False)
    print("  ✔ Guardados: 'validacion_5_scenicplus_eregulons_hg38.csv' y 'validacion_6_celloracle_base_grn.csv'")

    generar_grafica_validacion_apis(
        df_enrichr, df_encode,
        out_png="standard_models_cloud_validation_EN.png"
    )
    print("\n" + "=" * 78)
    print("  PROCESO EN WINDOWS FINALIZADO.")
    print("=" * 78)


if __name__ == "__main__":
    main()