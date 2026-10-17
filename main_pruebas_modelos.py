# main_pruebas_modelos.py
"""
ORQUESTADOR PRINCIPAL EN WINDOWS
Ejecuta GSEApy, ARCHS4, ENCODE REST API, SCENIC+/CellOracle y la gráfica.
El análisis pyBigWig se ejecuta en Ubuntu mediante WSL.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from prueba_gseapy_archs4 import TranscriptomicCloudValidator
from prueba_encode_bigwig import EncodeEpigenomicValidator
from prueba_scenic_celloracle import SingleCellGRNAdapter

BASE_DIR = Path(__file__).resolve().parent
WSL_DISTRIBUTION = "Ubuntu"
WSL_PYTHON = "/home/precision_5560/encode-env/bin/python"
WORKER_TIMEOUT = 600


def guardar_csv(df: pd.DataFrame, nombre: str) -> None:
    ruta = BASE_DIR / nombre
    df.to_csv(ruta, index=False)
    print(f"  ✔ Guardado: '{ruta.name}'", flush=True)


def generar_grafica_validacion_apis(
    df_enrichr: pd.DataFrame,
    df_encode: pd.DataFrame,
    out_png: str,
) -> None:
    """Genera la figura resumen en inglés con los datos disponibles."""
    fig, axs = plt.subplots(1, 2, figsize=(15.5, 6.2))
    fig.subplots_adjust(wspace=0.32)

    if (
        df_enrichr is not None
        and not df_enrichr.empty
        and {"Term", "Adjusted P-value"}.issubset(df_enrichr.columns)
    ):
        top_terms = df_enrichr.head(8).iloc[::-1]
        terms = [
            str(t)[:38] + "..." if len(str(t)) > 38 else str(t)
            for t in top_terms["Term"]
        ]
        pvals = pd.to_numeric(
            top_terms["Adjusted P-value"], errors="coerce"
        ).fillna(1.0).clip(lower=1e-30)
        scores = -np.log10(pvals)
        axs[0].barh(terms, scores, color="#1f77b4", edgecolor="black")
        axs[0].axvline(
            -np.log10(0.05), color="red", linestyle="--", label="FDR = 0.05"
        )
        axs[0].set_xlabel("-log10(Adjusted P-value) [GSEApy / Enrichr]")
        axs[0].set_title(
            "A. Functional Pathway Enrichment (AR/ESR1/CYP19A1 Network)",
            fontweight="bold", fontsize=10.5,
        )
        axs[0].legend(loc="lower right")
    else:
        axs[0].text(
            0.5, 0.5, "No Enrichr terms retrieved",
            ha="center", va="center", transform=axs[0].transAxes,
        )
        axs[0].set_title("A. Functional Pathway Enrichment (GSEApy)", fontweight="bold")

    required = {"Assay_Type", "Clinical_Category", "Released_Experiments_Found"}
    if df_encode is not None and not df_encode.empty and required.issubset(df_encode.columns):
        piv = (
            df_encode.groupby(["Assay_Type", "Clinical_Category"])[
                "Released_Experiments_Found"
            ].sum().unstack(fill_value=0)
        )
        assays = piv.index.tolist()
        x = np.arange(len(assays))
        w = 0.35
        col_mamm = piv.get(
            "Mammary_Oncology_Epithelium", pd.Series(0, index=assays)
        ).to_numpy()
        col_adip = piv.get(
            "Adipose_Stromal_Reconstruction", pd.Series(0, index=assays)
        ).to_numpy()

        axs[1].bar(
            x - w / 2, col_mamm, w, color="#306089", edgecolor="black",
            label="Mammary Epithelium / MCF-7 / T47D",
        )
        axs[1].bar(
            x + w / 2, col_adip, w, color="#388E3C", edgecolor="black",
            label="Adipose / Stromal Tissue",
        )
        for i in range(len(x)):
            axs[1].text(
                x[i] - w / 2, col_mamm[i] + 1, str(int(col_mamm[i])),
                ha="center", va="bottom", fontsize=8.5, fontweight="bold",
            )
            axs[1].text(
                x[i] + w / 2, col_adip[i] + 1, str(int(col_adip[i])),
                ha="center", va="bottom", fontsize=8.5, fontweight="bold",
                color="#2e7d32",
            )
        axs[1].set_xticks(x)
        axs[1].set_xticklabels(assays, rotation=15)
        ymax = max(
            float(np.max(col_mamm)) if len(col_mamm) else 0.0,
            float(np.max(col_adip)) if len(col_adip) else 0.0,
            1.0,
        )
        axs[1].set_ylim(0, ymax * 1.18 + 1)
        axs[1].set_ylabel("Released ENCODE Experiments Queried")
        axs[1].set_title(
            "B. ENCODE DCC API Experimental Coverage", fontweight="bold", fontsize=10.5
        )
        axs[1].legend(loc="upper right")
    else:
        axs[1].text(
            0.5, 0.5, "No ENCODE inventory data retrieved",
            ha="center", va="center", transform=axs[1].transAxes,
        )
        axs[1].set_title("B. ENCODE DCC API Experimental Coverage", fontweight="bold")

    plt.tight_layout()
    out_path = BASE_DIR / out_png
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔ Figura exportada en Windows: '{out_path.name}'", flush=True)


def main() -> int:
    print("=" * 78, flush=True)
    print("  SUITE HÍBRIDA (WINDOWS MAIN + UBUNTU/WSL PYBIGWIG WORKER)", flush=True)
    print("=" * 78, flush=True)

    df_enrichr = pd.DataFrame()
    df_encode = pd.DataFrame()

    # 1. GSEApy + ARCHS4
    print(
        "\n[1/4] Ejecutando en Windows: GSEApy (Enrichr API) y ARCHS4 Cloud Endpoints...",
        flush=True,
    )
    try:
        val_trans = TranscriptomicCloudValidator()
        df_enrichr = val_trans.run_gseapy_enrichr()
        if df_enrichr is not None and not df_enrichr.empty:
            guardar_csv(df_enrichr, "validacion_1_gseapy_enrichr.csv")
        else:
            print("  [AVISO] GSEApy no devolvió resultados.", flush=True)

        df_corr_cyp = val_trans.query_archs4_gene_correlation(
            "CYP19A1", meta_filter="adipose", k=20
        )
        df_corr_esr1 = val_trans.query_archs4_gene_correlation(
            "ESR1", meta_filter="breast", k=20
        )
        df_corr_all = pd.concat([df_corr_cyp, df_corr_esr1], ignore_index=True)
        if not df_corr_all.empty:
            guardar_csv(df_corr_all, "validacion_2_archs4_correlaciones.csv")
            print(df_corr_all.head(6).to_string(index=False), flush=True)
        else:
            print("  [AVISO] ARCHS4 no devolvió correlaciones.", flush=True)

        knn_res = val_trans.query_archs4_knn_signature()
        if knn_res:
            guardar_csv(
                pd.DataFrame({
                    "Matching_GSM_Sample": knn_res.get("top_matching_gsm_samples", []),
                    "Cosine_Distance": knn_res.get("top_distances", []),
                }),
                "validacion_3_archs4_knn_pacientes_gsm.csv",
            )
        else:
            print("  [AVISO] ARCHS4 k-NN no devolvió resultados.", flush=True)
    except Exception as exc:
        print(f"  [ERROR] Bloque GSEApy/ARCHS4: {type(exc).__name__}: {exc}", flush=True)

    # 2. ENCODE REST API. Se configura WSL, no Docker.
    print("\n[2/4] Ejecutando en Windows: Inventario ENCODE REST API...", flush=True)
    val_encode = EncodeEpigenomicValidator(
        container_mode="wsl",
        wsl_distribution=WSL_DISTRIBUTION,
        conda_env_Or_python=WSL_PYTHON,
    )
    try:
        df_encode = val_encode.build_encode_inventory()
        if df_encode is not None and not df_encode.empty:
            guardar_csv(df_encode, "validacion_4_encode_api_inventory.csv")
        else:
            print("  [AVISO] El inventario ENCODE está vacío.", flush=True)
    except Exception as exc:
        print(f"  [ERROR] Inventario ENCODE: {type(exc).__name__}: {exc}", flush=True)

    # 3. BigWig en Ubuntu/WSL
    print(
        "\n[3/4] Obteniendo URLs BigWig (GRCh38) para el worker de Ubuntu/WSL...",
        flush=True,
    )
    for biosample_target, tag in [
        ("MCF-7", "MCF7"),
        ("subcutaneous adipose tissue", "Adipose"),
    ]:
        try:
            print(f"\n  Buscando BigWig para: {biosample_target}", flush=True)
            bw_info = val_encode.get_released_bigwig_url(
                biosample_term=biosample_target,
                assembly="GRCh38",
            )
            if not bw_info or not bw_info.get("url"):
                print(
                    f"  [AVISO] ENCODE no devolvió una URL BigWig para {biosample_target}.",
                    flush=True,
                )
                continue

            accession = bw_info.get("accession", "unknown")
            out_bw_csv = BASE_DIR / f"validacion_4b_pybigwig_{tag}_{accession}.csv"
            print(
                f"  Archivo ENCODE: {accession}\n"
                f"  Ensayo: {bw_info.get('assay_title', '')}\n"
                f"  URL: {bw_info['url']}",
                flush=True,
            )
            df_bw = val_encode.run_bigwig_via_ubuntu_container(
                bigwig_url=bw_info["url"],
                output_csv=str(out_bw_csv),
                window_bp=150,
                force_refresh=True,
                worker_timeout=WORKER_TIMEOUT,
            )
            if df_bw is not None and not df_bw.empty:
                print(df_bw.to_string(index=False), flush=True)
            else:
                print(f"  [AVISO] No se obtuvieron señales para {biosample_target}.", flush=True)
        except Exception as exc:
            print(
                f"  [ERROR] BigWig para {biosample_target}: {type(exc).__name__}: {exc}",
                flush=True,
            )

    try:
        val_encode.close()
    except Exception:
        pass

    # 4. SCENIC+, CellOracle y gráfica
    print(
        "\n[4/4] Generando tablas SCENIC+/CellOracle y figura final en Windows...",
        flush=True,
    )
    try:
        sc_adapter = SingleCellGRNAdapter(enhancer_half_width=250)
        guardar_csv(
            sc_adapter.export_scenicplus_eregulon_metadata(),
            "validacion_5_scenicplus_eregulons_hg38.csv",
        )
        guardar_csv(
            sc_adapter.export_celloracle_base_grn_bed(),
            "validacion_6_celloracle_base_grn.csv",
        )
    except Exception as exc:
        print(f"  [ERROR] SCENIC+/CellOracle: {type(exc).__name__}: {exc}", flush=True)

    try:
        generar_grafica_validacion_apis(
            df_enrichr, df_encode, "standard_models_cloud_validation_EN.png"
        )
    except Exception as exc:
        print(f"  [ERROR] Gráfica: {type(exc).__name__}: {exc}", flush=True)

    print("\n" + "=" * 78, flush=True)
    print("  PROCESO EN WINDOWS FINALIZADO.", flush=True)
    print("=" * 78, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
