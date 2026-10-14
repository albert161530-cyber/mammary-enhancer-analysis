#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analisis_posthoc_mamario.py
===========================
Pipeline de auditoría estadística post-hoc, filtrado tejido-específico (Mamario/Adiposo),
cuantificación de epistasis (no-aditividad), pruebas exactas binomiales y de Wilcoxon,
y generación de figuras de calidad de publicación a partir de las salidas de AlphaGenome.

Entradas requeridas (generadas por automatizar_bloques.py):
  - resultados_impacto_conjunto_regiones.csv
  - resultados_impacto_conjunto_resumen_modalidad.csv
  - resultados_impacto_conjunto_comparacion.csv
  - resultados_impacto_conjunto_por_track.csv
  - resultados_impacto_individual_por_track.csv

Salidas reproducibles generadas por este script:
  - tabla_1_resumen_global_y_mamario_por_modalidad.csv
  - tabla_2_epistasis_no_aditividad_por_locus.csv
  - tabla_3_top_variantes_drivers_mamario_adiposo.csv
  - tabla_4_cross_talk_factores_directos_MCF7.csv
  - resumen_validacion_alphagenome_mamario.png (300 DPI)
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon
import matplotlib.pyplot as plt

# Configuración de rutas
CARPETA = Path(__file__).resolve().parent

ARCHIVO_REGIONES = CARPETA / "resultados_impacto_conjunto_regiones.csv"
ARCHIVO_MODALIDAD = CARPETA / "resultados_impacto_conjunto_resumen_modalidad.csv"
ARCHIVO_COMPARACION = CARPETA / "resultados_impacto_conjunto_comparacion.csv"
ARCHIVO_CONJUNTO = CARPETA / "resultados_impacto_conjunto_por_track.csv"
ARCHIVO_INDIVIDUAL = CARPETA / "resultados_impacto_individual_por_track.csv"

LOCUS_MAP = {
    "chrX": "AR (chrX)",
    "chr6": "ESR1 (chr6)",
    "chr14": "ESR2 (chr14)",
    "chr15": "CYP19A1 (chr15)",
}

ORDEN_MODALIDADES = ["ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF", "CAGE", "RNA_SEQ"]


def clasificar_tejido(biosample_name: str) -> str:
    """
    Clasifica cada biosample_name en compartimentos biológicos relevantes
    para separar el efecto mamario/adiposo del resto de tejidos endocrinos.
    """
    n = str(biosample_name).lower()
    if any(k in n for k in ["preadipocyte of the breast", "adipocyte of breast"]):
        return "Mammary_Adipose"
    if any(k in n for k in ["breast", "mammary", "mcf-7", "t47d"]):
        return "Mammary_Epithelial_Stromal"
    if any(k in n for k in ["adipose", "adipocyte", "preadipocyte", "fat pad"]):
        return "Adipose_General"
    return "Other_Reproductive_Endocrine"


def calcular_estadisticas_comparacion(df_sub: pd.DataFrame, etiqueta_conjunto: str) -> pd.DataFrame:
    """
    Calcula conteos de victorias AR_loss vs ESR1_gain, medianas de ratio,
    prueba binomial exacta (dos colas, p0=0.5) y prueba pareada de Wilcoxon.
    """
    filas = []
    for mod in ORDEN_MODALIDADES:
        grp = df_sub[df_sub["modalidad"] == mod]
        if grp.empty:
            continue
        n_total = len(grp)
        ar_wins = int((grp["ganador_por_impacto_conjunto"] == "AR_loss").sum())
        esr1_wins = int((grp["ganador_por_impacto_conjunto"] == "ESR1_gain").sum())
        empates = n_total - ar_wins - esr1_wins

        n_efectivo = ar_wins + esr1_wins
        p_binom = binomtest(ar_wins, n_efectivo, 0.5, alternative="two-sided").pvalue if n_efectivo > 0 else np.nan

        dif = grp["impacto_conjunto_abs_AR_loss"] - grp["impacto_conjunto_abs_ESR1_gain"]
        try:
            if np.all(dif == 0) or len(dif) < 2:
                p_wilcoxon = np.nan
            else:
                p_wilcoxon = wilcoxon(dif, alternative="two-sided").pvalue
        except Exception:
            p_wilcoxon = np.nan

        filas.append({
            "subconjunto": etiqueta_conjunto,
            "modalidad": mod,
            "n_tracks": n_total,
            "AR_loss_wins": ar_wins,
            "ESR1_gain_wins": esr1_wins,
            "empates": empates,
            "pct_AR_loss_wins": 100.0 * ar_wins / n_total if n_total > 0 else np.nan,
            "mediana_impacto_abs_AR_loss": grp["impacto_conjunto_abs_AR_loss"].median(),
            "mediana_impacto_abs_ESR1_gain": grp["impacto_conjunto_abs_ESR1_gain"].median(),
            "mediana_ratio_AR_sobre_ESR1": grp["ratio_abs_AR_sobre_ESR1"].median(),
            "media_ratio_AR_sobre_ESR1": grp["ratio_abs_AR_sobre_ESR1"].mean(),
            "mediana_epistasis_abs_AR_loss": grp["interaccion_abs_AR_loss"].median(),
            "mediana_epistasis_abs_ESR1_gain": grp["interaccion_abs_ESR1_gain"].median(),
            "ratio_mediana_epistasis_AR_sobre_ESR1": (
                grp["interaccion_abs_AR_loss"].median() / grp["interaccion_abs_ESR1_gain"].median()
                if grp["interaccion_abs_ESR1_gain"].median() > 0 else np.nan
            ),
            "p_valor_binomial_exacto": p_binom,
            "p_valor_wilcoxon_pareado": p_wilcoxon,
        })
    return pd.DataFrame(filas)


def generar_figura_publicacion(
    df_comp: pd.DataFrame,
    mam_adip: pd.DataFrame,
    conj_ma: pd.DataFrame,
    top_ar: pd.DataFrame,
    ruta_salida: Path,
) -> None:
    """
    Genera la figura compuesta de 4 paneles (A-D) lista para el manuscrito y el README.
    """
    plt.rcParams["font.sans-serif"] = "DejaVu Sans"
    plt.rcParams["axes.edgecolor"] = "#333333"
    plt.rcParams["axes.linewidth"] = 0.8

    fig, axes = plt.subplots(2, 2, figsize=(15, 11), dpi=300)
    fig.patch.set_facecolor("white")

    # Panel A: Dominancia en Tracks Mamarios y Adiposos (n=178)
    mod_stats = mam_adip.groupby("modalidad").agg(
        ar_wins=("ganador_por_impacto_conjunto", lambda x: (x == "AR_loss").sum()),
        esr1_wins=("ganador_por_impacto_conjunto", lambda x: (x == "ESR1_gain").sum()),
        total=("track_id", "count"),
    ).loc[ORDEN_MODALIDADES]

    mod_stats["ar_pct"] = 100.0 * mod_stats["ar_wins"] / mod_stats["total"]
    mod_stats["esr1_pct"] = 100.0 * mod_stats["esr1_wins"] / mod_stats["total"]

    ax1 = axes[0, 0]
    y_pos = np.arange(len(mod_stats))
    ax1.barh(
        y_pos,
        mod_stats["ar_pct"],
        color="#1f77b4",
        edgecolor="black",
        height=0.6,
        label="AR_loss (Pérdida de freno)",
    )
    ax1.barh(
        y_pos,
        mod_stats["esr1_pct"],
        left=mod_stats["ar_pct"],
        color="#ff7f0e",
        edgecolor="black",
        height=0.6,
        label="ESR1_gain (Ganancia estrogénica)",
    )
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(mod_stats.index, fontweight="bold")
    ax1.set_xlabel("% de Tracks ganados (Tejido Mamario y Adiposo, n=178)", fontweight="bold")
    ax1.set_title("A. Dominancia en Tracks Mamarios y Adiposos (164 vs 14)", fontweight="bold", fontsize=12)
    ax1.set_xlim(0, 108)
    for i, (ar_w, es_w, pct) in enumerate(zip(mod_stats["ar_wins"], mod_stats["esr1_wins"], mod_stats["ar_pct"])):
        ax1.text(
            pct / 2.0,
            i,
            f"{ar_w}/{ar_w + es_w} ({pct:.1f}%)",
            va="center",
            ha="center",
            color="white",
            fontweight="bold",
            fontsize=9.5,
        )
    ax1.legend(loc="lower right", frameon=True)
    ax1.invert_yaxis()

    # Panel B: Ratio Mediano AR_loss / ESR1_gain por Grupo Tisular y Modalidad
    ax2 = axes[0, 1]
    tg_mod = df_comp.groupby(["modalidad", "tissue_group"])["ratio_abs_AR_sobre_ESR1"].median().unstack()
    tg_mod = tg_mod.loc[ORDEN_MODALIDADES]
    x = np.arange(len(tg_mod))
    width = 0.25
    ax2.bar(
        x - width,
        tg_mod["Mammary_Epithelial_Stromal"],
        width,
        label="Epitelio/Estroma Mamario (n=149)",
        color="#2b5c8f",
        edgecolor="black",
    )
    ax2.bar(
        x,
        tg_mod["Adipose_General"],
        width,
        label="Tejido Adiposo (n=25)",
        color="#388e3c",
        edgecolor="black",
    )
    ax2.bar(
        x + width,
        tg_mod["Other_Reproductive_Endocrine"],
        width,
        label="Otros Endocrinos/Reprod. (n=101)",
        color="#7b1fa2",
        edgecolor="black",
    )
    ax2.axhline(1.0, color="red", linestyle="--", linewidth=1.5, label="Equilibrio Neutro (Ratio = 1.0)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(tg_mod.index, fontweight="bold", rotation=15)
    ax2.set_ylabel("Mediana Ratio Impacto (AR_loss / ESR1_gain)", fontweight="bold")
    ax2.set_title("B. Magnitud del Efecto Conjunto por Tejido y Modalidad", fontweight="bold", fontsize=12)
    ax2.legend(fontsize=8.5, frameon=True)

    # Panel C: Epistasis (No-aditividad L1) por Locus en Tejido Mamario/Adiposo (ChIP-TF)
    ax3 = axes[1, 0]
    ep_locus = (
        conj_ma[conj_ma["modalidad"] == "CHIP_TF"]
        .groupby(["locus", "categoria"])["no_aditividad_interaccion_abs_L1"]
        .median()
        .unstack()
    )
    orden_loci = ["ESR1 (chr6)", "ESR2 (chr14)", "CYP19A1 (chr15)", "AR (chrX)"]
    ep_locus = ep_locus.loc[orden_loci]
    x_loc = np.arange(len(ep_locus))
    ax3.bar(
        x_loc - 0.18,
        ep_locus["AR_loss"],
        0.35,
        label="AR_loss (Epistasis del haplotipo)",
        color="#1f77b4",
        edgecolor="black",
    )
    ax3.bar(
        x_loc + 0.18,
        ep_locus["ESR1_gain"],
        0.35,
        label="ESR1_gain (Epistasis del haplotipo)",
        color="#ff7f0e",
        edgecolor="black",
    )
    ax3.set_xticks(x_loc)
    ax3.set_xticklabels(ep_locus.index, fontweight="bold")
    ax3.set_ylabel("Mediana No-Aditividad |Conjunto - Suma Ind.| (ChIP-TF)", fontweight="bold")
    ax3.set_title("C. Sinergia Epistática por Locus en Tejido Mamario/Adiposo", fontweight="bold", fontsize=12)
    for i, (v_ar, v_es) in enumerate(zip(ep_locus["AR_loss"], ep_locus["ESR1_gain"])):
        ratio_ep = v_ar / v_es if v_es > 0 else np.nan
        ax3.text(
            i - 0.18,
            v_ar + 80,
            f"{v_ar:.0f}\n({ratio_ep:.1f}x)",
            ha="center",
            va="bottom",
            fontsize=8.5,
            fontweight="bold",
            color="#1f77b4",
        )
        ax3.text(
            i + 0.18,
            v_es + 80,
            f"{v_es:.0f}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="#d95f02",
        )
    ax3.set_ylim(0, max(ep_locus["AR_loss"]) * 1.25)
    ax3.legend(frameon=True)

    # Panel D: Top 10 Variantes Individuales AR_loss en Mamario/Adiposo
    ax4 = axes[1, 1]
    top10_plot = top_ar.head(10).iloc[::-1]
    labels_v = [f"{r['variante']} ({r['locus'].split()[0]})" for _, r in top10_plot.iterrows()]
    ax4.barh(labels_v, top10_plot["mean_abs_L1"], color="#008080", edgecolor="black", height=0.65)
    ax4.set_xlabel("Impacto Medio Absoluto L1 por Track Mamario/Adiposo (n=178)", fontweight="bold")
    ax4.set_title("D. Top 10 Variantes 'Driver' de AR_loss (Coordenadas hg38)", fontweight="bold", fontsize=12)
    for i, v in enumerate(top10_plot["mean_abs_L1"]):
        ax4.text(v + 5, i, f"{v:.1f}", va="center", fontsize=9, fontweight="bold")
    ax4.set_xlim(0, max(top10_plot["mean_abs_L1"]) * 1.15)

    plt.tight_layout()
    plt.savefig(ruta_salida, dpi=300)
    plt.close(fig)


def main() -> None:
    print("=== INICIANDO AUDITORÍA POST-HOC Y FILTRADO TEJIDO-ESPECÍFICO ===")
    df_reg = pd.read_csv(ARCHIVO_REGIONES)
    df_comp = pd.read_csv(ARCHIVO_COMPARACION)
    df_conj = pd.read_csv(ARCHIVO_CONJUNTO)
    df_ind = pd.read_csv(ARCHIVO_INDIVIDUAL)

    print(f"Regiones cargadas: {len(df_reg)} | Filas conjuntas: {len(df_conj)} | Filas individuales: {len(df_ind)}")

    # 1. Clasificación tisular
    df_comp["tissue_group"] = df_comp["biosample_name"].apply(clasificar_tejido)
    df_conj["tissue_group"] = df_conj["biosample_name"].apply(clasificar_tejido)
    df_ind["tissue_group"] = df_ind["biosample_name"].apply(clasificar_tejido)

    df_conj["locus"] = df_conj["cromosoma"].map(LOCUS_MAP)
    df_ind["locus"] = df_ind["cromosoma"].map(LOCUS_MAP)

    grupos_mamario_adiposo = ["Mammary_Epithelial_Stromal", "Mammary_Adipose", "Adipose_General"]
    mam_adip = df_comp[df_comp["tissue_group"].isin(grupos_mamario_adiposo)].copy()
    solo_adiposo = df_comp[df_comp["tissue_group"].isin(["Mammary_Adipose", "Adipose_General"])].copy()

    # 2. Tabla 1: Estadísticas globales, Mamario+Adiposo y Solo Adiposo
    t1_global = calcular_estadisticas_comparacion(df_comp, "Global_Todos_Tejidos (n=279)")
    t1_ma = calcular_estadisticas_comparacion(mam_adip, "Mamario_y_Adiposo (n=178)")
    t1_adip = calcular_estadisticas_comparacion(solo_adiposo, "Solo_Adiposo_y_Adipocito_Mamario (n=29)")
    tabla_1 = pd.concat([t1_global, t1_ma, t1_adip], ignore_index=True)
    ruta_t1 = CARPETA / "tabla_1_resumen_global_y_mamario_por_modalidad.csv"
    tabla_1.to_csv(ruta_t1, index=False)

    # Prueba binomial global y mamario-adiposa
    ar_tot_glob = int((df_comp["ganador_por_impacto_conjunto"] == "AR_loss").sum())
    es_tot_glob = int((df_comp["ganador_por_impacto_conjunto"] == "ESR1_gain").sum())
    p_glob = binomtest(ar_tot_glob, ar_tot_glob + es_tot_glob, 0.5).pvalue

    ar_tot_ma = int((mam_adip["ganador_por_impacto_conjunto"] == "AR_loss").sum())
    es_tot_ma = int((mam_adip["ganador_por_impacto_conjunto"] == "ESR1_gain").sum())
    p_ma = binomtest(ar_tot_ma, ar_tot_ma + es_tot_ma, 0.5).pvalue

    ar_tot_ad = int((solo_adiposo["ganador_por_impacto_conjunto"] == "AR_loss").sum())
    es_tot_ad = int((solo_adiposo["ganador_por_impacto_conjunto"] == "ESR1_gain").sum())
    p_ad = binomtest(ar_tot_ad, ar_tot_ad + es_tot_ad, 0.5).pvalue

    print(f"\n[1] Dominancia Global (n={len(df_comp)}): AR_loss={ar_tot_glob} vs ESR1_gain={es_tot_glob} (p_binom = {p_glob:.3e})")
    print(f"[2] Dominancia Mamario+Adiposo (n={len(mam_adip)}): AR_loss={ar_tot_ma} vs ESR1_gain={es_tot_ma} (p_binom = {p_ma:.3e})")
    print(f"[3] Dominancia Compartimento Adiposo (n={len(solo_adiposo)}): AR_loss={ar_tot_ad} vs ESR1_gain={es_tot_ad} (p_binom = {p_ad:.3e})")

    # 3. Tabla 2: Epistasis / No-aditividad por Locus en tejido Mamario+Adiposo
    conj_ma = df_conj[df_conj["tissue_group"].isin(grupos_mamario_adiposo)].copy()
    tabla_2 = (
        conj_ma.groupby(["locus", "modalidad", "categoria"], as_index=False)
        .agg(
            n_observaciones=("track_id", "count"),
            mediana_efecto_conjunto_abs_L1=("efecto_conjunto_abs_L1", "median"),
            mediana_efecto_por_variante=("efecto_conjunto_abs_por_variante", "median"),
            mediana_no_aditividad_abs_L1=("no_aditividad_interaccion_abs_L1", "median"),
            mediana_fraccion_epistasis_sobre_suma_ind=("interaccion_abs_sobre_suma_individual_abs", "median"),
        )
    )
    ruta_t2 = CARPETA / "tabla_2_epistasis_no_aditividad_por_locus.csv"
    tabla_2.to_csv(ruta_t2, index=False)

    # 4. Tabla 3: Ranking completo de variantes individuales en tejido Mamario+Adiposo
    ind_ma = df_ind[df_ind["tissue_group"].isin(grupos_mamario_adiposo)].copy()
    var_summary = (
        ind_ma.groupby(["categoria", "locus", "cromosoma", "region_id", "variante"], as_index=False)
        .agg(
            n_tracks_mamario_adiposo=("track_id", "count"),
            sum_abs_L1=("efecto_individual_abs_L1", "sum"),
            mean_abs_L1=("efecto_individual_abs_L1", "mean"),
            median_abs_L1=("efecto_individual_abs_L1", "median"),
            max_abs_L1=("efecto_individual_abs_L1", "max"),
            mean_signed_humano_menos_ancestral=("efecto_individual_firmado_humano_menos_ancestral", "mean"),
        )
        .sort_values(["categoria", "sum_abs_L1"], ascending=[True, False])
    )
    var_summary["rank_dentro_de_categoria"] = var_summary.groupby("categoria")["sum_abs_L1"].rank(
        ascending=False, method="first"
    ).astype(int)
    var_summary["rank_dentro_de_locus_y_categoria"] = var_summary.groupby(["categoria", "locus"])["sum_abs_L1"].rank(
        ascending=False, method="first"
    ).astype(int)

    ruta_t3 = CARPETA / "tabla_3_top_variantes_drivers_mamario_adiposo.csv"
    var_summary.to_csv(ruta_t3, index=False)

    # 5. Tabla 4: Cross-Talk en tracks con factor directo (ESR1 ChIP-seq en MCF-7)
    tabla_4 = df_conj[
        df_conj["factor_directo_detectado"].notna() & (df_conj["factor_directo_detectado"] != "")
    ][
        [
            "region_id",
            "cromosoma",
            "inicio_region_hg38_1based",
            "fin_region_hg38_1based",
            "categoria",
            "modalidad",
            "biosample_name",
            "factor_directo_detectado",
            "n_variantes_conjuntas",
            "efecto_conjunto_firmado_humano_menos_ancestral",
            "efecto_conjunto_abs_L1",
            "no_aditividad_interaccion_firmada",
            "no_aditividad_interaccion_abs_L1",
            "direccion_concordante_en_track_factor_directo",
        ]
    ].copy()
    ruta_t4 = CARPETA / "tabla_4_cross_talk_factores_directos_MCF7.csv"
    tabla_4.to_csv(ruta_t4, index=False)

    # 6. Generación de figura compuesta de 4 paneles
    top_ar = var_summary[var_summary["categoria"] == "AR_loss"].sort_values("sum_abs_L1", ascending=False)
    ruta_fig = CARPETA / "resumen_validacion_alphagenome_mamario.png"
    generar_figura_publicacion(df_comp, mam_adip, conj_ma, top_ar, ruta_fig)

    print("\n=== ARCHIVOS GENERADOS PARA REPOSITORIO Y MANUSCRITO ===")
    print(f"1. {ruta_t1.name}")
    print(f"2. {ruta_t2.name}")
    print(f"3. {ruta_t3.name}")
    print(f"4. {ruta_t4.name}")
    print(f"5. {ruta_fig.name}")


if __name__ == "__main__":
    main()
