#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_multiclass_report.py
==========================================
Combina de forma 100% dinámica y reproducible los resultados de Sustituciones/Deleciones
(SNVs/Indels cortos) e Inserciones Específicas Humanas (HSIs) a nivel de track (n=178
tracks mamarios/adiposos y n=279 globales), evaluando las modalidades Conjunto, Aditivo,
No Aditivo (Epistasis) y el Gran Total Global Multiclase.
"""

from pathlib import Path
import numpy as np
import pandas as pd

CARPETA = Path(__file__).resolve().parent


def clasificar_tejido(name: str) -> str:
    n = str(name).lower()
    if any(k in n for k in ['preadipocyte of the breast', 'adipocyte of breast']):
        return 'Mammary_Adipose'
    if any(k in n for k in ['breast', 'mammary', 'mcf-7', 't47d']):
        return 'Mammary_Epithelial_Stromal'
    if any(k in n for k in ['adipose', 'adipocyte', 'preadipocyte', 'fat pad']):
        return 'Adipose_General'
    return 'Other_Reproductive_Endocrine'


def construir_fila(nombre: str, serie_ar: pd.Series, serie_esr1: pd.Series) -> dict:
    tot = len(serie_ar)
    ar_wins = int((serie_ar > serie_esr1).sum())
    esr1_wins = int((serie_esr1 > serie_ar).sum())
    pct_ar = round(100.0 * ar_wins / tot, 2) if tot > 0 else 0.0
    mediana_ratio = round(float((serie_ar / serie_esr1).median()), 3) if tot > 0 else np.nan
    suma_ar = round(float(serie_ar.sum()), 2)
    suma_esr1 = round(float(serie_esr1.sum()), 2)
    return {
        "Clase_y_Modalidad": nombre,
        "n_tracks_evaluados": tot,
        "AR_loss_victorias": ar_wins,
        "ESR1_gain_victorias": esr1_wins,
        "pct_dominancia_AR_loss": pct_ar,
        "mediana_ratio_AR_sobre_ESR1": mediana_ratio,
        "suma_total_L1_AR_loss": suma_ar,
        "suma_total_L1_ESR1_gain": suma_esr1,
        "Ganador_Global_Peso": "AR_loss" if suma_ar > suma_esr1 else "ESR1_gain"
    }


def main():
    print("=== CONSOLIDANDO ANÁLISIS MULTICLASE (SNVs/INDELS + INSERCIONES HSIs) ===")

    ruta_snv = CARPETA / "resultados_impacto_conjunto_comparacion.csv"
    ruta_ins = CARPETA / "tabla_5_impacto_inserciones_mamario.csv"

    if not ruta_snv.exists():
        raise FileNotFoundError(f"No se encontró el archivo de SNVs: {ruta_snv}")
    if not ruta_ins.exists():
        raise FileNotFoundError(f"No se encontró el archivo de Inserciones: {ruta_ins}")

    # 1. Cargar SNVs y clasificar tejido
    snv_comp = pd.read_csv(ruta_snv)
    snv_comp['tissue_group'] = snv_comp['biosample_name'].apply(clasificar_tejido)

    # Identificador interno de hebra (+ / -) para tracks con el mismo nombre (ej. CAGE, RNA_SEQ)
    cols_meta = ['modalidad', 'track_name', 'ontology_curie', 'biosample_name']
    snv_comp['strand_idx'] = snv_comp.groupby(cols_meta).cumcount()

    # 2. Cargar Inserciones y agregar las 12 regiones por cada uno de los 279 tracks únicos
    ins_df = pd.read_csv(ruta_ins)
    if 'tissue_group' not in ins_df.columns:
        ins_df['tissue_group'] = ins_df['biosample_name'].apply(clasificar_tejido)

    ins_df['strand_idx'] = ins_df.groupby(['region_id', 'categoria'] + cols_meta).cumcount()

    # Suma exacta de las 12 regiones por track y categoría (misma unidad que SNVs)
    ins_por_track = ins_df.pivot_table(
        index=cols_meta + ['tissue_group', 'strand_idx'],
        columns='categoria',
        values=['efecto_insercion_abs_L1', 'bp_totales_insertados'],
        aggfunc='sum'
    ).reset_index()

    ins_flat = pd.DataFrame({
        'modalidad': ins_por_track[('modalidad', '')],
        'track_name': ins_por_track[('track_name', '')],
        'ontology_curie': ins_por_track[('ontology_curie', '')],
        'biosample_name': ins_por_track[('biosample_name', '')],
        'tissue_group': ins_por_track[('tissue_group', '')],
        'strand_idx': ins_por_track[('strand_idx', '')],
        'ins_abs_AR_loss': ins_por_track[('efecto_insercion_abs_L1', 'AR_loss')],
        'ins_abs_ESR1_gain': ins_por_track[('efecto_insercion_abs_L1', 'ESR1_gain')],
        'bp_ins_AR_loss': ins_por_track[('bp_totales_insertados', 'AR_loss')],
        'bp_ins_ESR1_gain': ins_por_track[('bp_totales_insertados', 'ESR1_gain')],
    })

    # 3. Emparejar 1 a 1 cada track entre SNVs e Inserciones
    merged = pd.merge(
        snv_comp,
        ins_flat,
        on=cols_meta + ['tissue_group', 'strand_idx'],
        how='inner'
    )

    # Filtrar los 178 tracks de Tejido Mamario y Adiposo
    grupos_mamarios = ['Mammary_Epithelial_Stromal', 'Mammary_Adipose', 'Adipose_General']
    mam = merged[merged['tissue_group'].isin(grupos_mamarios)].copy()

    # 4. Calcular las 3 modalidades en SNVs, métricas de Inserciones y el Gran Total Global
    # En el pipeline: interaccion_abs = |Efecto_Conjunto - Suma_Aditiva_Individual|
    mam['snv_aditivo_AR'] = mam['interaccion_abs_AR_loss'] - mam['impacto_conjunto_abs_AR_loss']
    mam['snv_aditivo_ESR1'] = mam['interaccion_abs_ESR1_gain'] - mam['impacto_conjunto_abs_ESR1_gain']

    # Impacto de inserciones normalizado por par de base (bp) insertado
    mam['ins_per_bp_AR'] = mam['ins_abs_AR_loss'] / mam['bp_ins_AR_loss']
    mam['ins_per_bp_ESR1'] = mam['ins_abs_ESR1_gain'] / mam['bp_ins_ESR1_gain']

    # Impacto Global Multiclase por track (SNVs + Inserciones)
    mam['global_conjunto_AR'] = mam['impacto_conjunto_abs_AR_loss'] + mam['ins_abs_AR_loss']
    mam['global_conjunto_ESR1'] = mam['impacto_conjunto_abs_ESR1_gain'] + mam['ins_abs_ESR1_gain']

    mam['global_aditivo_AR'] = mam['snv_aditivo_AR'] + mam['ins_abs_AR_loss']
    mam['global_aditivo_ESR1'] = mam['snv_aditivo_ESR1'] + mam['ins_abs_ESR1_gain']

    # 5. Construir tabla consolidada
    consolidado_global = pd.DataFrame([
        construir_fila(
            "1. SNVs/Indels — Modelo Conjunto",
            mam['impacto_conjunto_abs_AR_loss'],
            mam['impacto_conjunto_abs_ESR1_gain']
        ),
        construir_fila(
            "2. SNVs/Indels — Modelo Aditivo (Suma Individual)",
            mam['snv_aditivo_AR'],
            mam['snv_aditivo_ESR1']
        ),
        construir_fila(
            "3. SNVs/Indels — Modelo No Aditivo (Epistasis)",
            mam['interaccion_abs_AR_loss'],
            mam['interaccion_abs_ESR1_gain']
        ),
        construir_fila(
            "4. Inserciones (HSIs) — Densidad por bp insertado",
            mam['ins_per_bp_AR'],
            mam['ins_per_bp_ESR1']
        ),
        construir_fila(
            "5. Inserciones (HSIs) — Volumen Bruto (12 regiones)",
            mam['ins_abs_AR_loss'],
            mam['ins_abs_ESR1_gain']
        ),
        construir_fila(
            "6. GLOBAL MULTICLASE — Conjunto (SNVs + HSIs)",
            mam['global_conjunto_AR'],
            mam['global_conjunto_ESR1']
        ),
        construir_fila(
            "7. GLOBAL MULTICLASE — Aditivo Total (SNVs + HSIs)",
            mam['global_aditivo_AR'],
            mam['global_aditivo_ESR1']
        ),
    ])

    ruta_cons = CARPETA / "tabla_6_consolidado_multiclase_evolutivo.csv"
    consolidado_global.to_csv(ruta_cons, index=False)
    print(f"\n[ÉXITO] Tabla consolidada multiclase guardada en: {ruta_cons.name}\n")
    print(consolidado_global.to_string(index=False))


if __name__ == "__main__":
    main()