#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
auditoria_inserciones.py
========================
Módulo complementario (ortogonal a automatizar_bloques.py) para detectar, clasificar
y evaluar en AlphaGenome las Inserciones Específicas del Linaje Humano (HSIs / Indels)
en los 4 loci candidatos (AR, ESR1, ESR2, CYP19A1) sin repetir la corrida de SNVs/deleciones.

Versión Auto-Resume Inteligente:
  - Lee tabla_5_impacto_inserciones_mamario.csv y procesa ÚNICAMENTE las regiones faltantes (ej. Región 6 si falló).
"""

import argparse
import gzip
import json
import os
import re
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

try:
    import genoma_kernels  # type: ignore
    TIENE_CPP_KERNEL = True
except ImportError:
    genoma_kernels = None
    TIENE_CPP_KERNEL = False

CARPETA = Path(__file__).resolve().parent

ARCHIVO_REGIONES = CARPETA / "resultados_impacto_conjunto_regiones.csv"
ARCHIVO_COMPARACION_SNV = CARPETA / "resultados_impacto_conjunto_comparacion.csv"
ARCHIVO_SALIDA_INSERCIONES_CSV = CARPETA / "variantes_inserciones_loci.csv"
ARCHIVO_SALIDA_TABLA5 = CARPETA / "tabla_5_impacto_inserciones_mamario.csv"
ARCHIVO_SALIDA_RESUMEN_INS = CARPETA / "tabla_5_resumen_inserciones_por_modalidad.csv"

CONTEXT_LENGTH = 1048576  # 1 Mb exacto requerido por AlphaGenome

LOCUS_MAP = {
    "chrX": "AR (chrX)",
    "chr6": "ESR1 (chr6)",
    "chr14": "ESR2 (chr14)",
    "chr15": "CYP19A1 (chr15)",
}

ONTOLOGY_TERMS_DEFAULT = [
    "CL:0002231", "CL:0002324", "CL:0002326", "CL:0002327", "CL:0002451",
    "CL:0002555", "CL:0002570", "CL:0002577", "CL:0002580", "CL:0002601",
    "CL:0002617", "CL:0002622", "CL:1000487", "CL:2000071", "CL:2000078",
    "EFO:0001203", "EFO:0001247", "NTR:0003013", "NTR:0003079",
    "UBERON:0000002", "UBERON:0000007", "UBERON:0000310", "UBERON:0000473",
    "UBERON:0000992", "UBERON:0000995", "UBERON:0001013", "UBERON:0001987",
    "UBERON:0002118", "UBERON:0002119", "UBERON:0002190", "UBERON:0002367",
    "UBERON:0005795", "UBERON:0008367", "UBERON:0010414", "UBERON:0014455",
]

PWM_ARE_MATRIX = np.array([
    [ 1.2, -1.5,  1.4,  1.3, -1.8,  1.5,  0.0,  0.0,  0.0, -1.5, -1.8, -1.5, -1.8, -1.5, -1.5],
    [-1.5, -1.8, -1.5, -1.5,  1.8, -1.5,  0.0,  0.0,  0.0, -1.2, -1.8, -1.5, -1.5,  1.7, -1.2],
    [ 0.8,  1.8, -1.5, -1.5, -1.8, -1.5,  0.0,  0.0,  0.0, -1.8,  1.9, -1.8, -1.8, -1.5, -1.8],
    [-1.5, -1.5,  0.5,  0.4, -1.5, -1.5,  0.0,  0.0,  0.0,  1.7, -1.8,  1.8,  1.8,  0.5,  1.7],
], dtype=np.float64)

PWM_ERE_MATRIX = np.array([
    [ 1.7, -1.8, -1.8, -1.5, -1.8,  1.8,  0.0,  0.0,  0.0, -1.8, -1.8,  1.8, -1.8, -1.8, -1.5],
    [-1.8, -1.8, -1.8, -1.5,  1.9, -1.8,  0.0,  0.0,  0.0, -1.5, -1.8, -1.8,  1.9,  1.9, -1.5],
    [-1.5,  1.9,  1.9, -1.8, -1.8, -1.8,  0.0,  0.0,  0.0, -1.8,  1.9, -1.8, -1.8, -1.8, -1.8],
    [-1.8, -1.8, -1.8,  1.8, -1.8, -1.8,  0.0,  0.0,  0.0,  1.8, -1.8, -1.8, -1.8, -1.8,  1.8],
], dtype=np.float64)

PATRONES_ARE = [
    re.compile(r"[AG]G[AT]ACA[ACGT]{3}TGTT[CT][T]", re.IGNORECASE),
    re.compile(r"AGAACA|TGTTCT|GGTACA|TGTACC", re.IGNORECASE),
]
PATRONES_ERE = [
    re.compile(r"AGGTCA[ACGT]{3}TGACCT", re.IGNORECASE),
    re.compile(r"AGGTCA|TGACCT|GGTCAN{3}TGACC", re.IGNORECASE),
]

LOOKUP_ASCII_BASE = np.full(256, -1, dtype=np.int8)
LOOKUP_ASCII_BASE[ord("A")] = 0
LOOKUP_ASCII_BASE[ord("a")] = 0
LOOKUP_ASCII_BASE[ord("C")] = 1
LOOKUP_ASCII_BASE[ord("c")] = 1
LOOKUP_ASCII_BASE[ord("G")] = 2
LOOKUP_ASCII_BASE[ord("g")] = 2
LOOKUP_ASCII_BASE[ord("T")] = 3
LOOKUP_ASCII_BASE[ord("t")] = 3


def clasificar_tejido(biosample_name: str) -> str:
    n = str(biosample_name).lower()
    if any(k in n for k in ["preadipocyte of the breast", "adipocyte of breast"]):
        return "Mammary_Adipose"
    if any(k in n for k in ["breast", "mammary", "mcf-7", "t47d"]):
        return "Mammary_Epithelial_Stromal"
    if any(k in n for k in ["adipose", "adipocyte", "preadipocyte", "fat pad"]):
        return "Adipose_General"
    return "Other_Reproductive_Endocrine"


def evaluar_afinidad_cpp_o_simd(seq: str, pwm: np.ndarray, patrones_regex: List[re.Pattern]) -> float:
    s = seq.upper().replace("-", "")
    if not s:
        return 0.0

    if TIENE_CPP_KERNEL and hasattr(genoma_kernels, "score_pwm_sequence"):
        score_cpp = float(genoma_kernels.score_pwm_sequence(s, pwm))
    else:
        k_len = pwm.shape[1]
        if len(s) >= k_len:
            buf = np.frombuffer(s.encode("ascii", errors="ignore"), dtype=np.uint8)
            idx = LOOKUP_ASCII_BASE[buf]
            ventanas = np.lib.stride_tricks.sliding_window_view(idx, k_len)
            mascara_valida = np.all(ventanas >= 0, axis=1)
            if np.any(mascara_valida):
                v_validas = ventanas[mascara_valida]
                pos_cols = np.arange(k_len)
                scores_ventanas = np.sum(pwm[v_validas, pos_cols], axis=1)
                score_cpp = float(np.max(scores_ventanas)) if scores_ventanas.size > 0 else 0.0
            else:
                score_cpp = 0.0
        else:
            score_cpp = 0.0

    score_regex = 3.0 * len(patrones_regex[0].findall(s)) + 1.0 * len(patrones_regex[1].findall(s))
    return max(0.0, score_cpp) + score_regex


def abrir_maf(ruta: Path):
    if str(ruta).endswith(".gz"):
        return gzip.open(ruta, "rt", encoding="utf-8", errors="replace")
    return open(ruta, "r", encoding="utf-8", errors="replace")


def _procesar_bloque_maf_vectorizado(lineas_s, regiones_por_chr, min_len_bp, max_len_bp, salida):
    registros = {}
    for l in lineas_s:
        partes = l.split()
        if len(partes) < 7:
            continue
        src, start, size, strand, srcSize, text = partes[1], int(partes[2]), int(partes[3]), partes[4], int(partes[5]), partes[6]
        especie = src.split(".")[0]
        registros[especie] = (src, start, size, strand, srcSize, text)

    if "hg38" not in registros:
        return

    src_hg, start_hg0, size_hg, _, _, seq_hg = registros["hg38"]
    chrom = src_hg.split(".", 1)[1] if "." in src_hg else src_hg
    if chrom not in regiones_por_chr:
        return

    end_hg0 = start_hg0 + size_hg
    regiones_activas = [(rid, ini, fin) for (rid, ini, fin) in regiones_por_chr[chrom] if not (end_hg0 < ini or start_hg0 > fin)]
    if not regiones_activas:
        return

    primates_nhp = [esp for esp in ["panTro6", "panTro5", "gorGor6", "ponAbe3", "rheMac10"] if esp in registros]
    if not primates_nhp:
        return

    GAP_BYTE = 45
    buf_hg = np.frombuffer(seq_hg.encode("ascii"), dtype=np.uint8)
    n_cols = buf_hg.size
    hg_es_base = buf_hg != GAP_BYTE
    nhp_todos_gap = np.ones(n_cols, dtype=bool)
    for esp in primates_nhp:
        buf_nhp = np.frombuffer(registros[esp][5].encode("ascii"), dtype=np.uint8)
        if buf_nhp.size == n_cols:
            nhp_todos_gap &= (buf_nhp == GAP_BYTE)

    es_insercion_humana = hg_es_base & nhp_todos_gap
    if not np.any(es_insercion_humana):
        return

    coords_hg_1based = start_hg0 + np.cumsum(hg_es_base, dtype=np.int64)
    padded = np.empty(n_cols + 2, dtype=np.int8)
    padded[0] = 0
    padded[1:-1] = es_insercion_humana.view(np.int8)
    padded[-1] = 0
    cambios = np.diff(padded)
    inicios_run = np.flatnonzero(cambios == 1)
    finales_run = np.flatnonzero(cambios == -1)

    for i_ini, j_fin in zip(inicios_run, finales_run):
        largo_ins = int(j_fin - i_ini)
        if not (min_len_bp <= largo_ins <= max_len_bp):
            continue
        pos_inicio_ins = int(coords_hg_1based[i_ini])
        izq = max(0, int(i_ini) - 15)
        der = min(n_cols, int(j_fin) + 15)
        ventana_humana = seq_hg[izq:der].replace("-", "").upper()
        ventana_ancestral = (seq_hg[izq:int(i_ini)] + seq_hg[int(j_fin):der]).replace("-", "").upper()

        are_hum = evaluar_afinidad_cpp_o_simd(ventana_humana, PWM_ARE_MATRIX, PATRONES_ARE)
        are_anc = evaluar_afinidad_cpp_o_simd(ventana_ancestral, PWM_ARE_MATRIX, PATRONES_ARE)
        ere_hum = evaluar_afinidad_cpp_o_simd(ventana_humana, PWM_ERE_MATRIX, PATRONES_ERE)
        ere_anc = evaluar_afinidad_cpp_o_simd(ventana_ancestral, PWM_ERE_MATRIX, PATRONES_ERE)
        seq_ins_str = seq_hg[int(i_ini):int(j_fin)].replace("-", "").upper()

        for rid, ini, fin in regiones_activas:
            if ini <= pos_inicio_ins <= fin:
                if are_anc > are_hum + 0.5:
                    salida.append({"region_id": rid, "cromosoma": chrom, "pos_hg38_1based": pos_inicio_ins, "longitud_ins_bp": largo_ins, "secuencia_insertada_humana": seq_ins_str, "categoria": "AR_loss", "delta_motivo": round(are_hum - are_anc, 4), "variante": f"{chrom}:{pos_inicio_ins}:INS_{seq_ins_str}"})
                if ere_hum > ere_anc + 0.5:
                    salida.append({"region_id": rid, "cromosoma": chrom, "pos_hg38_1based": pos_inicio_ins, "longitud_ins_bp": largo_ins, "secuencia_insertada_humana": seq_ins_str, "categoria": "ESR1_gain", "delta_motivo": round(ere_hum - ere_anc, 4), "variante": f"{chrom}:{pos_inicio_ins}:INS_{seq_ins_str}"})


def escanear_maf_en_busca_de_inserciones(df_reg, min_len_bp=1, max_len_bp=50):
    archivos_maf = sorted(list(CARPETA.glob("*.maf")) + list(CARPETA.glob("*.maf.gz")))
    if not archivos_maf:
        return pd.DataFrame()
    regiones_por_chr = {}
    for _, r in df_reg.iterrows():
        c = str(r["cromosoma"])
        regiones_por_chr.setdefault(c, []).append((int(r["region_id"]), int(r["inicio_region_hg38_1based"]), int(r["fin_region_hg38_1based"])))
    inserciones_detectadas = []
    for maf_path in archivos_maf:
        with abrir_maf(maf_path) as f:
            lineas_s = []
            for linea in f:
                if linea.startswith("s "):
                    lineas_s.append(linea.strip())
                elif linea.strip() == "" and lineas_s:
                    _procesar_bloque_maf_vectorizado(lineas_s, regiones_por_chr, min_len_bp, max_len_bp, inserciones_detectadas)
                    lineas_s = []
            if lineas_s:
                _procesar_bloque_maf_vectorizado(lineas_s, regiones_por_chr, min_len_bp, max_len_bp, inserciones_detectadas)
    df_ins = pd.DataFrame(inserciones_detectadas)
    if not df_ins.empty:
        df_ins = df_ins.drop_duplicates(subset=["region_id", "cromosoma", "pos_hg38_1based", "secuencia_insertada_humana", "categoria"]).sort_values(["region_id", "pos_hg38_1based"])
    return df_ins


def descargar_secuencia_ucsc(chrom, start_1based, end_1based, timeout=90):
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start_1based - 1};end={end_1based}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    ultimo_error = None
    for intento in range(1, 4):
        try:
            try:
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    data = json.loads(resp.read().decode("utf-8"))
            except urllib.error.URLError as e:
                if "CERTIFICATE_VERIFY_FAILED" in str(e) or isinstance(getattr(e, "reason", None), ssl.SSLError):
                    ctx_fallback = ssl._create_unverified_context()
                    with urllib.request.urlopen(req, timeout=timeout, context=ctx_fallback) as resp:
                        data = json.loads(resp.read().decode("utf-8"))
                else:
                    raise
            seq = data.get("dna", "").upper()
            if len(seq) != (end_1based - start_1based + 1):
                raise ValueError("Longitud inesperada")
            return seq
        except Exception as exc:
            ultimo_error = exc
            if intento < 3:
                time.sleep(3 * intento)
    raise RuntimeError(f"Fallo UCSC: {ultimo_error}")


def aplicar_reversion_inserciones_in_place(seq_humana_1mb, win_start_1based, df_ins_sub):
    buf = bytearray(seq_humana_1mb.encode("ascii"))
    n_total = len(buf)
    N_BYTE = ord("N")
    for pos1, largo in zip(df_ins_sub["pos_hg38_1based"].to_numpy(), df_ins_sub["longitud_ins_bp"].to_numpy()):
        idx0 = int(pos1) - win_start_1based
        if 0 <= idx0 < n_total:
            fin0 = min(n_total, idx0 + int(largo))
            buf[idx0:fin0] = bytes([N_BYTE]) * (fin0 - idx0)
    return buf.decode("ascii")


def cargar_api_key(ruta_archivo_key=""):
    candidatos = [Path(ruta_archivo_key) if ruta_archivo_key else None]
    for nombre in ["api_key.txt", "alphagenome_api_key.txt", "alphagenome_key.txt", "apikey.txt", "key.txt", ".env"]:
        candidatos.append(CARPETA / nombre)
    for extra in CARPETA.glob("*key*.txt"):
        candidatos.append(extra)
    for ruta in candidatos:
        if ruta and ruta.exists() and ruta.is_file():
            for linea in ruta.read_text(encoding="utf-8", errors="ignore").splitlines():
                l = linea.strip()
                if not l or l.startswith("#"):
                    continue
                val = l.split("=", 1)[1].strip().strip('"').strip("'") if "=" in l and "KEY" in l.upper() else l.strip('"').strip("'")
                if val:
                    print(f"[API KEY] Cargada desde: {ruta.name}")
                    return val
    env_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if env_key:
        return env_key
    raise RuntimeError("No se encontró la API key.")


def obtener_parametros_prediccion(dna_client_mod):
    if ARCHIVO_COMPARACION_SNV.exists():
        df_c = pd.read_csv(ARCHIVO_COMPARACION_SNV)
        ont_terms = sorted(df_c["ontology_curie"].dropna().astype(str).unique().tolist())
        mod_names = sorted(df_c["modalidad"].dropna().astype(str).unique().tolist())
    else:
        ont_terms = ONTOLOGY_TERMS_DEFAULT
        mod_names = ["ATAC", "CAGE", "CHIP_HISTONE", "CHIP_TF", "DNASE", "RNA_SEQ"]
    output_type_enum = getattr(dna_client_mod, "OutputType", None)
    req_outputs = [getattr(output_type_enum, m) for m in mod_names if hasattr(output_type_enum, m)] if output_type_enum else mod_names
    return req_outputs, ont_terms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ejecutar-alphagenome", action="store_true")
    parser.add_argument("--min-len", type=int, default=1)
    parser.add_argument("--max-len", type=int, default=50)
    parser.add_argument("--api-key-file", type=str, default="")
    args = parser.parse_args()

    df_reg = pd.read_csv(ARCHIVO_REGIONES)
    if ARCHIVO_SALIDA_INSERCIONES_CSV.exists():
        df_ins = pd.read_csv(ARCHIVO_SALIDA_INSERCIONES_CSV)
    else:
        df_ins = escanear_maf_en_busca_de_inserciones(df_reg, min_len_bp=args.min_len, max_len_bp=args.max_len)
        if not df_ins.empty:
            df_ins.to_csv(ARCHIVO_SALIDA_INSERCIONES_CSV, index=False)

    if df_ins.empty:
        return

    if not args.ejecutar_alphagenome:
        return

    api_key = cargar_api_key(args.api_key_file)
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    client = dna_client.create(api_key)
    requested_outputs, ontology_terms = obtener_parametros_prediccion(dna_client)

    filas_salida = []
    regiones_completadas = set()
    if ARCHIVO_SALIDA_TABLA5.exists():
        try:
            df_existente = pd.read_csv(ARCHIVO_SALIDA_TABLA5)
            if not df_existente.empty and "region_id" in df_existente.columns:
                filas_salida = df_existente.to_dict("records")
                regiones_completadas = set(df_existente["region_id"].unique().tolist())
                print(f"[AUTO-RESUME] Regiones ya completadas en tabla_5: {sorted(list(regiones_completadas))}")
        except Exception:
            pass

    regiones_pendientes = [r for _, r in df_reg.iterrows() if int(r["region_id"]) not in regiones_completadas]
    if not regiones_pendientes:
        print("[AUTO-RESUME] ¡Todas las regiones ya están procesadas!")
        return

    print(f"[AUTO-RESUME] Regiones pendientes a procesar: {[int(r['region_id']) for r in regiones_pendientes]}")
    print("\n=== INICIANDO FASE 2 (SOLO PENDIENTES) EN ALPHAGENOME ===")

    for reg in regiones_pendientes:
        rid = int(reg["region_id"])
        chrom, ini, fin = str(reg["cromosoma"]), int(reg["inicio_region_hg38_1based"]), int(reg["fin_region_hg38_1based"])
        centro = (ini + fin) // 2
        win_start, win_end = max(1, centro - CONTEXT_LENGTH // 2), max(1, centro - CONTEXT_LENGTH // 2) + CONTEXT_LENGTH - 1
        sub_ins = df_ins[df_ins["region_id"] == rid]
        if sub_ins.empty:
            continue

        exito = False
        for intento in range(1, 4):
            try:
                print(f"[REGIÓN {rid}] Descargando 1Mb y evaluando inserciones (intento {intento}/3)...")
                seq_humana = descargar_secuencia_ucsc(chrom, win_start, win_end)
                intervalo = genome.Interval(chromosome=chrom, start=win_start - 1, end=win_end)
                out_hum = client.predict_sequence(sequence=seq_humana, interval=intervalo, requested_outputs=requested_outputs, ontology_terms=ontology_terms)

                for cat in ["AR_loss", "ESR1_gain"]:
                    df_cat = sub_ins[sub_ins["categoria"] == cat]
                    if df_cat.empty:
                        continue
                    seq_anc_sin_ins = aplicar_reversion_inserciones_in_place(seq_humana, win_start, df_cat)
                    out_anc = client.predict_sequence(sequence=seq_anc_sin_ins, interval=intervalo, requested_outputs=requested_outputs, ontology_terms=ontology_terms)

                    for mod_name in ["atac", "dnase", "chip_histone", "chip_tf", "cage", "rna_seq"]:
                        track_h, track_a = getattr(out_hum, mod_name, None), getattr(out_anc, mod_name, None)
                        if track_h is None or track_a is None:
                            continue
                        vals_h, vals_a, meta = np.asarray(track_h.values, dtype=np.float64), np.asarray(track_a.values, dtype=np.float64), track_h.metadata
                        diff_l1 = np.sum(np.abs(vals_h - vals_a), axis=0)
                        diff_signed = np.sum(vals_h - vals_a, axis=0)

                        for idx_t in range(len(diff_l1)):
                            bname = str(meta.iloc[idx_t].get("biosample_name", "")) if idx_t < len(meta) else ""
                            tname = str(meta.iloc[idx_t].get("name", f"{mod_name}_{idx_t}")) if idx_t < len(meta) else ""
                            ocurie = str(meta.iloc[idx_t].get("ontology_curie", "")) if idx_t < len(meta) else ""
                            filas_salida.append({
                                "region_id": rid, "cromosoma": chrom, "locus": LOCUS_MAP.get(chrom, chrom),
                                "categoria": cat, "modalidad": mod_name.upper(), "track_name": tname,
                                "ontology_curie": ocurie, "biosample_name": bname, "tissue_group": clasificar_tejido(bname),
                                "n_inserciones_en_grupo": len(df_cat), "bp_totales_insertados": int(df_cat["longitud_ins_bp"].sum()),
                                "efecto_insercion_abs_L1": float(diff_l1[idx_t]), "efecto_insercion_firmado": float(diff_signed[idx_t])
                            })

                exito = True
                print(f"  [OK] Región {rid} evaluada exitosamente.")
                break
            except Exception as e:
                print(f"  [AVISO] Región {rid} intento {intento} falló: {e}")
                time.sleep(5 * intento)

        if exito:
            df_res_ins = pd.DataFrame(filas_salida)
            df_res_ins.to_csv(ARCHIVO_SALIDA_TABLA5, index=False)
            print(f"  [CHECKPOINT] {ARCHIVO_SALIDA_TABLA5.name} actualizado con Región {rid}.")

    if filas_salida:
        df_res_ins = pd.DataFrame(filas_salida)
        df_res_ins.to_csv(ARCHIVO_SALIDA_TABLA5, index=False)

        piv = (
            df_res_ins.groupby(["tissue_group", "modalidad", "track_name", "categoria"])["efecto_insercion_abs_L1"]
            .sum()
            .unstack("categoria")
            .dropna()
        )
        if "AR_loss" in piv.columns and "ESR1_gain" in piv.columns:
            piv["ratio_AR_sobre_ESR1"] = piv["AR_loss"] / np.where(piv["ESR1_gain"] > 0, piv["ESR1_gain"], np.nan)
            piv["ganador"] = np.where(piv["AR_loss"] > piv["ESR1_gain"], "AR_loss", "ESR1_gain")
            resumen_mod = piv.reset_index().groupby(["tissue_group", "modalidad"]).agg(
                n_tracks=("track_name", "count"),
                AR_loss_wins=("ganador", lambda x: (x == "AR_loss").sum()),
                ESR1_gain_wins=("ganador", lambda x: (x == "ESR1_gain").sum()),
                mediana_ratio_AR_sobre_ESR1=("ratio_AR_sobre_ESR1", "median"),
            ).reset_index()
            resumen_mod.to_csv(ARCHIVO_SALIDA_RESUMEN_INS, index=False)
            print("\n=== RESUMEN DE IMPACTO DE INSERCIONES POR TEJIDO Y MODALIDAD ===")
            print(resumen_mod.to_string(index=False))

        print(f"\n[¡COMPLETADO!]. Tabla 5 guardada en: {ARCHIVO_SALIDA_TABLA5.name}")


if __name__ == "__main__":
    main()
