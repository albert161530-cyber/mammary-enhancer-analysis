#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evalúa candidatos de deleción humana con AlphaGenome sin ejecutar ni modificar Genoma1.py.

Entrada esperada: deleciones_todos_loci.csv, generado previamente por Genoma1.py.
Por defecto valida/reconstruye candidatos y ejecuta el filtro de motivos C++.
La API de AlphaGenome se invoca únicamente con --run-alphagenome.

Para un hueco humano en el MAF, la secuencia ancestral se representa como una inserción
respecto a hg38: REF = base ancla humana; ALT = base ancla + secuencia ancestral.
Los eventos del MAF son candidatos evolutivos, no deleciones confirmadas experimentalmente.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import gzip
import pickle
import shutil
import json
import math
import os
import re
import importlib.util
from collections import defaultdict
import ssl
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any


CPP_KERNELS_SOURCE = r'''
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <cctype>

namespace py = pybind11;

inline char complement_char(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default: return c;
    }
}

std::string reverse_complement_cpp(const std::string& seq) {
    std::string rev(seq.rbegin(), seq.rend());
    for (char& c : rev) c = complement_char(c);
    return rev;
}

double puntuacion_pwm_cpp(const std::string& window,
                          const std::vector<std::unordered_map<char, double>>& matrix,
                          double min_score, double max_score) {
    if (window.size() != matrix.size()) return -999.0;
    double score = 0.0;
    for (size_t i = 0; i < window.size(); ++i) {
        char b = static_cast<char>(std::toupper(static_cast<unsigned char>(window[i])));
        auto it = matrix[i].find(b);
        if (it == matrix[i].end()) return -999.0;
        score += it->second;
    }
    if (max_score <= min_score) return 0.0;
    return (score - min_score) / (max_score - min_score);
}

py::dict mejor_pwm_en_intervalo_cpp(const std::string& sequence,
                                    int region_start, int region_end,
                                    const std::vector<std::unordered_map<char, double>>& matrix,
                                    double min_score, double max_score, int width) {
    if (width <= 0 || static_cast<int>(matrix.size()) != width || sequence.empty()) return py::dict();
    region_start = std::max(0, region_start);
    region_end = std::min(static_cast<int>(sequence.size()), region_end);
    if (region_end <= region_start) return py::dict();
    int lo = std::max(0, region_start - width + 1);
    int hi = std::min(region_end - 1, static_cast<int>(sequence.size()) - width);
    bool found = false;
    double best_score = -1.0;
    int best_start = 0;
    std::string best_orientation = "+", best_sequence;
    for (int start = lo; start <= hi; ++start) {
        std::string window = sequence.substr(start, width);
        for (const std::string& orientation : {"+", "-"}) {
            std::string oriented = orientation == "-" ? reverse_complement_cpp(window) : window;
            double score = puntuacion_pwm_cpp(oriented, matrix, min_score, max_score);
            if (score < -100.0) continue;
            if (!found || score > best_score) {
                found = true;
                best_score = score;
                best_start = start;
                best_orientation = orientation;
                best_sequence = oriented;
            }
        }
    }
    if (!found) return py::dict();
    py::dict result;
    result["score"] = best_score;
    result["start"] = best_start;
    result["orientacion"] = best_orientation;
    result["secuencia"] = best_sequence;
    return result;
}

PYBIND11_MODULE(genoma_indel_kernels, m) {
    m.def("puntuacion_pwm_cpp", &puntuacion_pwm_cpp);
    m.def("mejor_pwm_en_intervalo_cpp", &mejor_pwm_en_intervalo_cpp);
    m.def("reverse_complement_cpp", &reverse_complement_cpp);
}
'''

BASES = set("ACGT")
ESPECIES = (
    ("chimpance", "secuencia_chimpance"),
    ("gorila", "secuencia_gorila"),
    ("orangutan", "secuencia_orangutan"),
    ("bonobo", "secuencia_bonobo"),
)
MODALIDADES = (
    ("RNA_SEQ", "RNA_SEQ", "rna_seq", 100_000),
    ("ATAC", "ATAC", "atac", 5_000),
    ("DNASE", "DNASE", "dnase", 5_000),
    ("CAGE", "CAGE", "cage", 50_000),
    ("PROCAP", "PROCAP", "procap", 50_000),
    ("CHIP_TF", "CHIP_TF", "chip_tf", 5_000),
    ("CHIP_HISTONE", "CHIP_HISTONE", "chip_histone", 5_000),
)
CAMPOS_TRACK = (
    "name", "strand", "ontology_curie", "biosample_name", "biosample_type",
    "gtex_tissue", "target", "histone_mark", "Assay title", "assay_title",
    "data_source", "gene_id", "gene_symbol",
)
ONTOLOGIAS_PREDETERMINADAS = [
    "CL:0002231", "CL:0002324", "CL:0002326", "CL:0002327", "CL:0002451",
    "CL:0002555", "CL:0002570", "CL:0002577", "CL:0002580", "CL:0002601",
    "CL:0002617", "CL:0002622", "CL:1000487", "CL:2000071", "CL:2000078",
    "EFO:0001203", "EFO:0001247",
    "UBERON:0000002", "UBERON:0000007", "UBERON:0000310", "UBERON:0000473",
    "UBERON:0000992", "UBERON:0000995", "UBERON:0001013", "UBERON:0001987",
    "UBERON:0002118", "UBERON:0002119", "UBERON:0002190", "UBERON:0002367",
    "UBERON:0005795", "UBERON:0008367", "UBERON:0010414", "UBERON:0014455",
]
ONTOLOGIAS_MAMARIAS_DEFAULT = [
    "UBERON:0008367", "CL:0002327", "CL:0002326",
    "CL:0002324", "CL:0002451", "CL:0002555",
    "EFO:0001203", "EFO:0001247",
]
ALPHAGENOME_PREDICTION_METHOD_VERSION = "multi_indel_coordinate_mapped_profiles_v6_ontology_filtered_local_cache"
LONGITUD_SECUENCIA = 2**20
MAX_SPAN_HAPLOTIPO = LONGITUD_SECUENCIA - 2 * 100_000 - 25_000
UMBRAL_PWM = 0.70
DELTA_PWM = 0.05
UMBRAL_CONSENSO = 0.95
DELTA_CONSENSO = 0.05
UMBRAL_RESPUESTA_NO_CANONICA = 0.83  # >=10/12 bases del dímero canónico
UMBRAL_ARE_SELECTIVO = 0.83
FLANCO_MOTIVOS = 60
URL_JASPAR_MEME = "https://jaspar2026.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
CATEGORIAS_HORMONALES = ("AR_loss", "ESR1_gain", "ERE_gain_shared", "ESR2_change")
RADIO_POR_MODALIDAD = {
    "RNA_SEQ": 100_000, "ATAC": 5_000, "DNASE": 5_000,
    "CAGE": 50_000, "PROCAP": 50_000, "CHIP_TF": 5_000,
    "CHIP_HISTONE": 5_000,
}

CAMPOS_CANDIDATO = [
    "variant_id", "gen_locus", "cromosoma", "punto_union_interbase_hg38_0based",
    "position_alphagenome_1based", "longitud_columnas_maf", "longitud_secuencia_ancestral_bp",
    "secuencia_ancestral_inferida", "especies_con_secuencia", "n_especies_con_secuencia",
    "soporte_minimo_consenso", "consenso_minimo_fraccion", "estado_maf_original",
    "estado_validacion", "motivo_validacion", "fuente_candidato", "archivo_maf",
    "bloque_maf", "columna_inicio_bloque_1based", "columna_fin_bloque_1based",
    "flancos_humanos_consecutivos", "longitud_indel_fuente_bp",
]

MIN_INDEL_MAF_DEFAULT = 1
MAX_INDEL_MAF_DEFAULT = 9
MAF_LOCUS_CONFIG = {
    "AR": {
        "chrom": "chrX",
        "files": ("AR_region_distal.maf", "alineamiento_ar_ampliado.maf", "alineamiento_ar.maf"),
        "tokens": ("androgenreceptor", "ar_region", "alineamiento_ar", "_ar_"),
    },
    "ESR1": {
        "chrom": "chr6",
        "files": ("alineamiento_esr1.maf", "ESR1_region.maf", "esr1.maf"),
        "tokens": ("esr1", "eralfa", "eralpha", "estrogenreceptoralpha"),
    },
    "ESR2": {
        "chrom": "chr14",
        "files": ("alineamiento_esr2.maf", "ESR2_region.maf", "esr2.maf"),
        "tokens": ("esr2", "erbeta", "estrogenreceptorbeta"),
    },
    "CYP19A1": {
        "chrom": "chr15",
        "files": ("alineamiento_cyp19a1.maf", "CYP19A1_region.maf", "cyp19a1.maf"),
        "tokens": ("cyp19a1", "aromatase", "aromatasa"),
    },
}
META_CANDIDATO_FIELDS = (
    "fuente_candidato", "archivo_maf", "bloque_maf", "columna_inicio_bloque_1based",
    "columna_fin_bloque_1based", "flancos_humanos_consecutivos", "longitud_indel_fuente_bp",
)

CAMPOS_EVENTO = [
    "variant_id", "gen_locus", "cromosoma", "position_alphagenome_1based",
    "punto_union_interbase_hg38_0based", "longitud_secuencia_ancestral_bp",
    "secuencia_ancestral_inferida", "base_ancla_hg38", "ref_alphagenome",
    "alt_alphagenome", "n_modalidades_evaluadas", "n_tracks_evaluados",
    "estado", "error", "nota",
]

CAMPOS_TRACK = [
    "variant_id", "gen_locus", "cromosoma", "position_alphagenome_1based",
    "longitud_secuencia_ancestral_bp", "secuencia_ancestral_inferida", "modalidad",
    "track_index", "track_name", "strand", "ontology_curie", "biosample_name",
    "biosample_type", "gtex_tissue", "target", "histone_mark", "assay_title",
    "data_source", "gene_id", "gene_symbol", "radio_local_bp", "n_bins_locales",
    "senal_humana_suma", "senal_ancestral_suma", "efecto_abs_L1_humano_vs_ancestral",
    "efecto_firmado_humano_menos_ancestral", "efecto_firmado_ancestral_menos_humano",
    "diferencia_abs_media_por_bin", "diferencia_abs_max_bin", "direccion_interpretacion",
]


def limpiar(value: Any) -> str:
    if value is None:
        return ""
    result = str(value).strip()
    if result.lower() in {"nan", "none", "<na>", "null"}:
        return ""
    return result


def normalizar_cromosoma(value: Any) -> str:
    chrom = limpiar(value)
    if not chrom:
        raise ValueError("cromosoma vacío")
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    chrom = "chr" + chrom[3:].upper()
    if chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}:
        raise ValueError(f"cromosoma no admitido: {value!r}")
    return chrom


def entero(value: Any, nombre: str) -> int:
    try:
        number = int(float(limpiar(value)))
    except (TypeError, ValueError):
        raise ValueError(f"{nombre} no es un entero válido: {value!r}") from None
    return number


def leer_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(
            f"No existe {path}. Coloca junto al script deleciones_todos_loci.csv, "
            "generado previamente por Genoma1.py; este script no ejecuta Genoma1.py."
        )
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        required = {
            "gen_locus", "cromosoma", "punto_union_interbase_hg38_0based",
            "longitud_columnas", *(col for _, col in ESPECIES),
        }
        missing = sorted(required - fields)
        if missing:
            raise ValueError("Faltan columnas en el CSV de deleciones: " + ", ".join(missing))
        return [dict(row) for row in reader]


def inferir_consenso(row: dict[str, str], min_support: int = 2) -> tuple[str, int, float, list[str]]:
    seqs: list[tuple[str, str]] = []
    for especie, col in ESPECIES:
        seq = limpiar(row.get(col, "")).upper()
        if seq:
            seqs.append((especie, seq))
    if len(seqs) < min_support:
        raise ValueError(f"solo hay secuencias de {len(seqs)} especies; se requieren {min_support}")
    lengths = {len(seq) for _, seq in seqs}
    if len(lengths) != 1:
        raise ValueError(f"secuencias MAF de longitudes distintas: {sorted(lengths)}")
    expected = entero(row.get("longitud_columnas"), "longitud_columnas")
    aligned_length = next(iter(lengths))
    if expected != aligned_length:
        raise ValueError(
            f"longitud_columnas={expected}, pero las secuencias tienen {aligned_length} columnas"
        )

    consensus: list[str] = []
    support_counts: list[int] = []
    support_fractions: list[float] = []
    for col_index in range(aligned_length):
        bases = [seq[col_index] for _, seq in seqs if seq[col_index] in BASES]
        counts = Counter(bases)
        if len(bases) < min_support or not counts:
            raise ValueError(f"columna {col_index + 1}: soporte insuficiente para reconstrucción")
        ordered = counts.most_common()
        if ordered[0][1] < min_support:
            raise ValueError(f"columna {col_index + 1}: ninguna base tiene soporte suficiente")
        if len(ordered) > 1 and ordered[0][1] == ordered[1][1]:
            raise ValueError(f"columna {col_index + 1}: empate ancestral {dict(counts)}")
        consensus.append(ordered[0][0])
        support_counts.append(ordered[0][1])
        support_fractions.append(ordered[0][1] / len(bases))

    return (
        "".join(consensus),
        min(support_counts),
        min(support_fractions),
        [especie for especie, _ in seqs],
    )


def crear_id(gen: str, chrom: str, position: int, sequence: str) -> str:
    payload = f"{gen}|{chrom}|{position}|{sequence}".encode("utf-8")
    digest = hashlib.sha256(payload).hexdigest()[:16]
    return f"DEL_{chrom}_{position}_{digest}"


def preparar_candidatos(rows: list[dict[str, str]], min_support: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    validos: list[dict[str, Any]] = []
    validation: list[dict[str, Any]] = []
    seen: set[str] = set()

    for input_row, row in enumerate(rows, start=2):
        result: dict[str, Any] = {field: "" for field in CAMPOS_CANDIDATO}
        for meta_field in META_CANDIDATO_FIELDS:
            result[meta_field] = limpiar(row.get(meta_field, ""))
        result["fila_entrada"] = input_row
        result["gen_locus"] = limpiar(row.get("gen_locus", ""))
        result["estado_maf_original"] = limpiar(row.get("estado", ""))
        try:
            chrom = normalizar_cromosoma(row.get("cromosoma", ""))
            point = entero(row.get("punto_union_interbase_hg38_0based"), "punto_union_interbase_hg38_0based")
            if point < 2:
                raise ValueError("punto de unión demasiado cercano al inicio del cromosoma para crear una variante anclada")
            ancestral, support, fraction, species = inferir_consenso(row, min_support=min_support)
            if not ancestral or set(ancestral) - BASES:
                raise ValueError("secuencia ancestral reconstruida no válida")
            gen = limpiar(row.get("gen_locus", "")) or "unknown"
            variant_id = crear_id(gen, chrom, point, ancestral)
            result.update({
                "variant_id": variant_id,
                "cromosoma": chrom,
                "punto_union_interbase_hg38_0based": point,
                "position_alphagenome_1based": point,
                "longitud_columnas_maf": entero(row.get("longitud_columnas"), "longitud_columnas"),
                "longitud_secuencia_ancestral_bp": len(ancestral),
                "secuencia_ancestral_inferida": ancestral,
                "especies_con_secuencia": ";".join(species),
                "n_especies_con_secuencia": len(species),
                "soporte_minimo_consenso": support,
                "consenso_minimo_fraccion": round(fraction, 4),
                "estado_maf_original": limpiar(row.get("estado", "")),
                "estado_validacion": "lista_para_validar_referencia",
                "motivo_validacion": "consenso ancestral reconstruido; requiere comprobar base ancla hg38 y alineamiento",
            })
            if variant_id in seen:
                result["estado_validacion"] = "duplicado_omitido"
                result["motivo_validacion"] = "mismo locus, coordenada y secuencia ancestral ya presente"
            else:
                seen.add(variant_id)
                validos.append({**row, **result, "gen_locus": gen})
        except Exception as exc:
            result["estado_validacion"] = "rechazada_reconstruccion"
            result["motivo_validacion"] = f"{type(exc).__name__}: {exc}"
        validation.append(result)
    return validos, validation



def iterar_bloques_maf(path: Path):
    """Itera bloques MAF sin cargar el archivo completo en memoria."""
    block: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line.strip():
                if block:
                    yield block
                    block = []
            elif not line.startswith("#") and not line.startswith("track "):
                block.append(line)
    if block:
        yield block


def _maf_species(src: str) -> str | None:
    prefix = src.split(".", 1)[0].lower()
    full = src.lower()
    if prefix in {"hg38", "human", "homo_sapiens"} or "homo_sapiens" in full:
        return "humano"
    if prefix.startswith("pantro5") or prefix in {"chimp", "chimpanzee"} or "pan_troglodytes" in full:
        return "chimpance"
    if prefix.startswith("gorgor5") or prefix == "gorilla" or "gorilla_gorilla" in full:
        return "gorila"
    if prefix.startswith("ponabe2") or prefix == "orangutan" or "pongo_abelii" in full:
        return "orangutan"
    if prefix.startswith("panpan2") or prefix == "bonobo" or "pan_paniscus" in full:
        return "bonobo"
    return None


def _maf_chrom(src: str) -> str:
    parts = src.split(".", 1)
    if len(parts) < 2:
        return ""
    return parts[1].split(":", 1)[0]


def _parsear_bloque_maf_para_indels(lines: list[str], expected_chrom: str) -> dict[str, Any] | None:
    rows: list[dict[str, Any]] = []
    for line in lines:
        fields = line.split()
        if not fields or fields[0] != "s" or len(fields) < 7:
            continue
        try:
            row = {
                "src": fields[1], "chrom": _maf_chrom(fields[1]),
                "start": int(fields[2]), "size": int(fields[3]), "strand": fields[4],
                "src_size": int(fields[5]), "text": fields[6].upper(),
                "species": _maf_species(fields[1]),
            }
        except (ValueError, IndexError):
            continue
        rows.append(row)
    if not rows:
        return None
    aln_len = len(rows[0]["text"])
    if aln_len == 0 or any(len(row["text"]) != aln_len for row in rows):
        return None
    humans = [
        row for row in rows
        if row["species"] == "humano"
        and row["src"].split(".", 1)[0].lower() == "hg38"
        and row["chrom"].lower() == expected_chrom.lower()
    ]
    if len(humans) != 1:
        return None
    human = humans[0]
    species_rows: dict[str, dict[str, Any]] = {}
    for row in rows:
        species = row["species"]
        if species is None:
            continue
        if species not in species_rows or row["size"] > species_rows[species]["size"]:
            species_rows[species] = row
    if not all(x in species_rows for x in ("chimpance", "gorila", "orangutan")):
        return None

    coords: list[int | None] = [None] * aln_len
    offset = 0
    for i, base in enumerate(human["text"]):
        if base == "-":
            continue
        coords[i] = human["start"] + offset if human["strand"] == "+" else human["src_size"] - human["start"] - 1 - offset
        offset += 1

    complement = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn-", "TGCAYRKMVHDBNtgcayrkmvhdbn-")
    for row in species_rows.values():
        row["text_hg38"] = row["text"].translate(complement) if human["strand"] == "-" else row["text"]
    return {"human": human, "rows": species_rows, "coords": coords, "length": aln_len}


def _canonizar_locus(value: str) -> str | None:
    compact = re.sub(r"[^a-z0-9]+", "", limpiar(value).lower())
    if "cyp19a1" in compact or "aromatase" in compact or "aromatasa" in compact:
        return "CYP19A1"
    if "esr1" in compact or "estrogenreceptoralpha" in compact or "eralpha" in compact or "eralfa" in compact:
        return "ESR1"
    if "esr2" in compact or "estrogenreceptorbeta" in compact or "erbeta" in compact:
        return "ESR2"
    if compact == "ar" or "androgenreceptor" in compact:
        return "AR"
    tokens = set(re.split(r"[^a-z0-9]+", limpiar(value).lower()))
    if "ar" in tokens:
        return "AR"
    return None


def resolver_mafs_suplementarios(maf_dir: Path) -> tuple[list[tuple[str, Path, str]], list[dict[str, Any]]]:
    """Encuentra como máximo un MAF por locus; no ejecuta ni importa Genoma1.py."""
    all_mafs = sorted(maf_dir.glob("*.maf")) if maf_dir.exists() else []
    resolved: list[tuple[str, Path, str]] = []
    status_rows: list[dict[str, Any]] = []
    used: set[Path] = set()
    for locus, config in MAF_LOCUS_CONFIG.items():
        chosen = None
        for name in config["files"]:
            exact = [p for p in all_mafs if p.name.lower() == name.lower()]
            if exact:
                chosen = exact[0]
                break
        if chosen is None:
            candidates = []
            for path in all_mafs:
                stem = re.sub(r"[^a-z0-9]+", "_", path.stem.lower()).strip("_")
                compact = re.sub(r"[^a-z0-9]+", "", path.stem.lower())
                tokens = set(stem.split("_"))
                if locus == "AR":
                    matched = ("androgenreceptor" in compact or "ar_region" in stem or "alineamiento_ar" in stem or stem == "ar" or stem.startswith("ar_"))
                else:
                    matched = any(token in compact for token in config["tokens"])
                if matched:
                    candidates.append(path)
            if len(candidates) == 1:
                chosen = candidates[0]
            elif len(candidates) > 1:
                prioritized = [p for p in candidates if locus.lower() in p.stem.lower()]
                if len(prioritized) == 1:
                    chosen = prioritized[0]
                else:
                    status_rows.append({"gen_locus": locus, "archivo_maf": "", "estado": "MAF_AMBIGUO", "detalle": ";".join(p.name for p in candidates)})
                    continue
        if chosen is None:
            status_rows.append({"gen_locus": locus, "archivo_maf": "", "estado": "MAF_NO_ENCONTRADO", "detalle": f"No hallé MAF para {locus} en {maf_dir}"})
            continue
        chosen = chosen.resolve()
        if chosen in used:
            status_rows.append({"gen_locus": locus, "archivo_maf": chosen.name, "estado": "MAF_REUTILIZADO_OMITIDO", "detalle": "El mismo archivo ya fue asignado a otro locus"})
            continue
        used.add(chosen)
        resolved.append((locus, chosen, config["chrom"]))
        status_rows.append({"gen_locus": locus, "archivo_maf": chosen.name, "estado": "MAF_SELECCIONADO", "detalle": str(chosen)})
    return resolved, status_rows


def detectar_indels_pequenos_maf(
    maf_dir: Path, min_bp: int = MIN_INDEL_MAF_DEFAULT, max_bp: int = MAX_INDEL_MAF_DEFAULT,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Recupera huecos humanos alineados de 1–9 columnas con consenso de primates.

    Solo detecta el mismo tipo de evento que produce deleciones_todos_loci.csv:
    bases presentes en primates ancestrales y ausentes (gap) en hg38. No intenta
    llamar inserciones humanas ni segmentos MAF no alineados.
    """
    seeds: list[dict[str, Any]] = []
    scan_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    maf_inputs, resolver_status = resolver_mafs_suplementarios(maf_dir)
    stats = {row.get("gen_locus", ""): {"bloques": 0, "bloques_utiles": 0, "huecos_1_9": 0, "candidatos": 0, "rechazados": 0} for row in resolver_status}

    for locus, path, chrom in maf_inputs:
        item_stats = {"gen_locus": locus, "archivo_maf": path.name, "ruta_maf": str(path), "estado": "PROCESANDO", "bloques": 0, "bloques_utiles": 0, "huecos_1_9": 0, "candidatos_generados": 0, "rechazados_por_calidad": 0, "detalle": ""}
        try:
            for block_no, lines in enumerate(iterar_bloques_maf(path), start=1):
                item_stats["bloques"] += 1
                block = _parsear_bloque_maf_para_indels(lines, chrom)
                if block is None:
                    continue
                item_stats["bloques_utiles"] += 1
                human = block["human"]
                rows = block["rows"]
                coords = block["coords"]
                aligned_human = human["text"]
                length = block["length"]
                i = 0
                while i < length:
                    if aligned_human[i] != "-":
                        i += 1
                        continue
                    start_col = i
                    while i < length and aligned_human[i] == "-":
                        i += 1
                    end_col = i
                    run_len = end_col - start_col
                    if not (min_bp <= run_len <= max_bp):
                        continue

                    item_stats["huecos_1_9"] += 1
                    primary = ("chimpance", "gorila", "orangutan")
                    ancestral_columns: list[str] = []
                    support_ok = True
                    support_min = 99
                    consensus_fracs: list[float] = []
                    for col in range(start_col, end_col):
                        bases = [rows[sp]["text_hg38"][col] for sp in primary if rows[sp]["text_hg38"][col] in BASES]
                        counts = Counter(bases)
                        support_min = min(support_min, len(bases))
                        if len(bases) < 2 or not counts or counts.most_common(1)[0][1] < 2:
                            support_ok = False
                            break
                        ancestral_columns.append(counts.most_common(1)[0][0])
                        consensus_fracs.append(counts.most_common(1)[0][1] / len(bases))
                    left = next((j for j in range(start_col - 1, -1, -1) if coords[j] is not None), None)
                    right = next((j for j in range(end_col, length) if coords[j] is not None), None)
                    left_coord = coords[left] if left is not None else None
                    right_coord = coords[right] if right is not None else None
                    reason = ""
                    if not support_ok:
                        reason = "soporte_ancestral_insuficiente_o_sin_consenso_2_de_3"
                    elif left_coord is None or right_coord is None:
                        reason = "sin_dos_flancos_humanos_alineados"
                    elif abs(int(left_coord) - int(right_coord)) != 1:
                        reason = "flancos_humanos_no_consecutivos"

                    seq_by_species: dict[str, str] = {}
                    for species, field in ESPECIES:
                        row = rows.get(species)
                        if row is None:
                            seq_by_species[field] = ""
                            continue
                        segment = row["text_hg38"][start_col:end_col]
                        if human["strand"] == "-":
                            segment = segment[::-1]
                        seq_by_species[field] = segment

                    point0 = max(int(left_coord), int(right_coord)) if left_coord is not None and right_coord is not None else ""
                    record = {
                        "gen_locus": locus, "cromosoma": chrom,
                        "punto_union_interbase_hg38_0based": point0,
                        "longitud_columnas": run_len,
                        "secuencia_chimpance": seq_by_species.get("secuencia_chimpance", ""),
                        "secuencia_gorila": seq_by_species.get("secuencia_gorila", ""),
                        "secuencia_orangutan": seq_by_species.get("secuencia_orangutan", ""),
                        "secuencia_bonobo": seq_by_species.get("secuencia_bonobo", ""),
                        "estado": "hueco_humano_maf_1_9bp",
                        "fuente_candidato": "MAF_SUPLEMENTARIO_1_9BP",
                        "archivo_maf": path.name, "bloque_maf": block_no,
                        "columna_inicio_bloque_1based": start_col + 1,
                        "columna_fin_bloque_1based": end_col,
                        "flancos_humanos_consecutivos": (abs(int(left_coord) - int(right_coord)) == 1) if left_coord is not None and right_coord is not None else False,
                        "longitud_indel_fuente_bp": run_len,
                    }
                    scan_row = {
                        "gen_locus": locus, "cromosoma": chrom, "archivo_maf": path.name,
                        "bloque_maf": block_no, "columna_inicio_bloque_1based": start_col + 1,
                        "columna_fin_bloque_1based": end_col, "longitud_indel_bp": run_len,
                        "punto_union_interbase_hg38_0based": point0,
                        "flanco_izquierdo_hg38_0based": left_coord if left_coord is not None else "",
                        "flanco_derecho_hg38_0based": right_coord if right_coord is not None else "",
                        "soporte_minimo_primates": support_min if support_min != 99 else 0,
                        "consenso_minimo_fraccion": round(min(consensus_fracs), 4) if consensus_fracs else "",
                        "secuencia_ancestral_consenso_preliminar": "".join(ancestral_columns) if support_ok else "",
                        "secuencia_chimpance": seq_by_species.get("secuencia_chimpance", ""),
                        "secuencia_gorila": seq_by_species.get("secuencia_gorila", ""),
                        "secuencia_orangutan": seq_by_species.get("secuencia_orangutan", ""),
                        "secuencia_bonobo": seq_by_species.get("secuencia_bonobo", ""),
                        "estado": "LISTO_PARA_VALIDACION" if not reason else "RECHAZADO_CALIDAD_MAF",
                        "detalle": reason,
                    }
                    if reason:
                        item_stats["rechazados_por_calidad"] += 1
                    else:
                        record["secuencia_ancestral_consenso_preliminar"] = "".join(ancestral_columns)
                        seeds.append(record)
                        item_stats["candidatos_generados"] += 1
                    scan_rows.append(scan_row)
            item_stats["estado"] = "COMPLETADO"
        except Exception as exc:
            item_stats["estado"] = "ERROR_LECTURA_MAF"
            item_stats["detalle"] = f"{type(exc).__name__}: {exc}"
        file_rows.append(item_stats)

    for row in resolver_status:
        if row.get("estado") != "MAF_SELECCIONADO":
            file_rows.append({
                "gen_locus": row.get("gen_locus", ""), "archivo_maf": row.get("archivo_maf", ""),
                "ruta_maf": "", "estado": row.get("estado", ""), "bloques": 0,
                "bloques_utiles": 0, "huecos_1_9": 0, "candidatos_generados": 0,
                "rechazados_por_calidad": 0, "detalle": row.get("detalle", ""),
            })
    return seeds, scan_rows, file_rows


def combinar_candidatos_suplementarios(
    existentes: list[dict[str, Any]], suplementarios: list[dict[str, Any]],
    validacion_suplementaria: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, str]]:
    """Deduplica eventos y excluye posiciones con secuencias ancestrales conflictivas."""
    def key(row: dict[str, Any]) -> tuple[str, int, str]:
        return (
            normalizar_cromosoma(row["cromosoma"]),
            entero(row["punto_union_interbase_hg38_0based"], "punto_union_interbase_hg38_0based"),
            limpiar(row.get("secuencia_ancestral_inferida", "")).upper(),
        )

    existing_keys: set[tuple[str, int, str]] = set()
    existing_positions: dict[tuple[str, int], set[str]] = defaultdict(set)
    for row in existentes:
        try:
            k = key(row)
        except Exception:
            continue
        existing_keys.add(k)
        existing_positions[(k[0], k[1])].add(k[2])

    by_position: dict[tuple[str, int], set[str]] = defaultdict(set)
    for row in suplementarios:
        try:
            k = key(row)
        except Exception:
            continue
        by_position[(k[0], k[1])].add(k[2])

    validation_by_id: dict[str, dict[str, Any]] = {}
    for row in validacion_suplementaria:
        vid = limpiar(row.get("variant_id", ""))
        if vid:
            validation_by_id[vid] = row
    final_additions: list[dict[str, Any]] = []
    final_status: dict[str, str] = {}
    seen_keys: set[tuple[str, int, str]] = set()
    for row in suplementarios:
        vid = limpiar(row.get("variant_id", ""))
        k = key(row)
        pos_key = (k[0], k[1])
        status = ""
        if len(by_position[pos_key]) > 1:
            status = "conflicto_entre_alineamientos_MAF_secuencia_ancestral_distinta"
        elif k in existing_keys:
            status = "duplicado_de_CSV_existente"
        elif pos_key in existing_positions and any(seq != k[2] for seq in existing_positions[pos_key]):
            status = "conflicto_con_CSV_existente_en_misma_posicion"
        elif k in seen_keys:
            status = "duplicado_suplementario_omitido"
        else:
            final_additions.append(row)
            seen_keys.add(k)
            status = "candidato_suplementario_agregado"
        final_status[vid] = status
        vr = validation_by_id.get(vid)
        if vr is not None and status != "candidato_suplementario_agregado":
            vr["estado_validacion"] = status
            vr["motivo_validacion"] = status
    for row in suplementarios:
        if not row.get("variant_id"):
            continue
        final_status.setdefault(row["variant_id"], "candidato_suplementario_agregado")
    return final_additions, validacion_suplementaria, final_status

def escribir_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    temp = Path(temp_name)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: csv_value(row.get(key, "")) for key in fields})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temp), str(path))
    finally:
        try:
            temp.unlink(missing_ok=True)
        except OSError:
            pass


def csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False)
    return value


def guardar_json_atomico(path: Path, data: dict[str, Any]) -> None:
    """Guarda el checkpoint con temporal único y reintentos ante bloqueos transitorios de Windows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    temp = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            json.dump(data, handle, ensure_ascii=False, separators=(",", ":"))
            handle.flush()
            os.fsync(handle.fileno())
        delays = (0.2, 0.5, 1.0, 2.0, 3.0)
        for attempt in range(len(delays) + 1):
            try:
                os.replace(str(temp), str(path))
                return
            except PermissionError:
                if attempt >= len(delays):
                    raise
                time.sleep(delays[attempt])
    finally:
        try:
            temp.unlink(missing_ok=True)
        except OSError:
            pass


def cargar_checkpoint(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"schema_version": 1, "tasks": {}, "failure_history": []}
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
        if not isinstance(value, dict):
            raise ValueError("el contenido no es un objeto JSON")
        value.setdefault("schema_version", 1)
        value.setdefault("tasks", {})
        value.setdefault("failure_history", [])
        if not isinstance(value["tasks"], dict):
            value["tasks"] = {}
        if not isinstance(value["failure_history"], list):
            value["failure_history"] = []
        return value
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise RuntimeError(f"No se pudo leer el checkpoint {path}: {exc}") from exc


def guardar_estado_checkpoint(out_dir: Path, checkpoint: dict[str, Any]) -> None:
    checkpoint_path = out_dir / "evaluacion_deleciones_alphagenome_checkpoint.json"
    guardar_json_atomico(checkpoint_path, checkpoint)

    # El JSON es el checkpoint autoritativo. Los CSV son vistas auxiliares y
    # nunca deben impedir una predicción si Excel, un sincronizador o el antivirus
    # mantiene bloqueado el archivo de destino en Windows.
    tasks = []
    for task_id, value in checkpoint.get("tasks", {}).items():
        tasks.append({"task_id": task_id, **value})
    task_fields = ["task_id", "status", "attempts", "updated_at", "completed_at", "cache_file", "last_error", "blocked_by", "error_history"]
    failures = checkpoint.get("failure_history", [])
    failure_fields = ["timestamp", "task_id", "attempt", "error"]
    auxiliary_csvs = [
        (out_dir / "evaluacion_deleciones_alphagenome_tareas.csv", tasks, task_fields),
        (out_dir / "evaluacion_deleciones_alphagenome_historial_fallos.csv", failures, failure_fields),
    ]
    for csv_path, rows, fields in auxiliary_csvs:
        try:
            escribir_csv(csv_path, rows, fields)
        except OSError as exc:
            print(
                f"  [AVISO CHECKPOINT] No se pudo actualizar el CSV auxiliar {csv_path.name}: "
                f"{type(exc).__name__}: {exc}. El estado sigue guardado en el JSON; "
                "cierra ese CSV si está abierto en Excel para poder actualizarlo.",
                file=sys.stderr, flush=True,
            )


def registrar_tarea_bloqueada(
    out_dir: Path, checkpoint: dict[str, Any], task_id: str, blocked_by: list[str],
) -> None:
    now = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    old = checkpoint.setdefault("tasks", {}).get(task_id, {})
    old.update({
        "status": "blocked_by_failed_individuals",
        "blocked_by": blocked_by,
        "updated_at": now,
        "last_error": "",
    })
    checkpoint["tasks"][task_id] = old
    guardar_estado_checkpoint(out_dir, checkpoint)


def filtrar_bloques_por_ontologias(
    blocks: dict[str, dict[str, Any]], ontology_terms: list[str],
) -> dict[str, dict[str, Any]]:
    """Reduce una caché de predicción más amplia a los tracks solicitados actualmente."""
    filtered: dict[str, dict[str, Any]] = {}
    wanted = set(ontology_terms)
    for modality, block in blocks.items():
        metadata = block.get("metadata")
        if metadata is None or "ontology_curie" not in getattr(metadata, "columns", []):
            continue
        mask = metadata["ontology_curie"].astype(str).isin(wanted).to_numpy()
        if not mask.any():
            continue
        indices = [i for i, keep in enumerate(mask) if keep]
        new_metadata = metadata.iloc[indices].copy().reset_index(drop=True)
        new_values = block["values"][:, mask]
        filtered[modality] = {
            "values": new_values,
            "metadata": new_metadata,
            "signature": tuple(firmas_metadata(new_metadata)),
            "start": int(block["start"]),
            "resolution": int(block["resolution"]),
        }
    return filtered


def predecir_bloques_con_checkpoint(
    model: Any, sequence: str, interval: Any, requested: list[Any], ontology_terms: list[str],
    out_dir: Path, cache_dir: Path, checkpoint: dict[str, Any], task_id: str,
) -> dict[str, dict[str, Any]]:
    """Ejecuta una predicción o carga su resultado persistido. Cada llamada es reanudable."""
    sequence_hash = hashlib.sha256(sequence.encode("ascii", errors="ignore")).hexdigest()
    requested_names = [getattr(item, "name", str(item)) for item in requested]
    identity = {
        "task_id": task_id,
        "sequence_sha256": sequence_hash,
        "chromosome": str(getattr(interval, "chromosome", "")),
        "start": int(getattr(interval, "start", 0)),
        "end": int(getattr(interval, "end", 0)),
        "requested_outputs": requested_names,
        "ontology_terms": list(ontology_terms),
    }
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"prediction_{digest}.pkl.gz"
    tasks = checkpoint.setdefault("tasks", {})
    entry = tasks.setdefault(task_id, {"attempts": 0, "error_history": []})

    if cache_path.is_file():
        try:
            with gzip.open(cache_path, "rb") as handle:
                blocks = pickle.load(handle)
            if isinstance(blocks, dict):
                if entry.get("status") != "completed":
                    entry["status"] = "completed"
                    entry["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
                    entry["updated_at"] = entry["completed_at"]
                    entry["cache_file"] = str(cache_path.relative_to(out_dir))
                    entry["last_error"] = ""
                    entry["blocked_by"] = []
                    guardar_estado_checkpoint(out_dir, checkpoint)
                print(f"  [CHECKPOINT] reutilizada: {task_id}", flush=True)
                return blocks
        except Exception as exc:
            entry["status"] = "failed"
            entry["last_error"] = f"No se pudo leer caché: {type(exc).__name__}: {exc}"
            entry.setdefault("error_history", []).append(entry["last_error"])

    # Si el checkpoint ya tenía una predicción de referencia con más ontologías,
    # filtra sus tracks mamarios y crea la nueva clave de caché sin repetir la llamada.
    previous_cache_rel = entry.get("cache_file", "")
    if previous_cache_rel:
        previous_cache_path = out_dir / Path(previous_cache_rel)
        if previous_cache_path.is_file() and previous_cache_path != cache_path:
            try:
                with gzip.open(previous_cache_path, "rb") as handle:
                    previous_blocks = pickle.load(handle)
                if isinstance(previous_blocks, dict):
                    filtered_blocks = filtrar_bloques_por_ontologias(previous_blocks, ontology_terms)
                    if filtered_blocks:
                        temp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
                        with gzip.open(temp_path, "wb", compresslevel=3) as handle:
                            pickle.dump(filtered_blocks, handle, protocol=pickle.HIGHEST_PROTOCOL)
                            handle.flush()
                        temp_path.replace(cache_path)
                        entry["status"] = "completed"
                        entry["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
                        entry["updated_at"] = entry["completed_at"]
                        entry["cache_file"] = str(cache_path.relative_to(out_dir))
                        entry["last_error"] = ""
                        entry["blocked_by"] = []
                        guardar_estado_checkpoint(out_dir, checkpoint)
                        print(f"  [CHECKPOINT] referencia previa filtrada por ontologías: {task_id}", flush=True)
                        return filtered_blocks
            except Exception as exc:
                print(
                    f"  [AVISO CHECKPOINT] No se pudo reutilizar la referencia previa filtrando ontologías: "
                    f"{type(exc).__name__}: {exc}. Se solicitará una predicción nueva.",
                    file=sys.stderr, flush=True,
                )

    entry["attempts"] = int(entry.get("attempts", 0)) + 1
    entry["status"] = "running"
    entry["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    entry["cache_file"] = str(cache_path.relative_to(out_dir))
    entry["blocked_by"] = []
    entry["last_error"] = ""
    guardar_estado_checkpoint(out_dir, checkpoint)

    try:
        output = model.predict_sequence(
            sequence=sequence,
            interval=interval,
            requested_outputs=requested,
            ontology_terms=ontology_terms,
        )
        blocks = output_blocks(output)
        temp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
        with gzip.open(temp_path, "wb", compresslevel=3) as handle:
            pickle.dump(blocks, handle, protocol=pickle.HIGHEST_PROTOCOL)
            handle.flush()
        temp_path.replace(cache_path)
        entry["status"] = "completed"
        entry["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        entry["updated_at"] = entry["completed_at"]
        entry["cache_file"] = str(cache_path.relative_to(out_dir))
        entry["last_error"] = ""
        entry["blocked_by"] = []
        guardar_estado_checkpoint(out_dir, checkpoint)
        return blocks
    except Exception as exc:
        message = f"{type(exc).__name__}: {exc}"
        entry["status"] = "failed"
        entry["last_error"] = message
        entry["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        entry.setdefault("error_history", []).append(message)
        # Preserve a chronological record, including repeated failures of the same task.
        checkpoint.setdefault("failure_history", []).append({
            "timestamp": entry["updated_at"], "task_id": task_id,
            "attempt": entry["attempts"], "error": message,
        })
        guardar_estado_checkpoint(out_dir, checkpoint)
        raise

def valor_score_serializable(value: Any) -> Any:
    """Convierte valores escalares de pandas/numpy a tipos estables para caché/CSV."""
    if value is None:
        return ""
    if hasattr(value, "item"):
        try:
            value = value.item()
        except Exception:
            pass
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    text = str(value)
    if text in {"nan", "NaN", "<NA>", "None", "NaT"}:
        return ""
    if isinstance(value, (str, int, float, bool)):
        return value
    return text


def score_a_numero(value: Any) -> float | None:
    value = valor_score_serializable(value)
    if value == "":
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def archivar_resultados_metodo_anterior(out_dir: Path) -> None:
    """Preserva resultados calculados con una versión anterior del método de puntuación."""
    names = (
        "evaluacion_deleciones_alphagenome_tracks.csv",
        "evaluacion_deleciones_impacto_conjunto_por_track.csv",
        "evaluacion_deleciones_comparacion_AR_vs_ESR1.csv",
        "evaluacion_deleciones_comparacion_individual_AR_vs_ESR1.csv",
        "evaluacion_deleciones_haplotipos_regiones.csv",
        "evaluacion_deleciones_alphagenome_eventos.csv",
    )
    for name in names:
        source = out_dir / name
        if not source.exists():
            continue
        archive = source.with_name(source.stem + ".legacy_previous_scoring_method" + source.suffix)
        if archive.exists():
            archive = source.with_name(source.stem + f".legacy_previous_scoring_method_{time.strftime('%Y%m%d_%H%M%S')}" + source.suffix)
        try:
            source.replace(archive)
            print(f"  [MÉTODO] salida anterior preservada: {archive.name}", flush=True)
        except OSError as exc:
            print(
                f"  [AVISO] No se pudo archivar {source.name}: {exc}. Cierra el CSV si está abierto en Excel; "
                "se intentará sobrescribir al generar resultados nuevos.",
                file=sys.stderr, flush=True,
            )


def registrar_conjunto_no_evaluado(
    out_dir: Path, checkpoint: dict[str, Any], task_id: str, reason: str,
) -> None:
    entry = checkpoint.setdefault("tasks", {}).setdefault(task_id, {"attempts": 0, "error_history": []})
    entry.update({
        "status": "not_evaluated_requires_multivariant_indel_alignment",
        "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "last_error": reason,
        "blocked_by": [],
    })
    guardar_estado_checkpoint(out_dir, checkpoint)

def obtener_api_key(key_file: Path) -> str:
    if key_file.exists():
        value = key_file.read_text(encoding="utf-8").strip()
        if value:
            return value
    value = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if value:
        return value
    raise RuntimeError(
        f"No se encontró API key. Ponla en {key_file.name} o define ALPHAGENOME_API_KEY."
    )


def crear_contexto_ssl() -> ssl.SSLContext:
    """Usa certificados del sistema y añade el bundle CA de certifi si está instalado.

    No desactiva la verificación TLS. Esto ayuda en entornos Windows donde el
    almacén de certificados usado por Python no contiene toda la cadena pública.
    Si la conexión pasa por un proxy/antivirus que inspecciona HTTPS, su CA raíz
    deberá instalarse/configurarse por separado.
    """
    context = ssl.create_default_context()
    try:
        import certifi
        context.load_verify_locations(cafile=certifi.where())
    except ImportError:
        pass
    return context


def obtener_json(url: str, timeout: int = 60) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": "evaluar-deleciones/1.1"})
    with urllib.request.urlopen(request, timeout=timeout, context=crear_contexto_ssl()) as response:
        return json.loads(response.read().decode("utf-8"))


def obtener_base_ancla_hg38(chrom: str, position_1based: int) -> str:
    start0 = position_1based - 1
    end0 = position_1based
    query = urllib.parse.urlencode({"genome": "hg38", "chrom": chrom, "start": start0, "end": end0})
    url = "https://api.genome.ucsc.edu/getData/sequence?" + query
    data = obtener_json(url)
    base = limpiar(data.get("dna", "")).upper()
    if len(base) != 1 or base not in BASES:
        raise ValueError(f"UCSC no devolvió una base ancla A/C/G/T válida en {chrom}:{position_1based}: {base!r}")
    return base


def normalizar_array(track_data: Any) -> tuple[Any, Any, int, int]:
    import numpy as np

    # Mantener la precisión de origen (normalmente float32). La conversión
    # incondicional a float64 duplica la RAM de perfiles genómicos grandes.
    raw_values = np.asarray(track_data.values)
    values = raw_values if raw_values.dtype == np.float32 else raw_values.astype(np.float32, copy=False)
    metadata = track_data.metadata
    n_tracks = len(metadata)
    if values.ndim == 1:
        if n_tracks == 1:
            values = values.reshape(-1, 1)
        elif len(values) == n_tracks:
            values = values.reshape(1, -1)
        else:
            raise ValueError(f"vector inesperado; bins={len(values)}, tracks={n_tracks}")
    if values.ndim != 2:
        raise ValueError(f"se esperaba matriz bins×tracks y se recibió {values.shape}")
    if values.shape[1] != n_tracks and values.shape[0] == n_tracks:
        values = values.T
    if values.shape[1] != n_tracks:
        raise ValueError(f"no coinciden columnas de predicción {values.shape} con {n_tracks} tracks")
    if not np.isfinite(values).all():
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
    interval = getattr(track_data, "interval", None)
    if interval is None:
        raise ValueError("AlphaGenome devolvió una modalidad sin intervalo")
    resolution = int(getattr(track_data, "resolution", 1))
    if resolution < 1:
        raise ValueError(f"resolución inválida: {resolution}")
    return values, metadata, int(interval.start), resolution


def fila_metadata(metadata: Any, index: int) -> dict[str, str]:
    row = metadata.iloc[index]
    return {field: limpiar(row.get(field, "")) for field in CAMPOS_TRACK}


def firmas_metadata(metadata: Any) -> list[tuple[str, ...]]:
    fields = [field for field in CAMPOS_TRACK if field in getattr(metadata, "columns", [])]
    if not fields:
        fields = [field for field in ("name", "strand", "ontology_curie", "biosample_name") if field in getattr(metadata, "columns", [])]
    return [tuple(limpiar(metadata.iloc[i].get(field, "")) for field in fields) for i in range(len(metadata))]


def evaluar_modalidad(
    variant: dict[str, Any], modality: str, attr: str, radius: int, output: Any,
) -> list[dict[str, Any]]:
    import numpy as np

    ref_track = getattr(getattr(output, "reference", None), attr, None)
    alt_track = getattr(getattr(output, "alternate", None), attr, None)
    if ref_track is None or alt_track is None:
        return []

    ref, ref_meta, ref_start, ref_res = normalizar_array(ref_track)
    alt, alt_meta, alt_start, alt_res = normalizar_array(alt_track)
    if ref.shape != alt.shape:
        raise ValueError(f"{modality}: REF/ALT tienen formas distintas {ref.shape} vs {alt.shape}")
    if ref_start != alt_start or ref_res != alt_res:
        raise ValueError(
            f"{modality}: AlphaGenome no devolvió perfiles REF/ALT en la misma cuadrícula "
            f"(start {ref_start}/{alt_start}, resolución {ref_res}/{alt_res}); no se calcula diferencia"
        )
    if firmas_metadata(ref_meta) != firmas_metadata(alt_meta):
        raise ValueError(f"{modality}: metadata/orden de tracks distintos entre REF y ALT")

    centers = ref_start + np.arange(ref.shape[0], dtype=np.int64) * ref_res + (ref_res // 2)
    breakpoint0 = int(variant["position_alphagenome_1based"]) - 1
    mask = (centers >= breakpoint0 - radius) & (centers <= breakpoint0 + radius)
    if not mask.any():
        return []

    ref_local = ref[mask, :]
    alt_local = alt[mask, :]
    diff = ref_local - alt_local
    abs_diff = np.abs(diff)
    rows: list[dict[str, Any]] = []
    for index in range(ref.shape[1]):
        meta = fila_metadata(ref_meta, index)
        rows.append({
            "variant_id": variant["variant_id"],
            "gen_locus": variant.get("gen_locus", ""),
            "cromosoma": variant["cromosoma"],
            "position_alphagenome_1based": variant["position_alphagenome_1based"],
            "longitud_secuencia_ancestral_bp": variant["longitud_secuencia_ancestral_bp"],
            "secuencia_ancestral_inferida": variant["secuencia_ancestral_inferida"],
            "modalidad": modality,
            "track_index": index,
            "track_name": meta.get("name", ""),
            "strand": meta.get("strand", ""),
            "ontology_curie": meta.get("ontology_curie", ""),
            "biosample_name": meta.get("biosample_name", ""),
            "biosample_type": meta.get("biosample_type", ""),
            "gtex_tissue": meta.get("gtex_tissue", ""),
            "target": meta.get("target", ""),
            "histone_mark": meta.get("histone_mark", ""),
            "assay_title": meta.get("Assay title", "") or meta.get("assay_title", ""),
            "data_source": meta.get("data_source", ""),
            "gene_id": meta.get("gene_id", ""),
            "gene_symbol": meta.get("gene_symbol", ""),
            "radio_local_bp": radius,
            "n_bins_locales": int(mask.sum()),
            "senal_humana_suma": float(np.sum(ref_local[:, index])),
            "senal_ancestral_suma": float(np.sum(alt_local[:, index])),
            "efecto_abs_L1_humano_vs_ancestral": float(np.sum(abs_diff[:, index])),
            "efecto_firmado_humano_menos_ancestral": float(np.sum(diff[:, index])),
            "efecto_firmado_ancestral_menos_humano": float(-np.sum(diff[:, index])),
            "diferencia_abs_media_por_bin": float(np.mean(abs_diff[:, index])),
            "diferencia_abs_max_bin": float(np.max(abs_diff[:, index])),
            "direccion_interpretacion": "REF=hg38 humano; ALT=secuencia ancestral insertada en la unión; signo=humano−ancestral",
        })
    return rows


def obtener_secuencia_hg38(chrom: str, start0: int, end0: int) -> str:
    query = urllib.parse.urlencode({"genome": "hg38", "chrom": chrom, "start": start0, "end": end0})
    data = obtener_json("https://api.genome.ucsc.edu/getData/sequence?" + query, timeout=120)
    seq = limpiar(data.get("dna", "")).upper()
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError(f"UCSC no devolvió una secuencia hg38 válida para {chrom}:{start0}-{end0}")
    return seq


def obtener_ventana_referencia(chrom: str, center0: int, length: int = LONGITUD_SECUENCIA) -> tuple[int, str]:
    start0 = max(0, int(center0) - length // 2)
    for _ in range(5):
        end0 = start0 + length
        seq = obtener_secuencia_hg38(chrom, start0, end0)
        if len(seq) == length:
            return start0, seq
        start0 = max(0, start0 - (length - len(seq)) - 1)
    raise ValueError(f"No se pudo obtener una ventana hg38 de {length:,} pb en {chrom} cerca de {center0}")


def cargar_kernels_cpp():
    """Carga o compila los kernels C++ cuyo código fuente vive en este mismo .py.

    La compilación se hace dentro del intérprete actual y guarda el binario en una
    caché del sistema temporal. No requiere un .cpp, setup.py ni .pyd en la carpeta
    del proyecto y no usa un proceso Python secundario.
    """
    # Busca pybind11 también en las instalaciones de Python del usuario. El
    # lanzador `py` de Windows puede seleccionar un intérprete distinto de aquel
    # donde se instaló pybind11, aunque el paquete siga presente en el equipo.
    try:
        import pybind11
    except ImportError:
        site_candidates: list[Path] = []
        local_python_root = Path.home() / "AppData" / "Local" / "Programs" / "Python"
        if local_python_root.is_dir():
            site_candidates.extend(local_python_root.glob("Python*/Lib/site-packages"))
        roaming_root = Path.home() / "AppData" / "Roaming" / "Python"
        if roaming_root.is_dir():
            site_candidates.extend(roaming_root.glob("Python*/site-packages"))
        current_tag = f"Python{sys.version_info.major}{sys.version_info.minor}"
        site_candidates.sort(key=lambda path: (0 if current_tag in str(path) else 1, str(path)))
        loaded_pybind11 = False
        for site_dir in site_candidates:
            if (site_dir / "pybind11").is_dir() and str(site_dir) not in sys.path:
                # pybind11 es principalmente Python puro. Importarlo temporalmente desde
                # otra instalación es aceptable para compilar, pero DEJAR ese site-packages
                # completo en sys.path contamina el intérprete actual: por ejemplo, Python
                # 3.11 acaba intentando cargar los .pyd de NumPy compilados para Python 3.13.
                site_text = str(site_dir)
                sys.path.insert(0, site_text)
                importlib.invalidate_caches()
                try:
                    import pybind11
                    loaded_pybind11 = True
                except ImportError:
                    loaded_pybind11 = False
                finally:
                    while site_text in sys.path:
                        sys.path.remove(site_text)
                    importlib.invalidate_caches()
                if loaded_pybind11:
                    break
        if not loaded_pybind11:
            raise RuntimeError(
                "No se pudo importar pybind11. Intérprete activo: "
                f"{sys.executable}. Rutas de site-packages revisadas: "
                + "; ".join(str(x) for x in site_candidates)
                + ". Comprueba con `py -0p` qué Python selecciona el lanzador y, "
                  "si hace falta, instala en ese intérprete con `py -m pip install pybind11`."
            )

    # Si pybind11 ya estaba disponible por PYTHONPATH desde otra versión de Python,
    # quita su site-packages externo del buscador global antes de importar cualquier
    # dependencia binaria (en particular NumPy). El submódulo setup_helpers sigue
    # siendo localizable mediante pybind11.__path__.
    pybind_package_file = Path(getattr(pybind11, "__file__", "")).resolve()
    pybind_site_dir = pybind_package_file.parent.parent if pybind_package_file.name else None
    expected_python_dir = f"Python{sys.version_info.major}{sys.version_info.minor}"
    if pybind_site_dir is not None and pybind_site_dir.name.lower() == "site-packages" and expected_python_dir.lower() not in str(pybind_site_dir).lower():
        retained_paths = []
        for item in sys.path:
            try:
                if Path(item or ".").resolve() == pybind_site_dir:
                    continue
            except OSError:
                pass
            retained_paths.append(item)
        sys.path[:] = retained_paths

    try:
        from pybind11.setup_helpers import Pybind11Extension, build_ext
        from setuptools import Distribution
    except ImportError as exc:
        raise RuntimeError(
            "Se encontró pybind11 en " + str(getattr(pybind11, "__file__", "ubicación desconocida"))
            + ", pero faltan componentes de compilación para el intérprete "
            + sys.executable + ": " + str(exc)
        ) from exc

    source_hash = hashlib.sha256(
        (CPP_KERNELS_SOURCE + sys.version + pybind11.__version__).encode("utf-8")
    ).hexdigest()[:20]
    cache_dir = (
        Path(tempfile.gettempdir())
        / "evaluar_deleciones_cpp_cache"
        / f"{sys.implementation.cache_tag or 'python'}-{source_hash}"
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    cpp_path = cache_dir / "genoma_indel_kernels.cpp"
    if not cpp_path.exists() or cpp_path.read_text(encoding="utf-8") != CPP_KERNELS_SOURCE:
        cpp_path.write_text(CPP_KERNELS_SOURCE, encoding="utf-8")

    extension = Pybind11Extension(
        "genoma_indel_kernels",
        [str(cpp_path)],
        cxx_std=17,
    )
    distribution = Distribution(
        {
            "name": "genoma_indel_kernels",
            "ext_modules": [extension],
            "cmdclass": {"build_ext": build_ext},
        }
    )
    builder = build_ext(distribution)
    builder.build_lib = str(cache_dir)
    builder.build_temp = str(cache_dir / "build")
    builder.inplace = False
    builder.ensure_finalized()
    module_path = Path(builder.get_ext_fullpath("genoma_indel_kernels"))

    if not module_path.is_file():
        try:
            builder.run()
        except Exception as exc:
            raise RuntimeError(
                "No se pudo compilar el C++ integrado. Verifica que este Python tenga "
                "pybind11 y un compilador C++ compatible (MSVC Build Tools en Windows). "
                "El código C++ ya está incluido en evaluar_deleciones.py. Detalle: "
                + f"{type(exc).__name__}: {exc}"
            ) from exc

    if not module_path.is_file():
        raise RuntimeError(
            "La compilación terminó sin producir el módulo nativo esperado: "
            + str(module_path)
        )

    spec = importlib.util.spec_from_file_location("genoma_indel_kernels", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"No se pudo crear el cargador del módulo C++: {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["genoma_indel_kernels"] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop("genoma_indel_kernels", None)
        raise

    for name in (
        "puntuacion_pwm_cpp",
        "mejor_pwm_en_intervalo_cpp",
        "reverse_complement_cpp",
    ):
        if not hasattr(module, name):
            raise RuntimeError(f"El módulo C++ integrado no expone {name}.")
    print(f"Kernels PWM C++ cargados desde el código integrado: {module_path}", flush=True)
    return module

def _bounds_pwm(matrix: list[dict[str, float]]) -> tuple[float, float]:
    return (
        sum(min(row.values()) for row in matrix),
        sum(max(row.values()) for row in matrix),
    )


def _consensus_response_model(
    motif_id: str, name: str, pattern: str, element: str,
    threshold: float, kind: str,
) -> dict[str, Any]:
    matrix: list[dict[str, float]] = []
    n_fixed = 0
    for char in pattern.upper():
        if char in "ACGT":
            matrix.append({base: (1.0 if base == char else 0.0) for base in "ACGT"})
            n_fixed += 1
        else:
            matrix.append({base: 0.0 for base in "ACGT"})
    return {
        "id": motif_id, "name": name, "width": len(pattern), "matrix": matrix,
        "min_score": 0.0, "max_score": float(n_fixed),
        "is_ar": element == "ARE", "is_esr1": False, "is_esr2": False,
        "tipo_motivo": kind, "elemento_respuesta": element,
        "familia_motivo": "AR" if element == "ARE" else "ERE_compartido_ERalpha_ERbeta",
        "umbral_modelo": threshold, "delta_min_modelo": DELTA_CONSENSO,
    }


def ampliar_motivos_respuesta(motifs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Construye paneles ARE/ERE canónicos, degenerados y de arquitectura alternativa.

    Incluye el ARE selectivo tipo repetición directa (sARE/DR3), ERE completos
    no canónicos con hasta dos discrepancias en los semisitios, medios sitios,
    arquitecturas de semisitios en tándem y espaciados alternativos. Son candidatos
    de secuencia in silico: no equivalen a unión experimental ni a exclusividad absoluta.
    """
    expanded: list[dict[str, Any]] = []
    native = {"AR": None, "ESR1": None, "ESR2": None}
    for motif in motifs:
        if motif.get("is_ar") and native["AR"] is None:
            native["AR"] = motif
        if motif.get("is_esr1") and native["ESR1"] is None:
            native["ESR1"] = motif
        if motif.get("is_esr2") and native["ESR2"] is None:
            native["ESR2"] = motif
        motif.setdefault("tipo_motivo", "PWM_JASPAR")
        motif.setdefault("umbral_modelo", UMBRAL_PWM)
        motif.setdefault("delta_min_modelo", DELTA_PWM)
        motif.setdefault("elemento_respuesta", "ARE" if motif.get("is_ar") else ("ERE" if motif.get("is_esr1") or motif.get("is_esr2") else ""))
        motif.setdefault("familia_motivo", "AR" if motif.get("is_ar") else ("ESR1" if motif.get("is_esr1") else ("ESR2" if motif.get("is_esr2") else "")))
        motif.setdefault("interpretacion_selectividad", "PWM JASPAR del receptor; preferencia in silico, no exclusividad experimental")
        expanded.append(motif)

    # Medios sitios PWM derivados de las matrices nativas de AR, ERalpha y ERbeta.
    half_pwm_models: list[dict[str, Any]] = []
    for receptor, element in (("AR", "ARE"), ("ESR1", "ERE"), ("ESR2", "ERE")):
        source = native.get(receptor)
        if not source:
            continue
        width = int(source["width"])
        if width < 12:
            continue
        for side, matrix_slice in (("izquierda", source["matrix"][:6]), ("derecha", source["matrix"][-6:])):
            matrix = [dict(row) for row in matrix_slice]
            mn, mx = _bounds_pwm(matrix)
            item = {
                "id": f"{source['id']}_HALF_{side.upper()}",
                "name": f"{receptor} medio-sitio PWM ({side}; derivado de {source['id']})",
                "width": 6, "matrix": matrix, "min_score": mn, "max_score": mx,
                "is_ar": receptor == "AR", "is_esr1": receptor == "ESR1", "is_esr2": receptor == "ESR2",
                "tipo_motivo": "PWM_HALF_SITE", "elemento_respuesta": element,
                "familia_motivo": receptor, "mitad_pwm": side,
                "umbral_modelo": 0.70, "delta_min_modelo": DELTA_PWM,
                "interpretacion_selectividad": "Medio sitio derivado de PWM; no demuestra unión aislada del receptor",
            }
            expanded.append(item)
            half_pwm_models.append(item)

    # Medios sitios de consenso. Se informan como candidatos, no como elementos completos confirmados.
    expanded.append(_consensus_response_model(
        "ARE_HALF_SITE_CONSENSUS", "ARE medio-sitio consenso", "AGAACA", "ARE", 0.80, "HALF_SITE_CONSENSUS"
    ))
    expanded[-1]["interpretacion_selectividad"] = "Medio sitio ARE candidato; no prueba exclusividad de AR"
    expanded.append(_consensus_response_model(
        "ERE_HALF_SITE_CONSENSUS", "ERE medio-sitio consenso compartido", "AGGTCA", "ERE", 0.80, "HALF_SITE_CONSENSUS"
    ))
    expanded[-1]["interpretacion_selectividad"] = "Medio sitio ERE candidato; puede requerir otros sitios/cofactores"

    # Consensus exacto frente a modelos tolerantes de secuencia: el exacto se exige perfecto;
    # el no canónico admite hasta dos discrepancias en los 12 nucleótidos de los semisitios.
    for element, left, right in (("ARE", "AGAACA", "TGTTCT"), ("ERE", "AGGTCA", "TGACCT")):
        canonical_pattern = left + "NNN" + right
        noncanonical = _consensus_response_model(
            f"{element}_NONCANONICAL_IR3_2MISMATCH",
            f"{element} completo no canónico; repetición invertida con 3 pb de separación; hasta 2 diferencias en semisitios",
            canonical_pattern, element, UMBRAL_RESPUESTA_NO_CANONICA, "NONCANONICAL_IR3",
        )
        noncanonical["interpretacion_selectividad"] = (
            "Arquitectura de repetición invertida con secuencia degenerada; candidato de motivo, no prueba funcional"
        )
        expanded.append(noncanonical)

    # sARE selectivo tipo repetición directa: AGAACA-NNN-AGAACA (y su hebra inversa,
    # cubierta por el escaneo reverse-complement del kernel C++). Otros espaciados se etiquetan
    # solo como arquitecturas alternativas exploratorias, no como sARE establecidos.
    selective_are = _consensus_response_model(
        "ARE_SELECTIVE_DR3",
        "ARE selectivo tipo sARE; repetición directa AGAACA-NNN-AGAACA",
        "AGAACA" + "NNN" + "AGAACA", "ARE", UMBRAL_ARE_SELECTIVO, "ARE_SELECTIVE_DIRECT_REPEAT",
    )
    selective_are["interpretacion_selectividad"] = (
        "Arquitectura sARE asociada a reconocimiento selectivo/preferencial de AR; no garantiza exclusión absoluta de PR u otros receptores"
    )
    expanded.append(selective_are)
    for spacer in (0, 1, 2, 4, 5, 6):
        model = _consensus_response_model(
            f"ARE_DIRECT_REPEAT_ALT_SPACER_{spacer}",
            f"ARE-like repetición directa; espaciador alternativo={spacer} pb (exploratorio)",
            "AGAACA" + ("N" * spacer) + "AGAACA", "ARE", UMBRAL_RESPUESTA_NO_CANONICA,
            "ALT_SPACER_DIRECT",
        )
        model["interpretacion_selectividad"] = "Arquitectura directa alternativa exploratoria; no asumir especificidad de AR"
        expanded.append(model)

    # ERE-like ERalpha-selective extended half-site (SFRE-like) reported as bound by ERalpha,
    # but not ERbeta in the cited binding literature. Included as an additional candidate class,
    # not renamed as a canonical ERE.
    alpha_sfre = _consensus_response_model(
        "ERE_ALPHA_SELECTIVE_SFRE_LIKE",
        "Sitio extendido SFRE-like asociado a reconocimiento de ERalpha: TCAAGGTCA",
        "TCAAGGTCA", "ERE", 0.88, "ERE_ALPHA_SELECTIVE_SFRE",
    )
    alpha_sfre["is_esr1"] = True
    alpha_sfre["familia_motivo"] = "ESR1"
    alpha_sfre["interpretacion_selectividad"] = (
        "Motivo SFRE-like asociado a unión de ERalpha y no ERbeta en estudios de unión; la coincidencia de secuencia no prueba ocupación en este locus"
    )
    expanded.append(alpha_sfre)

    # Variaciones de espaciado en la repetición invertida: se conservan explícitamente como
    # candidatos exploratorios; no se equiparan al ERE funcional de espaciado canónico de 3 pb.
    for element, left, right in (("ARE", "AGAACA", "TGTTCT"), ("ERE", "AGGTCA", "TGACCT")):
        for spacer in (0, 1, 2, 4, 5, 6):
            pattern = left + ("N" * spacer) + right
            model = _consensus_response_model(
                f"{element}_ALT_INVERTED_SPACER_{spacer}",
                f"{element}-like repetición invertida; espaciador alternativo={spacer} pb (exploratorio)",
                pattern, element, UMBRAL_RESPUESTA_NO_CANONICA, "ALT_SPACER_INVERTED",
            )
            model["interpretacion_selectividad"] = (
                "Geometría alternativa exploratoria; coincidencia de secuencia no confirma unión funcional del receptor"
            )
            expanded.append(model)

    # Repeticiones directas de semisitios ERE (incluidas en algunos contextos reguladores descritos).
    # El espaciado 3 se separa de otras geometrías, que se mantienen como exploratorias.
    for spacer in (0, 1, 2, 3, 4, 5, 6):
        kind = "ERE_TANDEM_HALF_SITE_DR3" if spacer == 3 else "ERE_TANDEM_HALF_SITE_ALT_DR"
        model = _consensus_response_model(
            f"ERE_TANDEM_HALF_DR_{spacer}",
            f"ERE-like dos medios sitios en repetición directa; espaciador={spacer} pb",
            "AGGTCA" + ("N" * spacer) + "AGGTCA", "ERE", UMBRAL_RESPUESTA_NO_CANONICA, kind,
        )
        model["interpretacion_selectividad"] = (
            "Arquitectura ERE-like de medios sitios en tándem; debe distinguirse de un ERE palindrómico canónico"
        )
        expanded.append(model)

    return expanded


def cargar_motivos_jaspar_cpp(carpeta: Path) -> list[dict[str, Any]]:
    path = carpeta / "JASPAR2026_CORE_vertebrates_nonredundant.meme.txt"
    if not path.exists() or path.stat().st_size < 1000:
        request = urllib.request.Request(URL_JASPAR_MEME, headers={"User-Agent": "indel-hormone-analysis/1.0"})
        with urllib.request.urlopen(request, timeout=180, context=crear_contexto_ssl()) as response:
            payload = response.read()
        tmp = path.with_suffix(path.suffix + ".part")
        tmp.write_bytes(payload)
        tmp.replace(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    motifs: list[dict[str, Any]] = []
    i = 0
    while i < len(lines):
        if not lines[i].startswith("MOTIF "):
            i += 1
            continue
        parts = lines[i].split(maxsplit=2)
        motif_id = parts[1] if len(parts) > 1 else "unknown"
        name = parts[2] if len(parts) > 2 else motif_id
        j = i + 1
        while j < len(lines) and not lines[j].startswith("letter-probability matrix:"):
            if lines[j].startswith("MOTIF "):
                break
            j += 1
        if j >= len(lines) or not lines[j].startswith("letter-probability matrix:"):
            i += 1
            continue
        match = re.search(r"\bw\s*=\s*(\d+)", lines[j])
        if not match:
            i = j + 1
            continue
        width = int(match.group(1))
        if width < 5 or width > 30:
            i = j + 1 + width
            continue
        probs = []
        for line in lines[j + 1:j + 1 + width]:
            try:
                row = [float(x) for x in line.split()]
            except ValueError:
                row = []
            if len(row) != 4:
                probs = []
                break
            probs.append(row)
        if len(probs) != width:
            i = j + 1
            continue
        matrix, min_score, max_score = [], 0.0, 0.0
        for row in probs:
            total = sum(row)
            vals = [(max(0.0, v) + 0.001) / (total + 0.004) for v in row]
            scores = {base: math.log2(vals[k] / 0.25) for k, base in enumerate("ACGT")}
            matrix.append(scores)
            min_score += min(scores.values())
            max_score += max(scores.values())
        normalized = re.sub(r"[^a-z0-9]+", "", (name + " " + motif_id).lower())
        base_id = motif_id.split(".")[0]
        is_ar = base_id == "MA0007" or bool(re.search(r"(^|[^a-z])ar([^a-z]|$)", name.lower())) or "androgenreceptor" in normalized
        is_esr1 = base_id == "MA0112" or "esr1" in normalized or "estrogenreceptoralpha" in normalized or "estrogenreceptoralfa" in normalized
        is_esr2 = base_id == "MA0258" or "esr2" in normalized or "estrogenreceptorbeta" in normalized
        if max_score > min_score and (is_ar or is_esr1 or is_esr2):
            element = "ARE" if is_ar else ("ERE" if is_esr1 or is_esr2 else "")
            family = "AR" if is_ar else ("ESR1" if is_esr1 else ("ESR2" if is_esr2 else ""))
            motifs.append({
                "id": motif_id, "name": name, "width": width, "matrix": matrix,
                "min_score": min_score, "max_score": max_score,
                "is_ar": is_ar, "is_esr1": is_esr1, "is_esr2": is_esr2,
                "tipo_motivo": "PWM_JASPAR", "elemento_respuesta": element,
                "familia_motivo": family, "umbral_modelo": UMBRAL_PWM,
                "delta_min_modelo": DELTA_PWM,
            })
        i = j + 1 + width
    if not motifs:
        raise RuntimeError("No se pudieron cargar matrices JASPAR PWM para AR/ESR1/ESR2.")
    native_counts = (sum(x["tipo_motivo"] == "PWM_JASPAR" and x["is_ar"] for x in motifs),
                     sum(x["tipo_motivo"] == "PWM_JASPAR" and x["is_esr1"] for x in motifs),
                     sum(x["tipo_motivo"] == "PWM_JASPAR" and x["is_esr2"] for x in motifs))
    expanded = ampliar_motivos_respuesta(motifs)
    # Conserva un indicador de modelos half-site para comparar preferencias de AR/ERalpha/ERbeta.
    half_models = [x for x in expanded if x.get("tipo_motivo") == "PWM_HALF_SITE"]
    print(
        f"Matrices PWM JASPAR: AR={native_counts[0]} | ESR1={native_counts[1]} | ESR2={native_counts[2]}; "
        f"modelos adicionales ARE/ERE={len(expanded)-len(motifs)} "
        f"(sARE tipo DR3, ARE/ERE no canónicos, medios sitios, SFRE-like ERalpha y arquitecturas alternativas).",
        flush=True,
    )
    return expanded


def consensus_matrix(kind: str) -> tuple[list[dict[str, float]], float, float]:
    if kind == "ARE":
        fixed = {0: "A", 1: "G", 2: "A", 3: "A", 4: "C", 5: "A", 9: "T", 10: "G", 11: "T", 12: "T", 13: "C", 14: "T"}
    else:
        pattern = "AGGTCANNN TGACCT".replace(" ", "")
        fixed = {i: base for i, base in enumerate(pattern) if base != "N"}
    matrix = []
    for i in range(15):
        if i in fixed:
            matrix.append({base: (1.0 if base == fixed[i] else 0.0) for base in "ACGT"})
        else:
            matrix.append({base: 0.0 for base in "ACGT"})
    return matrix, 0.0, float(len(fixed))


def mejor_par_cpp(kernel: Any, ancestral: str, human: str, center: int,
                  inserted_length: int, motif: dict[str, Any]) -> dict[str, Any] | None:
    width = int(motif["width"])
    a_lo, a_hi = max(0, center), min(len(ancestral), center + max(1, inserted_length))
    h_lo, h_hi = max(0, center), min(len(human), center + 1)
    ancestral_result = kernel.mejor_pwm_en_intervalo_cpp(
        ancestral, int(a_lo), int(a_hi), motif["matrix"], float(motif["min_score"]),
        float(motif["max_score"]), width,
    )
    human_result = kernel.mejor_pwm_en_intervalo_cpp(
        human, int(h_lo), int(h_hi), motif["matrix"], float(motif["min_score"]),
        float(motif["max_score"]), width,
    )
    if not ancestral_result or not human_result:
        return None
    ancestral_result, human_result = dict(ancestral_result), dict(human_result)
    return {
        "score_ancestral": float(ancestral_result["score"]),
        "score_humano": float(human_result["score"]),
        "delta_ancestral_menos_humano": float(ancestral_result["score"]) - float(human_result["score"]),
        "start_ancestral": int(ancestral_result["start"]),
        "start_humano": int(human_result["start"]),
        "orientacion_ancestral": ancestral_result["orientacion"],
        "orientacion_humana": human_result["orientacion"],
        "secuencia_ancestral": ancestral_result["secuencia"],
        "secuencia_humana": human_result["secuencia"],
    }

def _clasificar_preferencia(scores: dict[str, float], min_score: float = 0.70, margin: float = 0.10) -> str:
    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    if not ranked or ranked[0][1] < min_score:
        return "SIN_MOTIVO_FUERTE"
    if len(ranked) > 1 and ranked[0][1] - ranked[1][1] < margin:
        return "COMPARTIDO_O_AMBIGUO"
    return {
        "AR": "PREFERENCIA_AR_IN_SILICO",
        "ESR1": "PREFERENCIA_ERALPHA_IN_SILICO",
        "ESR2": "PREFERENCIA_ERBETA_IN_SILICO",
    }.get(ranked[0][0], "COMPARTIDO_O_AMBIGUO")


def _puntuar_secuencia_medio_sitio(
    kernel: Any, sequence: str, half_pwm_models: list[dict[str, Any]],
) -> tuple[dict[str, float], str]:
    """Compara el mismo medio sitio con modelos AR, ERalpha y ERbeta.

    Es una preferencia de motivos in silico, no una medida experimental de afinidad.
    """
    if len(sequence) != 6 or set(sequence.upper()) - set("ACGT"):
        return {"AR": 0.0, "ESR1": 0.0, "ESR2": 0.0}, "SIN_SECUENCIA_COMPARABLE"
    rev = kernel.reverse_complement_cpp(sequence.upper())
    scores = {"AR": -1.0, "ESR1": -1.0, "ESR2": -1.0}
    for model in half_pwm_models:
        receptor = model.get("familia_motivo", "")
        if receptor not in scores:
            continue
        forward = float(kernel.puntuacion_pwm_cpp(
            sequence.upper(), model["matrix"], float(model["min_score"]), float(model["max_score"])
        ))
        reverse = float(kernel.puntuacion_pwm_cpp(
            rev, model["matrix"], float(model["min_score"]), float(model["max_score"])
        ))
        scores[receptor] = max(scores[receptor], forward, reverse)
    scores = {key: max(0.0, value) for key, value in scores.items()}
    return scores, _clasificar_preferencia(scores)


def _puntuar_contexto_receptores_completo(
    kernel: Any, sequence: str, start: int, motif_width: int,
    full_pwm_models: list[dict[str, Any]],
) -> tuple[dict[str, float], str]:
    """Compara AR/ERalpha/ERbeta sobre la misma ventana local del hit PWM."""
    lo = max(0, int(start) - 2)
    hi = min(len(sequence), int(start) + int(motif_width) + 2)
    context = sequence[lo:hi].upper()
    scores = {"AR": 0.0, "ESR1": 0.0, "ESR2": 0.0}
    for model in full_pwm_models:
        receptor = model.get("familia_motivo", "")
        if receptor not in scores:
            continue
        width = int(model["width"])
        if len(context) < width:
            continue
        best = 0.0
        for offset in range(0, len(context) - width + 1):
            window = context[offset:offset + width]
            if set(window) - set("ACGT"):
                continue
            rc = kernel.reverse_complement_cpp(window)
            best = max(
                best,
                float(kernel.puntuacion_pwm_cpp(window, model["matrix"], float(model["min_score"]), float(model["max_score"]))),
                float(kernel.puntuacion_pwm_cpp(rc, model["matrix"], float(model["min_score"]), float(model["max_score"]))),
            )
        scores[receptor] = max(scores[receptor], best)
    return scores, _clasificar_preferencia(scores)


def puntuar_motivos_evento(kernel: Any, candidate: dict[str, Any], ref_seq: str,
                            region_start0: int, motifs: list[dict[str, Any]]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    pos0 = int(candidate["punto_union_interbase_hg38_0based"])
    center = pos0 - region_start0
    flank = max(FLANCO_MOTIVOS, max((m["width"] for m in motifs), default=15))
    left, right = max(0, center - flank), min(len(ref_seq), center + flank)
    human = ref_seq[left:right]
    center_local = center - left
    ancestral = human[:center_local] + candidate["secuencia_ancestral_inferida"] + human[center_local:]
    if len(ancestral) < 31:
        return {"variant_id": candidate["variant_id"], "estado_motivos": "contexto_insuficiente"}, []
    all_results: list[dict[str, Any]] = []
    family_rows: dict[str, list[dict[str, Any]]] = {"AR": [], "ESR1": [], "ESR2": []}
    half_pwm_models = [m for m in motifs if m.get("tipo_motivo") == "PWM_HALF_SITE"]
    full_pwm_models = [m for m in motifs if m.get("tipo_motivo") == "PWM_JASPAR"]
    for motif in motifs:
        result = mejor_par_cpp(kernel, ancestral, human, center_local,
                               len(candidate["secuencia_ancestral_inferida"]), motif)
        if result is None:
            continue
        family = motif.get("familia_motivo") or (
            "AR" if motif.get("is_ar") else ("ESR1" if motif.get("is_esr1") else ("ESR2" if motif.get("is_esr2") else ""))
        )
        row = {
            "variant_id": candidate["variant_id"], "gen_locus": candidate.get("gen_locus", ""),
            "cromosoma": candidate["cromosoma"], "punto_union_interbase_hg38_0based": pos0,
            "motif_id": motif["id"], "TF_motif": motif["name"], "familia": family,
            "elemento_respuesta": motif.get("elemento_respuesta", ""),
            "tipo_motivo": motif.get("tipo_motivo", "PWM_JASPAR"),
            "interpretacion_selectividad": motif.get("interpretacion_selectividad", "Candidato de secuencia in silico; no prueba función ni exclusividad"),
            "umbral_modelo": float(motif.get("umbral_modelo", UMBRAL_PWM)),
            "delta_min_modelo": float(motif.get("delta_min_modelo", DELTA_PWM)),
            "score_ancestral": float(result["score_ancestral"]),
            "score_humano_union": float(result["score_humano"]),
            "delta_ancestral_menos_humano": float(result["delta_ancestral_menos_humano"]),
            "orientacion_ancestral": result.get("orientacion_ancestral", ""),
            "orientacion_humana": result.get("orientacion_humana", ""),
            "posicion_inicio_ancestral_local": int(result.get("start_ancestral", -1)),
            "posicion_inicio_humana_local": int(result.get("start_humano", -1)),
            "secuencia_motivo_ancestral": result.get("secuencia_ancestral", ""),
            "secuencia_motivo_humana": result.get("secuencia_humana", ""),
            "preferencia_receptor_ancestral_in_silico": "",
            "preferencia_receptor_humana_in_silico": "",
            "score_AR_medio_sitio_ancestral": "", "score_ERalpha_medio_sitio_ancestral": "",
            "score_ERbeta_medio_sitio_ancestral": "", "score_AR_medio_sitio_humano": "",
            "score_ERalpha_medio_sitio_humano": "", "score_ERbeta_medio_sitio_humano": "",
            "score_AR_sitio_completo_ancestral": "", "score_ERalpha_sitio_completo_ancestral": "",
            "score_ERbeta_sitio_completo_ancestral": "", "score_AR_sitio_completo_humano": "",
            "score_ERalpha_sitio_completo_humano": "", "score_ERbeta_sitio_completo_humano": "",
        }
        if (motif.get("tipo_motivo") == "PWM_JASPAR"
                or (int(motif.get("width", 0)) >= 15 and motif.get("tipo_motivo") not in ("PWM_HALF_SITE", "HALF_SITE_CONSENSUS"))):
            full_a, pref_a = _puntuar_contexto_receptores_completo(
                kernel, ancestral, int(result.get("start_ancestral", -1)), int(motif["width"]), full_pwm_models
            )
            full_h, pref_h = _puntuar_contexto_receptores_completo(
                kernel, human, int(result.get("start_humano", -1)), int(motif["width"]), full_pwm_models
            )
            row["preferencia_receptor_ancestral_in_silico"] = pref_a
            row["preferencia_receptor_humana_in_silico"] = pref_h
            row["score_AR_sitio_completo_ancestral"] = full_a["AR"]
            row["score_ERalpha_sitio_completo_ancestral"] = full_a["ESR1"]
            row["score_ERbeta_sitio_completo_ancestral"] = full_a["ESR2"]
            row["score_AR_sitio_completo_humano"] = full_h["AR"]
            row["score_ERalpha_sitio_completo_humano"] = full_h["ESR1"]
            row["score_ERbeta_sitio_completo_humano"] = full_h["ESR2"]
        if motif.get("tipo_motivo") in ("ARE_SELECTIVE_DIRECT_REPEAT", "ALT_SPACER_DIRECT"):
            row["preferencia_receptor_ancestral_in_silico"] = (
                "ARE_SARE_LIKE_SELECTIVO_AR_CANDIDATO"
                if float(result["score_ancestral"]) >= float(motif.get("umbral_modelo", UMBRAL_ARE_SELECTIVO)) else "SIN_MOTIVO_FUERTE"
            )
            row["preferencia_receptor_humana_in_silico"] = (
                "ARE_SARE_LIKE_SELECTIVO_AR_CANDIDATO"
                if float(result["score_humano"]) >= float(motif.get("umbral_modelo", UMBRAL_ARE_SELECTIVO)) else "SIN_MOTIVO_FUERTE"
            )
        if motif.get("tipo_motivo") == "ERE_ALPHA_SELECTIVE_SFRE":
            row["preferencia_receptor_ancestral_in_silico"] = (
                "SITIO_SFRE_LIKE_ERalpha_SELECTIVO_CANDIDATO"
                if float(result["score_ancestral"]) >= float(motif.get("umbral_modelo", 0.88)) else "SIN_MOTIVO_FUERTE"
            )
            row["preferencia_receptor_humana_in_silico"] = (
                "SITIO_SFRE_LIKE_ERalpha_SELECTIVO_CANDIDATO"
                if float(result["score_humano"]) >= float(motif.get("umbral_modelo", 0.88)) else "SIN_MOTIVO_FUERTE"
            )
        if motif.get("tipo_motivo") == "PWM_HALF_SITE":
            scores_a, pref_a = _puntuar_secuencia_medio_sitio(
                kernel, row["secuencia_motivo_ancestral"], half_pwm_models
            )
            scores_h, pref_h = _puntuar_secuencia_medio_sitio(
                kernel, row["secuencia_motivo_humana"], half_pwm_models
            )
            row["preferencia_receptor_ancestral_in_silico"] = pref_a
            row["preferencia_receptor_humana_in_silico"] = pref_h
            row["score_AR_medio_sitio_ancestral"] = scores_a["AR"]
            row["score_ERalpha_medio_sitio_ancestral"] = scores_a["ESR1"]
            row["score_ERbeta_medio_sitio_ancestral"] = scores_a["ESR2"]
            row["score_AR_medio_sitio_humano"] = scores_h["AR"]
            row["score_ERalpha_medio_sitio_humano"] = scores_h["ESR1"]
            row["score_ERbeta_medio_sitio_humano"] = scores_h["ESR2"]
        all_results.append(row)
        # Las categorías originales de cada receptor se basan en su matriz JASPAR
        # nativa; los medios sitios se incorporan después mediante reglas explícitas.
        if motif.get("tipo_motivo") == "PWM_JASPAR":
            for receptor, flag in (("AR", "is_ar"), ("ESR1", "is_esr1"), ("ESR2", "is_esr2")):
                if motif.get(flag):
                    family_rows[receptor].append(row)

    consensus_results: dict[str, dict[str, Any] | None] = {}
    for kind in ("ARE", "ERE"):
        matrix, low, high = consensus_matrix(kind)
        result = mejor_par_cpp(kernel, ancestral, human, center_local,
                               len(candidate["secuencia_ancestral_inferida"]),
                               {"width": 15, "matrix": matrix, "min_score": low, "max_score": high})
        consensus_results[kind] = result
        if result:
            all_results.append({
                "variant_id": candidate["variant_id"], "gen_locus": candidate.get("gen_locus", ""),
                "cromosoma": candidate["cromosoma"], "punto_union_interbase_hg38_0based": pos0,
                "motif_id": f"{kind}_CONSENSUS", "TF_motif": f"{kind} consenso canónico",
                "familia": "AR" if kind == "ARE" else "ERE_compartido_ERalpha_ERbeta",
                "elemento_respuesta": kind, "tipo_motivo": "CONSENSO_CANONICO",
                "interpretacion_selectividad": "Consenso canónico estricto; requiere coincidencia completa; no prueba funcionalidad",
                "umbral_modelo": UMBRAL_CONSENSO, "delta_min_modelo": DELTA_CONSENSO,
                "score_ancestral": float(result["score_ancestral"]),
                "score_humano_union": float(result["score_humano"]),
                "delta_ancestral_menos_humano": float(result["score_ancestral"])-float(result["score_humano"]),
                "orientacion_ancestral": result.get("orientacion_ancestral", ""),
                "orientacion_humana": result.get("orientacion_humana", ""),
                "posicion_inicio_ancestral_local": int(result.get("start_ancestral", -1)),
                "posicion_inicio_humana_local": int(result.get("start_humano", -1)),
                "secuencia_motivo_ancestral": result.get("secuencia_ancestral", ""),
                "secuencia_motivo_humana": result.get("secuencia_humana", ""),
            })

    ar_loss_candidates = [
        row for row in family_rows["AR"]
        if row["score_ancestral"] >= UMBRAL_PWM and row["delta_ancestral_menos_humano"] >= DELTA_PWM
    ]
    esr1_gain_candidates = [
        row for row in family_rows["ESR1"]
        if row["score_humano_union"] >= UMBRAL_PWM and row["delta_ancestral_menos_humano"] <= -DELTA_PWM
    ]
    esr2_change_candidates = [
        row for row in family_rows["ESR2"]
        if max(row["score_ancestral"], row["score_humano_union"]) >= UMBRAL_PWM
        and abs(row["delta_ancestral_menos_humano"]) >= DELTA_PWM
    ]
    are_other_rows = [r for r in all_results if r.get("elemento_respuesta") == "ARE" and r.get("tipo_motivo") not in ("CONSENSO_CANONICO", "PWM_JASPAR")]
    ere_other_rows = [r for r in all_results if r.get("elemento_respuesta") == "ERE" and r.get("tipo_motivo") not in ("CONSENSO_CANONICO", "PWM_JASPAR")]
    are_alt_loss = [r for r in are_other_rows if r["score_ancestral"] >= r.get("umbral_modelo", UMBRAL_PWM)
                    and r["delta_ancestral_menos_humano"] >= r.get("delta_min_modelo", DELTA_PWM)]
    ere_alpha_selective_gain = [r for r in ere_other_rows if r.get("tipo_motivo") == "ERE_ALPHA_SELECTIVE_SFRE"
                                and r["score_humano_union"] >= r.get("umbral_modelo", UMBRAL_PWM)
                                and r["delta_ancestral_menos_humano"] <= -r.get("delta_min_modelo", DELTA_PWM)]
    ere_alt_gain = [r for r in ere_other_rows if r.get("tipo_motivo") in (
                        "NONCANONICAL_IR3", "HALF_SITE_CONSENSUS", "ALT_SPACER_INVERTED",
                        "ALT_SPACER_DIRECT", "ERE_TANDEM_HALF_SITE_DR3", "ERE_TANDEM_HALF_SITE_ALT_DR")
                    and r["score_humano_union"] >= r.get("umbral_modelo", UMBRAL_PWM)
                    and r["delta_ancestral_menos_humano"] <= -r.get("delta_min_modelo", DELTA_PWM)]
    # Los modelos half-site PWM aportan indicios adicionales de preferencia por receptor.
    half_ar_loss = [r for r in are_other_rows if r.get("familia") == "AR" and r.get("tipo_motivo") == "PWM_HALF_SITE"
                    and r["score_ancestral"] >= r.get("umbral_modelo", UMBRAL_PWM)
                    and r["delta_ancestral_menos_humano"] >= r.get("delta_min_modelo", DELTA_PWM)]
    half_esr1_gain = [r for r in ere_other_rows if r.get("familia") == "ESR1" and r.get("tipo_motivo") == "PWM_HALF_SITE"
                      and r["score_humano_union"] >= r.get("umbral_modelo", UMBRAL_PWM)
                      and r["delta_ancestral_menos_humano"] <= -r.get("delta_min_modelo", DELTA_PWM)]
    half_esr2_change = [r for r in ere_other_rows if r.get("familia") == "ESR2" and r.get("tipo_motivo") == "PWM_HALF_SITE"
                        and max(r["score_ancestral"], r["score_humano_union"]) >= r.get("umbral_modelo", UMBRAL_PWM)
                        and abs(r["delta_ancestral_menos_humano"]) >= r.get("delta_min_modelo", DELTA_PWM)]
    ar = max(ar_loss_candidates, key=lambda row: (row["delta_ancestral_menos_humano"], row["score_ancestral"]), default=None)
    esr1 = min(esr1_gain_candidates, key=lambda row: (row["delta_ancestral_menos_humano"], -row["score_humano_union"]), default=None)
    esr2 = max(esr2_change_candidates, key=lambda row: (abs(row["delta_ancestral_menos_humano"]), max(row["score_ancestral"], row["score_humano_union"])), default=None)
    are, ere = consensus_results["ARE"], consensus_results["ERE"]
    ar_loss_pwm = ar is not None or bool(half_ar_loss)
    ar_loss_are = bool(are and are["score_ancestral"] >= UMBRAL_CONSENSO and (are["score_ancestral"]-are["score_humano"]) >= DELTA_CONSENSO)
    esr1_gain_pwm = esr1 is not None or bool(half_esr1_gain) or bool(ere_alpha_selective_gain)
    ere_gain = bool(ere and ere["score_humano"] >= UMBRAL_CONSENSO and (ere["score_ancestral"]-ere["score_humano"]) <= -DELTA_CONSENSO) or bool(ere_alt_gain)
    esr2_change = esr2 is not None or bool(half_esr2_change)
    categories = []
    if ar_loss_pwm or ar_loss_are or bool(are_alt_loss):
        categories.append("AR_loss")
    if esr1_gain_pwm:
        categories.append("ESR1_gain")
    if ere_gain:
        categories.append("ERE_gain_shared")
    if esr2_change:
        categories.append("ESR2_change")

    changed_pwm_rows = [r for r in all_results
                        if r.get("preferencia_receptor_ancestral_in_silico") or r.get("preferencia_receptor_humana_in_silico")
                        if abs(r.get("delta_ancestral_menos_humano", 0.0)) >= r.get("delta_min_modelo", DELTA_PWM)]
    ancestral_preferences = sorted({r.get("preferencia_receptor_ancestral_in_silico", "") for r in changed_pwm_rows
                                    if r.get("preferencia_receptor_ancestral_in_silico") not in ("", "SIN_MOTIVO_FUERTE", "SIN_SECUENCIA_COMPARABLE")
                                    and r["score_ancestral"] >= r.get("umbral_modelo", UMBRAL_PWM)})
    human_preferences = sorted({r.get("preferencia_receptor_humana_in_silico", "") for r in changed_pwm_rows
                                if r.get("preferencia_receptor_humana_in_silico") not in ("", "SIN_MOTIVO_FUERTE", "SIN_SECUENCIA_COMPARABLE")
                                and r["score_humano_union"] >= r.get("umbral_modelo", UMBRAL_PWM)})
    summary = {
        "variant_id": candidate["variant_id"], "gen_locus": candidate.get("gen_locus", ""),
        "cromosoma": candidate["cromosoma"], "punto_union_interbase_hg38_0based": pos0,
        "longitud_ancestral_bp": len(candidate["secuencia_ancestral_inferida"]),
        "AR_PWM_ancestral": ar["score_ancestral"] if ar else "",
        "AR_PWM_humano_union": ar["score_humano_union"] if ar else "",
        "AR_PWM_delta_ancestral_menos_humano": ar["delta_ancestral_menos_humano"] if ar else "",
        "ARE_consenso_ancestral": are["score_ancestral"] if are else "",
        "ARE_consenso_humano_union": are["score_humano"] if are else "",
        "ARE_delta_ancestral_menos_humano": (are["score_ancestral"]-are["score_humano"]) if are else "",
        "ARE_modelos_no_canonicos_con_perdida": ";".join(sorted({r["motif_id"] for r in are_alt_loss if r.get("tipo_motivo") in ("NONCANONICAL_IR3", "ARE_SELECTIVE_DIRECT_REPEAT", "ALT_SPACER_DIRECT", "ALT_SPACER_INVERTED")})),
        "ARE_selectivo_sARE_DR3_con_perdida": ";".join(sorted({r["motif_id"] for r in are_alt_loss if r.get("tipo_motivo") == "ARE_SELECTIVE_DIRECT_REPEAT"})),
        "ARE_repeticiones_directas_alternativas_con_perdida": ";".join(sorted({r["motif_id"] for r in are_alt_loss if r.get("tipo_motivo") == "ALT_SPACER_DIRECT"})),
        "ARE_IR3_no_canonicos_con_perdida": ";".join(sorted({r["motif_id"] for r in are_alt_loss if r.get("tipo_motivo") == "NONCANONICAL_IR3"})),
        "ESR1_PWM_ancestral": esr1["score_ancestral"] if esr1 else "",
        "ESR1_PWM_humano_union": esr1["score_humano_union"] if esr1 else "",
        "ESR1_PWM_delta_ancestral_menos_humano": esr1["delta_ancestral_menos_humano"] if esr1 else "",
        "ESR2_PWM_ancestral": esr2["score_ancestral"] if esr2 else "",
        "ESR2_PWM_humano_union": esr2["score_humano_union"] if esr2 else "",
        "ESR2_PWM_delta_ancestral_menos_humano": esr2["delta_ancestral_menos_humano"] if esr2 else "",
        "ERE_compartido_ancestral": ere["score_ancestral"] if ere else "",
        "ERE_compartido_humano_union": ere["score_humano"] if ere else "",
        "ERE_compartido_delta_ancestral_menos_humano": (ere["score_ancestral"]-ere["score_humano"]) if ere else "",
        "ERE_modelos_no_canonicos_con_ganancia": ";".join(sorted({r["motif_id"] for r in ere_alt_gain if r.get("tipo_motivo") == "NONCANONICAL_IR3"})),
        "ERE_espaciados_invertidos_exploratorios_con_ganancia": ";".join(sorted({r["motif_id"] for r in ere_alt_gain if r.get("tipo_motivo") == "ALT_SPACER_INVERTED"})),
        "ERE_repeticiones_directas_de_medio_sitio_con_ganancia": ";".join(sorted({r["motif_id"] for r in ere_alt_gain if r.get("tipo_motivo") in ("ERE_TANDEM_HALF_SITE_DR3", "ERE_TANDEM_HALF_SITE_ALT_DR")})),
        "ERE_SFRE_like_selectivo_ERALPHA_con_ganancia": ";".join(sorted({r["motif_id"] for r in ere_alpha_selective_gain})),
        "ESR1_modelos_medio_sitio_con_ganancia": ";".join(sorted({r["motif_id"] for r in half_esr1_gain})),
        "ESR2_modelos_medio_sitio_con_cambio": ";".join(sorted({r["motif_id"] for r in half_esr2_change})),
        "preferencias_de_receptor_en_motivos_ancestral_in_silico": ";".join(ancestral_preferences),
        "preferencias_de_receptor_en_motivos_humanos_in_silico": ";".join(human_preferences),
        "nota_selectividad_receptor": "Las preferencias proceden de puntuaciones PWM comparativas y modelos de secuencia; no son prueba experimental de exclusividad o unión." ,
        "categorias_hormonales": ";".join(categories),
        "clasificacion_AR_loss": ar_loss_pwm or ar_loss_are or bool(are_alt_loss),
        "clasificacion_ESR1_gain_especifica": esr1_gain_pwm,
        "clasificacion_ERE_ganancia_compartida": ere_gain,
        "clasificacion_ESR2_cambio": esr2_change,
        "metodo_motivos": "PWM JASPAR AR/ESR1/ESR2; ARE/ERE canónicos exactos; ARE selectivo sARE-like DR3; ARE/ERE no canónicos con hasta 2 diferencias; medios sitios; ERE-like SFRE asociado a ERalpha; repeticiones directas e invertidas y espaciados alternativos exploratorios. Son candidatos in silico, no prueba de unión ni exclusividad experimental.",
    }
    return summary, all_results


def output_blocks(output: Any) -> dict[str, dict[str, Any]]:
    import numpy as np
    blocks = {}
    for modality, _, attr, _ in MODALIDADES:
        track = getattr(output, attr, None)
        if track is None or getattr(track, "values", None) is None:
            continue
        values, metadata, start, resolution = normalizar_array(track)
        blocks[modality] = {
            "values": values.astype(np.float32, copy=False), "metadata": metadata.copy(),
            "signature": tuple(firmas_metadata(metadata)), "start": start, "resolution": resolution,
        }
    return blocks


def agrupar_regiones(candidates: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    by_chrom = defaultdict(list)
    for item in candidates:
        by_chrom[item["cromosoma"]].append(item)
    regions = []
    for chrom in sorted(by_chrom):
        items = sorted(by_chrom[chrom], key=lambda x: (int(x["punto_union_interbase_hg38_0based"]), x["variant_id"]))
        current = []
        first = None
        insertion_lengths = {}
        for item in items:
            pos = int(item["punto_union_interbase_hg38_0based"])
            seq = item["secuencia_ancestral_inferida"]
            proposed_lengths = dict(insertion_lengths)
            proposed_lengths[pos] = max(proposed_lengths.get(pos, 0), len(seq))
            prior_insertions = sum(length for p, length in proposed_lengths.items() if p < pos)
            effective_span = (pos - first + prior_insertions) if current else 0
            if current and effective_span > MAX_SPAN_HAPLOTIPO:
                regions.append(current)
                current, first, insertion_lengths = [], None, {}
                proposed_lengths = {pos: len(seq)}
            if not current:
                first = pos
            current.append(item)
            insertion_lengths[pos] = max(insertion_lengths.get(pos, 0), len(seq))
        if current:
            regions.append(current)
    return regions

def haplotipo_inserciones(human_sequence: str, region_start0: int, events: list[dict[str, Any]]) -> str:
    sequence = human_sequence
    total = 0
    unique_positions = {}
    for event in events:
        pos = int(event["punto_union_interbase_hg38_0based"])
        ancestral = event["secuencia_ancestral_inferida"]
        if pos in unique_positions and unique_positions[pos] != ancestral:
            raise ValueError(f"Indels incompatibles en el mismo punto de unión {pos}")
        unique_positions[pos] = ancestral
    for pos, ancestral in sorted(unique_positions.items(), reverse=True):
        offset = pos - region_start0
        if offset < 1 or offset >= len(sequence):
            raise ValueError(f"Punto {pos} fuera de la ventana de secuencia")
        sequence = sequence[:offset] + ancestral + sequence[offset:]
        total += len(ancestral)
    if total:
        sequence = sequence[:-total]
    if len(sequence) != len(human_sequence):
        raise ValueError(f"El haplotipo editado tiene longitud {len(sequence)} y se esperaba {len(human_sequence)}")
    return sequence


def centro_en_secuencia(event: dict[str, Any], region_start0: int, inserted_events: list[dict[str, Any]]) -> int:
    pos = int(event["punto_union_interbase_hg38_0based"])
    shift = sum(len(x["secuencia_ancestral_inferida"]) for x in inserted_events
                if int(x["punto_union_interbase_hg38_0based"]) < pos)
    return pos - region_start0 + shift


def extraer_ventana(block: dict[str, Any], center_offset: int, radius: int, region_start0: int):
    values = block["values"]
    res = int(block["resolution"])
    center_coordinate_offset = region_start0 + center_offset - int(block["start"])
    center_bin = int(round((center_coordinate_offset - res / 2) / res))
    half_bins = max(1, int(math.ceil(radius / res)))
    lo, hi = center_bin - half_bins, center_bin + half_bins + 1
    if lo < 0 or hi > values.shape[0]:
        raise ValueError("Ventana local fuera del perfil AlphaGenome")
    return values[lo:hi, :]


def concatenar_ventanas(block: dict[str, Any], sites: list[dict[str, Any]], region_start0: int,
                        inserted_events: list[dict[str, Any]], radius: int):
    import numpy as np
    pieces = []
    for site in sites:
        center = centro_en_secuencia(site, region_start0, inserted_events)
        pieces.append(extraer_ventana(block, center, radius, region_start0))
    return np.concatenate(pieces, axis=0) if pieces else np.zeros((0, block["values"].shape[1]), dtype=np.float32)


def ventanas_por_sitio(block: dict[str, Any], sites: list[dict[str, Any]], region_start0: int,
                       inserted_events: list[dict[str, Any]], radius: int) -> list[Any]:
    """Devuelve vistas por sitio sin concatenar cientos de miles de bins en una sola matriz."""
    return [
        extraer_ventana(block, centro_en_secuencia(site, region_start0, inserted_events), radius, region_start0)
        for site in sites
    ]


def crear_acumulador_memmap(path: Path, baseline_windows: list[Any], n_tracks: int):
    """Crea en disco el acumulador de efectos aditivos para no reservar varios GiB de RAM."""
    import numpy as np
    offsets = []
    total_rows = 0
    for window in baseline_windows:
        offsets.append(total_rows)
        total_rows += int(window.shape[0])
    accumulator = np.memmap(path, dtype=np.float32, mode="w+", shape=(total_rows, n_tracks))
    # Inicialización por bloques; no construir una matriz temporal del tamaño completo.
    for start in range(0, total_rows, 16384):
        accumulator[start:min(start + 16384, total_rows), :] = 0.0
    accumulator.flush()
    return accumulator, offsets, total_rows


def cerrar_acumuladores_memmap(group_state: dict[str, Any], work_dir: Path) -> None:
    for state in group_state.values():
        for accumulator in state.get("sum_delta", {}).values():
            try:
                accumulator.flush()
            except Exception:
                pass
            mapped_file = getattr(accumulator, "_mmap", None)
            if mapped_file is not None:
                try:
                    mapped_file.close()
                except Exception:
                    pass
    shutil.rmtree(work_dir, ignore_errors=True)


def metadata_fields(meta: dict[str, str]) -> dict[str, str]:
    return {
        "track_name": meta.get("name", ""), "strand": meta.get("strand", ""),
        "ontology_curie": meta.get("ontology_curie", ""), "biosample_name": meta.get("biosample_name", ""),
        "biosample_type": meta.get("biosample_type", ""), "gtex_tissue": meta.get("gtex_tissue", ""),
        "target": meta.get("target", ""), "histone_mark": meta.get("histone_mark", ""),
        "assay_title": meta.get("Assay title", "") or meta.get("assay_title", ""),
        "data_source": meta.get("data_source", ""), "gene_id": meta.get("gene_id", ""),
        "gene_symbol": meta.get("gene_symbol", ""),
    }


def comparar_track_rows(joint_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Agrega efectos por track y compara total y media por variante entre categorías."""
    by_key: dict[tuple[Any, Any, Any], dict[str, Any]] = {}
    valid_states = {"completed", "single_variant_identity"}
    summed_fields = (
        "efecto_conjunto_firmado_humano_menos_ancestral", "efecto_conjunto_abs_L1",
        "suma_efectos_individuales_firmados", "suma_efectos_individuales_abs_L1",
        "no_aditividad_interaccion_firmada", "no_aditividad_interaccion_abs_L1",
    )
    for row in joint_rows:
        if row.get("estado") not in valid_states or not row.get("track_id"):
            continue
        if score_a_numero(row.get("efecto_conjunto_abs_L1")) is None:
            continue
        key = (row.get("categoria"), row.get("modalidad"), row.get("track_id"))
        if key not in by_key:
            first = dict(row)
            first["n_regiones"] = 1
            by_key[key] = first
        else:
            dst = by_key[key]
            for field in summed_fields:
                dst[field] = float(dst.get(field, 0) or 0) + float(row.get(field, 0) or 0)
            dst["n_regiones"] = int(dst.get("n_regiones", 1)) + 1
            dst["n_variantes_conjuntas"] = int(dst.get("n_variantes_conjuntas", 0)) + int(row.get("n_variantes_conjuntas", 0))
            dst["n_bins_evaluados"] = int(dst.get("n_bins_evaluados", 0) or 0) + int(row.get("n_bins_evaluados", 0) or 0)
    indexed = {(row["categoria"], row["modalidad"], row["track_id"]): row for row in by_key.values()}
    output = []
    all_tracks = {(modality, track_id) for _, modality, track_id in indexed}

    def total_value(row: dict[str, Any] | None, field: str) -> Any:
        return row.get(field, "") if row else ""

    def mean_per_variant(row: dict[str, Any] | None) -> Any:
        if not row:
            return ""
        count = int(row.get("n_variantes_conjuntas", 0) or 0)
        return float(row.get("efecto_conjunto_abs_L1", 0) or 0) / count if count else ""

    def mean_per_bin(row: dict[str, Any] | None) -> Any:
        if not row:
            return ""
        count = int(row.get("n_bins_evaluados", 0) or 0)
        return float(row.get("efecto_conjunto_abs_L1", 0) or 0) / count if count else ""

    for modality, track_id in sorted(all_tracks):
        ar = indexed.get(("AR_loss", modality, track_id))
        esr = indexed.get(("ESR1_gain", modality, track_id))
        ere = indexed.get(("ERE_gain_shared", modality, track_id))
        es2 = indexed.get(("ESR2_change", modality, track_id))
        ar_total = score_a_numero(total_value(ar, "efecto_conjunto_abs_L1"))
        esr_total = score_a_numero(total_value(esr, "efecto_conjunto_abs_L1"))
        ar_mean, esr_mean = mean_per_variant(ar), mean_per_variant(esr)
        ar_bin, esr_bin = mean_per_bin(ar), mean_per_bin(esr)
        comparable = ar_total is not None and esr_total is not None
        comparable_mean = ar_mean != "" and esr_mean != ""
        comparable_bin = ar_bin != "" and esr_bin != ""
        winner_total = "NO_COMPARABLE"
        winner_mean = "NO_COMPARABLE"
        winner_bin = "NO_COMPARABLE"
        if comparable:
            if math.isclose(ar_total, esr_total, rel_tol=1e-9, abs_tol=1e-12):
                winner_total = "EMPATE"
            else:
                winner_total = "AR_loss" if ar_total > esr_total else "ESR1_gain"
        if comparable_mean:
            if math.isclose(float(ar_mean), float(esr_mean), rel_tol=1e-9, abs_tol=1e-12):
                winner_mean = "EMPATE"
            else:
                winner_mean = "AR_loss" if float(ar_mean) > float(esr_mean) else "ESR1_gain"
        if comparable_bin:
            if math.isclose(float(ar_bin), float(esr_bin), rel_tol=1e-9, abs_tol=1e-12):
                winner_bin = "EMPATE"
            else:
                winner_bin = "AR_loss" if float(ar_bin) > float(esr_bin) else "ESR1_gain"
        output.append({
            "modalidad": modality, "track_id": track_id,
            "track_name": (ar or esr or ere or es2 or {}).get("track_name", ""),
            "ontology_curie": (ar or esr or ere or es2 or {}).get("ontology_curie", ""),
            "biosample_name": (ar or esr or ere or es2 or {}).get("biosample_name", ""),
            "impacto_conjunto_abs_AR_loss": total_value(ar, "efecto_conjunto_abs_L1"),
            "n_variantes_conjuntas_AR_loss": int(ar.get("n_variantes_conjuntas", 0) or 0) if ar else 0,
            "n_regiones_AR_loss": int(ar.get("n_regiones", 0) or 0) if ar else 0,
            "impacto_abs_medio_por_variante_AR_loss": ar_mean,
            "n_bins_evaluados_AR_loss": int(ar.get("n_bins_evaluados", 0) or 0) if ar else 0,
            "impacto_abs_medio_por_bin_AR_loss": ar_bin,
            "impacto_conjunto_abs_ESR1_gain": total_value(esr, "efecto_conjunto_abs_L1"),
            "n_variantes_conjuntas_ESR1_gain": int(esr.get("n_variantes_conjuntas", 0) or 0) if esr else 0,
            "n_regiones_ESR1_gain": int(esr.get("n_regiones", 0) or 0) if esr else 0,
            "impacto_abs_medio_por_variante_ESR1_gain": esr_mean,
            "n_bins_evaluados_ESR1_gain": int(esr.get("n_bins_evaluados", 0) or 0) if esr else 0,
            "impacto_abs_medio_por_bin_ESR1_gain": esr_bin,
            "impacto_conjunto_abs_ERE_compartido": total_value(ere, "efecto_conjunto_abs_L1"),
            "impacto_conjunto_abs_ESR2_change": total_value(es2, "efecto_conjunto_abs_L1"),
            "n_variantes_conjuntas_ESR2_change": int(es2.get("n_variantes_conjuntas", 0) or 0) if es2 else 0,
            "impacto_abs_medio_por_variante_ESR2_change": mean_per_variant(es2),
            "no_aditividad_abs_AR_loss": total_value(ar, "no_aditividad_interaccion_abs_L1"),
            "no_aditividad_abs_ESR1_gain": total_value(esr, "no_aditividad_interaccion_abs_L1"),
            "ratio_AR_loss_sobre_ESR1_gain_total": (ar_total / (esr_total + 1e-8)) if comparable else "",
            "ratio_AR_loss_sobre_ESR1_gain_por_variante": (float(ar_mean) / (float(esr_mean) + 1e-8)) if comparable_mean else "",
            "ratio_AR_loss_sobre_ESR1_gain_por_bin": (float(ar_bin) / (float(esr_bin) + 1e-8)) if comparable_bin else "",
            "ganador_por_impacto_total": winner_total,
            "ganador_por_impacto_medio_por_variante": winner_mean,
            "ganador_por_impacto_medio_por_bin": winner_bin,
            "nota": "Comparación descriptiva por el mismo track/modalidad. Se muestran total, media por variante y media por bin para reducir el sesgo por distinto número de candidatos o bins; no implica causalidad.",
        })
    return output



def comparar_efectos_individuales(individual_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compara deltas de perfiles REF/ALT por track compartido entre categorías.

    Es una comparación descriptiva de predicciones individuales, no una prueba causal ni
    una estimación de interacción entre variantes. Sólo compara el mismo track/modalidad.
    """
    grouped: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in individual_rows:
        modality = str(row.get("modalidad", ""))
        track_id = str(row.get("track_id", ""))
        if not modality or not track_id:
            continue
        abs_value = score_a_numero(row.get("magnitud_abs_L1_individual", row.get("magnitud_absoluta_score")))
        signed_value = score_a_numero(row.get("efecto_individual_firmado_humano_menos_ancestral"))
        if abs_value is None:
            continue
        categories = [x for x in str(row.get("categorias_hormonales", "")).split(";") if x]
        for category in categories:
            if category not in {"AR_loss", "ESR1_gain", "ERE_gain_shared", "ESR2_change"}:
                continue
            key = (category, modality, track_id)
            if key not in grouped:
                grouped[key] = {
                    "abs": [], "signed": [], "meta": row,
                    "variants": set(),
                }
            grouped[key]["abs"].append(abs_value)
            if signed_value is not None:
                grouped[key]["signed"].append(signed_value)
            grouped[key]["variants"].add(str(row.get("variant_id", "")))

    def summary(category: str, modality: str, track_id: str) -> dict[str, Any] | None:
        obj = grouped.get((category, modality, track_id))
        if not obj:
            return None
        values = sorted(obj["abs"])
        median = values[len(values) // 2] if len(values) % 2 else (values[len(values)//2 - 1] + values[len(values)//2]) / 2
        signed = obj["signed"]
        return {
            "n_variantes": len(obj["variants"]),
            "media_abs_score": sum(values) / len(values),
            "mediana_abs_score": median,
            "max_abs_score": max(values),
            "media_score_firmado_humano_menos_ancestral": (sum(signed) / len(signed)) if signed else "",
            "meta": obj["meta"],
        }

    output = []
    tracks = sorted({(modality, track_id) for _, modality, track_id in grouped})
    for modality, track_id in tracks:
        ar = summary("AR_loss", modality, track_id)
        esr = summary("ESR1_gain", modality, track_id)
        ere = summary("ERE_gain_shared", modality, track_id)
        es2 = summary("ESR2_change", modality, track_id)
        chosen = ar or esr or ere or es2 or {}
        meta = (chosen.get("meta") or {}) if chosen else {}
        comparable = bool(ar and esr)
        ratio = (ar["media_abs_score"] / (esr["media_abs_score"] + 1e-12)) if comparable else ""
        winner = "NO_COMPARABLE"
        if comparable:
            if math.isclose(ar["media_abs_score"], esr["media_abs_score"], rel_tol=1e-9, abs_tol=1e-12):
                winner = "EMPATE"
            else:
                winner = "AR_loss" if ar["media_abs_score"] > esr["media_abs_score"] else "ESR1_gain"
        output.append({
            "modalidad": modality,
            "track_id": track_id,
            "track_name": meta.get("track_name", ""),
            "ontology_curie": meta.get("ontology_curie", ""),
            "biosample_name": meta.get("biosample_name", ""),
            "gene_id": meta.get("gene_id", ""),
            "gene_symbol": meta.get("gene_symbol", ""),
            "n_variantes_AR_loss": ar["n_variantes"] if ar else 0,
            "media_abs_score_AR_loss": ar["media_abs_score"] if ar else "",
            "mediana_abs_score_AR_loss": ar["mediana_abs_score"] if ar else "",
            "max_abs_score_AR_loss": ar["max_abs_score"] if ar else "",
            "media_firmada_humano_menos_ancestral_AR_loss": ar["media_score_firmado_humano_menos_ancestral"] if ar else "",
            "n_variantes_ESR1_gain": esr["n_variantes"] if esr else 0,
            "media_abs_score_ESR1_gain": esr["media_abs_score"] if esr else "",
            "mediana_abs_score_ESR1_gain": esr["mediana_abs_score"] if esr else "",
            "max_abs_score_ESR1_gain": esr["max_abs_score"] if esr else "",
            "media_firmada_humano_menos_ancestral_ESR1_gain": esr["media_score_firmado_humano_menos_ancestral"] if esr else "",
            "n_variantes_ERE_gain_shared": ere["n_variantes"] if ere else 0,
            "media_abs_score_ERE_gain_shared": ere["media_abs_score"] if ere else "",
            "n_variantes_ESR2_change": es2["n_variantes"] if es2 else 0,
            "media_abs_score_ESR2_change": es2["media_abs_score"] if es2 else "",
            "ratio_media_abs_AR_loss_sobre_ESR1_gain": ratio,
            "ganador_por_media_abs_score": winner,
            "nota": "Comparación descriptiva entre scores individuales del mismo track y modalidad; no implica causalidad ni no-aditividad.",
        })
    return output


def prefiltrar_por_motivos_cpp(
    candidates: list[dict[str, Any]], out_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Escanea TODOS los candidatos con C++ antes de crear el cliente AlphaGenome.

    Devuelve candidatos con cambio de motivo hormonal y regiones con secuencia hg38
    ya descargada. Los candidatos sin una categoría objetivo no llegan a AlphaGenome.
    """
    started = time.monotonic()
    kernel = cargar_kernels_cpp()
    motifs = cargar_motivos_jaspar_cpp(Path(__file__).resolve().parent)
    all_summaries: list[dict[str, Any]] = []
    all_details: list[dict[str, Any]] = []
    screened_contexts: list[dict[str, Any]] = []
    candidate_by_id = {x["variant_id"]: x for x in candidates}
    regions = agrupar_regiones(candidates)
    print(f"\n=== FILTRO PREVIO C++: {len(candidates)} candidatos en {len(regions)} regiones ===", flush=True)
    print("AlphaGenome todavía no se ha inicializado; solo se recuperan secuencias hg38 y se puntúan motivos.", flush=True)

    for region_no, region in enumerate(regions, start=1):
        chrom = region[0]["cromosoma"]
        positions = [int(x["punto_union_interbase_hg38_0based"]) for x in region]
        insertion_by_position: dict[int, int] = {}
        for event in region:
            event_pos = int(event["punto_union_interbase_hg38_0based"])
            insertion_by_position[event_pos] = max(
                insertion_by_position.get(event_pos, 0), len(event["secuencia_ancestral_inferida"])
            )
        effective_max = max(positions) + sum(
            length for pos, length in insertion_by_position.items() if pos < max(positions)
        )
        midpoint = (min(positions) + effective_max) // 2
        try:
            start0, human_seq = obtener_ventana_referencia(chrom, midpoint)
        except Exception as exc:
            error_text = f"{type(exc).__name__}: {exc}"
            print(f"  [FILTRO {region_no}/{len(regions)}] ERROR recuperando hg38: {error_text}", flush=True)
            for event in region:
                event["categorias_hormonales"] = []
                all_summaries.append({
                    "variant_id": event["variant_id"],
                    "gen_locus": event.get("gen_locus", ""),
                    "cromosoma": event["cromosoma"],
                    "punto_union_interbase_hg38_0based": event["punto_union_interbase_hg38_0based"],
                    "estado_motivos": "ERROR_REFERENCIA",
                    "categorias_hormonales": "",
                    "estado_filtro_alphagenome": "EXCLUIDO_ERROR_DE_PRESELECCION",
                    "n_categorias_motivo": 0,
                    "motivo_filtro": error_text,
                })
            continue

        passing_region: list[dict[str, Any]] = []
        for event in region:
            try:
                summary, detail_rows = puntuar_motivos_evento(kernel, event, human_seq, start0, motifs)
                categories = [
                    value for value in limpiar(summary.get("categorias_hormonales", "")).split(";")
                    if value in CATEGORIAS_HORMONALES
                ]
                event["categorias_hormonales"] = categories
                summary["estado_filtro_alphagenome"] = (
                    "SELECCIONADO_POR_MOTIVO" if categories else "EXCLUIDO_SIN_CAMBIO_DE_MOTIVO_OBJETIVO"
                )
                summary["n_categorias_motivo"] = len(categories)
                summary["motivo_filtro"] = ";".join(categories) if categories else (
                    "Ningún cambio de motivo hormonal alcanzó los umbrales configurados"
                )
                all_summaries.append(summary)
                all_details.extend(detail_rows)
                if categories:
                    passing_region.append(event)
            except Exception as exc:
                error_text = f"{type(exc).__name__}: {exc}"
                event["categorias_hormonales"] = []
                all_summaries.append({
                    "variant_id": event["variant_id"],
                    "gen_locus": event.get("gen_locus", ""),
                    "cromosoma": event["cromosoma"],
                    "punto_union_interbase_hg38_0based": event["punto_union_interbase_hg38_0based"],
                    "estado_motivos": "ERROR_ESCANEO_PWM",
                    "categorias_hormonales": "",
                    "estado_filtro_alphagenome": "EXCLUIDO_ERROR_DE_PRESELECCION",
                    "n_categorias_motivo": 0,
                    "motivo_filtro": error_text,
                })
                print(f"    [OMITIDO] {event['variant_id']} no pudo clasificarse: {error_text}", flush=True)
        if passing_region:
            screened_contexts.append({
                "events": passing_region,
                "start0": start0,
                "human_seq": human_seq,
            })
        print(
            f"  [FILTRO {region_no}/{len(regions)}] {chrom}:{min(positions)}-{max(positions)} | "
            f"{len(region)} evaluados por C++; {len(passing_region)} pasan",
            flush=True,
        )

    summary_fields = list(dict.fromkeys(k for row in all_summaries for k in row.keys()))
    detail_fields = list(dict.fromkeys(k for row in all_details for k in row.keys()))
    escribir_csv(out_dir / "evaluacion_deleciones_motivos_hormonales.csv", all_summaries, summary_fields or ["variant_id", "estado_filtro_alphagenome"])
    escribir_csv(out_dir / "evaluacion_deleciones_motivos_detalle.csv", all_details, detail_fields or ["variant_id", "motif_id", "familia"])

    selected_ids = {x["variant_id"] for x in candidates if x.get("categorias_hormonales")}
    selected = [candidate_by_id[x["variant_id"]] for x in candidates if x["variant_id"] in selected_ids]
    failed = sum(row.get("estado_filtro_alphagenome") == "EXCLUIDO_ERROR_DE_PRESELECCION" for row in all_summaries)
    print(
        f"\n=== FILTRO C++ TERMINADO en {time.monotonic()-started:.1f} s ===\n"
        f"Candidatos evaluados: {len(candidates)} | con motivos objetivo: {len(selected)} | "
        f"sin motivos objetivo: {len(candidates)-len(selected)-failed} | con error de preselección: {failed}",
        flush=True,
    )
    print(f"Resumen del filtro: {out_dir / 'evaluacion_deleciones_motivos_hormonales.csv'}", flush=True)
    return selected, screened_contexts


def crear_salidas_alphagenome_vacias(out_dir: Path, note: str) -> None:
    """Escribe CSV con encabezados cuando ningún evento pasa el filtro previo."""
    escribir_csv(out_dir / "evaluacion_deleciones_alphagenome_tracks.csv", [], [
        "variant_id", "gen_locus", "cromosoma", "modalidad", "track_id",
        "efecto_individual_firmado_humano_menos_ancestral", "efecto_individual_abs_L1", "nota",
    ])
    escribir_csv(out_dir / "evaluacion_deleciones_impacto_conjunto_por_track.csv", [], [
        "categoria", "cromosoma", "modalidad", "track_id", "n_variantes_conjuntas",
        "efecto_conjunto_abs_L1", "no_aditividad_interaccion_abs_L1", "nota",
    ])
    escribir_csv(out_dir / "evaluacion_deleciones_comparacion_AR_vs_ESR1.csv", [], [
        "modalidad", "track_id", "impacto_conjunto_abs_AR_loss", "impacto_conjunto_abs_ESR1_gain",
        "impacto_conjunto_abs_ERE_compartido", "ganador_AR_vs_ESR1_por_impacto_abs", "nota",
    ])
    escribir_csv(out_dir / "evaluacion_deleciones_haplotipos_regiones.csv", [], [
        "region_id", "cromosoma", "n_eventos", "estado", "nota",
    ])
    escribir_csv(out_dir / "evaluacion_deleciones_alphagenome_eventos.csv", [], [
        "variant_id", "gen_locus", "cromosoma", "categorias_hormonales", "estado", "nota",
    ])
    print("No se generaron predicciones AlphaGenome: " + note, flush=True)


def instalar_timeout_rpc_alphagenome(timeout_seconds: int = 300) -> None:
    """Aplica un deadline real a las llamadas gRPC del cliente de AlphaGenome."""
    from alphagenome.protos.v1 import dna_model_service_pb2_grpc

    stub_class = dna_model_service_pb2_grpc.DnaModelServiceStub
    if getattr(stub_class, "_ar_timeout_patch_seconds", None) == timeout_seconds:
        return
    original_init = stub_class.__init__
    rpc_names = (
        "GetMetadata", "PredictSequence", "PredictInterval", "PredictVariant",
        "ScoreInterval", "ScoreVariant", "ScoreIsmVariant",
    )

    def init_with_deadline(self, channel):
        original_init(self, channel)
        for rpc_name in rpc_names:
            rpc = getattr(self, rpc_name, None)
            if not callable(rpc):
                continue

            def call_with_deadline(*args, _rpc=rpc, **kwargs):
                kwargs.setdefault("timeout", timeout_seconds)
                return _rpc(*args, **kwargs)

            setattr(self, rpc_name, call_with_deadline)

    stub_class.__init__ = init_with_deadline
    stub_class._ar_timeout_patch_seconds = timeout_seconds


def es_error_rpc_que_debe_detener_el_run(exc: BaseException) -> bool:
    """No repite masivamente una solicitud cuando cualquier RPC de AlphaGenome falla."""
    if isinstance(exc, MemoryError):
        return True
    code_method = getattr(exc, "code", None)
    if callable(code_method):
        try:
            code = code_method()
            code_name = getattr(code, "name", str(code)).upper()
            return code_name not in {"OK", "STATUSCODE.OK"}
        except Exception:
            return True
    return False


def eventos_insercion_unicos(events: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[int]]:
    """Deduplica las inserciones ancestrales por coordenada y excluye sitios conflictivos."""
    by_position: dict[int, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for event in sorted(events, key=lambda item: (int(item["punto_union_interbase_hg38_0based"]), str(item.get("variant_id", "")))):
        pos = int(event["punto_union_interbase_hg38_0based"])
        ancestral = limpiar(event.get("secuencia_ancestral_inferida", "")).upper()
        if not ancestral or set(ancestral) - BASES:
            continue
        by_position[pos][ancestral].append(event)
    clean_events: list[dict[str, Any]] = []
    conflicts: list[int] = []
    for pos in sorted(by_position):
        alternatives = by_position[pos]
        if len(alternatives) != 1:
            conflicts.append(pos)
            continue
        ancestral, records = next(iter(alternatives.items()))
        clean_events.append(records[0])
    return clean_events, conflicts


def construir_secuencia_alternativa(
    human_sequence: str, region_start0: int, events: list[dict[str, Any]],
) -> tuple[str, list[dict[str, Any]], int]:
    """Aplica todas las inserciones ancestrales en coordenadas de referencia y recorta solo el extremo derecho.

    Devuelve una secuencia del tamaño exacto aceptado por AlphaGenome. El recorte del extremo
    no se interpreta como alineamiento: la correspondencia ALT→REF se calcula después, bin a bin,
    usando la suma de longitudes insertadas antes de cada coordenada de referencia.
    """
    clean_events, conflicts = eventos_insercion_unicos(events)
    if conflicts:
        raise ValueError(
            "No se puede construir el haplotipo: hay secuencias ancestrales incompatibles en "
            + ", ".join(map(str, conflicts[:12]))
        )
    if not clean_events:
        raise ValueError("No hay eventos válidos para construir el haplotipo")

    edits: list[tuple[int, str, dict[str, Any]]] = []
    for event in clean_events:
        pos = int(event["punto_union_interbase_hg38_0based"])
        ancestral = limpiar(event.get("secuencia_ancestral_inferida", "")).upper()
        offset = pos - int(region_start0)
        if offset < 1 or offset >= len(human_sequence):
            raise ValueError(
                f"El punto de inserción {pos} queda fuera de la ventana con ancla válida "
                f"({region_start0}-{region_start0 + len(human_sequence)})."
            )
        if set(ancestral) - BASES:
            raise ValueError(f"Secuencia ancestral no válida en {pos}: {ancestral!r}")
        edits.append((offset, ancestral, event))

    edits.sort(key=lambda item: item[0])
    pieces: list[str] = []
    cursor = 0
    total_inserted = 0
    for offset, ancestral, _event in edits:
        if offset < cursor:
            raise ValueError("Los indels se solapan o no están ordenados de forma compatible")
        pieces.append(human_sequence[cursor:offset])
        pieces.append(ancestral)
        cursor = offset
        total_inserted += len(ancestral)
    pieces.append(human_sequence[cursor:])
    full_alternate = "".join(pieces)
    if len(full_alternate) != len(human_sequence) + total_inserted:
        raise AssertionError("El tamaño de la secuencia alternativa no coincide con el número de bases insertadas")
    if total_inserted >= len(human_sequence) // 4:
        raise ValueError("El haplotipo insertado es demasiado grande para conservar un contexto válido")
    alternate = full_alternate[:len(human_sequence)]
    if len(alternate) != len(human_sequence):
        raise AssertionError("No se pudo mantener la longitud de entrada requerida por AlphaGenome")
    return alternate, clean_events, total_inserted


def centros_union_eventos(centers: Any, events: list[dict[str, Any]], radius: int) -> Any:
    import numpy as np
    mask = np.zeros(len(centers), dtype=bool)
    for event in events:
        pos = int(event["punto_union_interbase_hg38_0based"])
        mask |= (centers >= pos - radius) & (centers <= pos + radius)
    return mask


def offsets_insertados_antes_de_centros(centers: Any, events: list[dict[str, Any]]) -> Any:
    import numpy as np
    ordered = sorted(
        (int(event["punto_union_interbase_hg38_0based"]), len(limpiar(event.get("secuencia_ancestral_inferida", ""))))
        for event in events
    )
    positions = np.asarray([pos for pos, _length in ordered], dtype=np.int64)
    lengths = np.asarray([length for _pos, length in ordered], dtype=np.int64)
    if not len(positions):
        return np.zeros(len(centers), dtype=np.int64)
    cumulative = np.cumsum(lengths, dtype=np.int64)
    indexes = np.searchsorted(positions, np.asarray(centers, dtype=np.int64), side="right") - 1
    output = np.zeros(len(centers), dtype=np.int64)
    valid = indexes >= 0
    output[valid] = cumulative[indexes[valid]]
    return output


def muestrear_alt_en_coordenadas_ref(
    values: Any, alt_start: int, resolution: int, target_centers: Any,
    alignment_events: list[dict[str, Any]], chunk_rows: int = 8192,
) -> Any:
    """Remapea bins ALT a centros REF usando offsets acumulados; interpola si el indel cae dentro de un bin."""
    import numpy as np
    centers = np.asarray(target_centers, dtype=np.int64)
    shifts = offsets_insertados_antes_de_centros(centers, alignment_events)
    source_indexes = (centers.astype(np.float64) + shifts - int(alt_start) - (int(resolution) // 2)) / float(resolution)
    n_bins = int(values.shape[0])
    if len(source_indexes) and (
        float(source_indexes.min()) < -1e-6 or float(source_indexes.max()) > n_bins - 1 + 1e-6
    ):
        raise ValueError(
            "El contexto ALT no cubre todas las ventanas en coordenadas hg38; "
            "la región necesita más flanco o un grupo de indels más pequeño."
        )
    source_indexes = np.clip(source_indexes, 0.0, max(0, n_bins - 1))
    result = np.empty((len(centers), values.shape[1]), dtype=np.float32)
    for start in range(0, len(centers), chunk_rows):
        end = min(start + chunk_rows, len(centers))
        idx = source_indexes[start:end]
        lo = np.floor(idx).astype(np.int64)
        hi = np.minimum(lo + 1, n_bins - 1)
        frac = (idx - lo).astype(np.float32)[:, None]
        left = values[lo, :]
        right = values[hi, :]
        result[start:end, :] = left + (right - left) * frac
    return result


def perfiles_locales_lote_con_checkpoint(
    model: Any,
    interval: Any,
    sequence: str,
    modality_specs: dict[str, tuple[str, Any, int]],
    ontology_terms: list[str],
    out_dir: Path,
    cache_dir: Path,
    checkpoint: dict[str, Any],
    task_prefix: str,
    *,
    kind: str,
    target_sites: list[dict[str, Any]] | None = None,
    target_centers_by_modality: dict[str, Any] | None = None,
    alignment_events: list[dict[str, Any]] | None = None,
    baseline_profiles: dict[str, dict[str, Any]] | None = None,
    failed_modalities_this_run: set[str] | None = None,
) -> dict[str, dict[str, Any]]:
    """Predice y guarda cada modalidad en una RPC independiente.

    Una modalidad lenta o fallida no invalida las demás. Los perfiles locales se
    cachean por modalidad; en modo delta, ALT se remapea a coordenadas hg38 usando
    offsets acumulados de los indels.
    """
    import numpy as np

    checkpoint_path = out_dir / "evaluacion_deleciones_alphagenome_checkpoint.json"
    sequence_hash = hashlib.sha256(sequence.encode("ascii", errors="strict")).hexdigest()
    alignment_events = alignment_events or []
    baseline_profiles = baseline_profiles or {}
    target_centers_by_modality = target_centers_by_modality or {}
    failed_modalities_this_run = failed_modalities_this_run if failed_modalities_this_run is not None else set()
    cache_dir.mkdir(parents=True, exist_ok=True)
    payloads: dict[str, dict[str, Any]] = {}
    pending: dict[str, dict[str, Any]] = {}

    # Primero resuelve caché e identidad de cada modalidad sin lanzar solicitudes.
    for modality, (attr, output_type, radius) in modality_specs.items():
        baseline_profile = baseline_profiles.get(modality)
        if kind == "baseline":
            if not target_sites:
                raise ValueError("La predicción base requiere al menos un sitio de interés")
            target_signature: Any = sorted(
                (int(event["punto_union_interbase_hg38_0based"]), int(radius))
                for event in target_sites
            )
        elif kind == "delta":
            if baseline_profile is None or modality not in target_centers_by_modality or not alignment_events:
                raise ValueError(f"Faltan centros, base o eventos de alineamiento para {modality}")
            target_signature = hashlib.sha256(
                np.asarray(target_centers_by_modality[modality], dtype=np.int64).tobytes()
            ).hexdigest()
        else:
            raise ValueError(f"kind de perfil desconocido: {kind}")

        align_signature = [
            (int(event["punto_union_interbase_hg38_0based"]), limpiar(event.get("secuencia_ancestral_inferida", "")).upper())
            for event in alignment_events
        ]
        baseline_sequence_sha256 = ""
        baseline_centers_sha256 = ""
        if baseline_profile:
            baseline_sequence_sha256 = str(baseline_profile.get("sequence_sha256", ""))
            if not baseline_sequence_sha256 and baseline_profile.get("sequence"):
                baseline_sequence_sha256 = hashlib.sha256(
                    baseline_profile["sequence"].encode("ascii", errors="strict")
                ).hexdigest()
            baseline_centers_sha256 = hashlib.sha256(
                np.asarray(baseline_profile["centers"], dtype=np.int64).tobytes()
            ).hexdigest()

        identity = {
            "method": ALPHAGENOME_PREDICTION_METHOD_VERSION,
            "sequence_sha256": sequence_hash,
            "baseline_sequence_sha256": baseline_sequence_sha256,
            "chromosome": str(getattr(interval, "chromosome", "")),
            "interval_start": int(getattr(interval, "start", 0)),
            "interval_end": int(getattr(interval, "end", 0)),
            "modality": modality,
            "ontology_terms": sorted(set(ontology_terms)),
            "kind": kind,
            "target_signature": target_signature,
            "alignment_events": align_signature,
            "radius": int(radius),
            "baseline_centers_sha256": baseline_centers_sha256,
        }
        digest = hashlib.sha256(json.dumps(identity, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()
        cache_path = cache_dir / f"profile_v5_{digest}.pkl.gz"
        task_id = f"{task_prefix}:{modality}"
        task = checkpoint.setdefault("tasks", {}).setdefault(task_id, {"attempts": 0, "error_history": []})
        pending[modality] = {
            "attr": attr, "output_type": output_type, "radius": int(radius), "digest": digest,
            "cache_path": cache_path, "task_id": task_id, "task": task,
            "baseline_profile": baseline_profile,
        }
        if cache_path.is_file():
            try:
                with gzip.open(cache_path, "rb") as handle:
                    cached = pickle.load(handle)
                if (
                    isinstance(cached, dict)
                    and cached.get("method") == ALPHAGENOME_PREDICTION_METHOD_VERSION
                    and cached.get("kind") == kind
                    and cached.get("modality") == modality
                    and cached.get("identity_sha256") == digest
                ):
                    payloads[modality] = cached
                    task.update({
                        "status": "completed", "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                        "completed_at": task.get("completed_at") or time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                        "cache_file": str(cache_path.relative_to(out_dir)), "last_error": "", "blocked_by": [],
                    })
                    print(f"      [CACHÉ PERFIL] {modality}: {task_id}", flush=True)
                    continue
            except Exception as exc:
                task.setdefault("error_history", []).append(f"Caché ilegible ({type(exc).__name__}): {exc}")

    missing = [name for name in modality_specs if name not in payloads]
    if not missing:
        guardar_json_atomico(checkpoint_path, checkpoint)
        return payloads

    # Cada modalidad se solicita y se persiste por separado. Si una RPC agota el
    # plazo, las siguientes modalidades aún se intentan y las completadas quedan cacheadas.
    for modality in missing:
        item = pending[modality]
        attr = item["attr"]
        if modality in failed_modalities_this_run:
            task = item["task"]
            task.update({
                "status": "not_attempted_due_to_modality_timeout",
                "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "last_error": f"Se omite en esta ejecución porque una RPC de {modality} falló previamente; se reintentará al iniciar otra ejecución.",
                "blocked_by": [f"modality_timeout:{modality}"],
            })
            guardar_json_atomico(checkpoint_path, checkpoint)
            print(f"      [OMITIDA {modality}] ya agotó el timeout en esta ejecución; no se repetirá por cada variante", flush=True)
            continue

        radius = item["radius"]
        cache_path = item["cache_path"]
        task = item["task"]
        task["attempts"] = int(task.get("attempts", 0)) + 1
        task.update({
            "status": "running", "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "cache_file": str(cache_path.relative_to(out_dir)), "last_error": "", "blocked_by": [],
            "modality": modality, "kind": kind,
        })
        guardar_json_atomico(checkpoint_path, checkpoint)
        print(
            f"      [PREDICCIÓN {kind.upper()}] modalidad={modality}; "
            f"secuencia={sequence_hash[:12]}; tarea={item['task_id']}",
            flush=True,
        )
        try:
            output = model.predict_sequence(
                sequence=sequence,
                interval=interval,
                requested_outputs=[item["output_type"]],
                ontology_terms=ontology_terms,
            )
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            now = time.strftime("%Y-%m-%dT%H:%M:%S%z")
            task.update({"status": "failed", "updated_at": now, "last_error": message})
            task.setdefault("error_history", []).append(message)
            checkpoint.setdefault("failure_history", []).append({
                "timestamp": now, "task_id": item["task_id"],
                "attempt": task["attempts"], "error": message,
            })
            if es_error_rpc_que_debe_detener_el_run(exc):
                failed_modalities_this_run.add(modality)
            guardar_json_atomico(checkpoint_path, checkpoint)
            print(
                f"      [FALLO SOLICITUD {modality}] {message}; "
                f"se intentarán las demás modalidades; no se repetirá {modality} en este run",
                flush=True,
            )
            continue

        now = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        try:
            track = getattr(output, attr, None)
            if track is None or getattr(track, "values", None) is None:
                raise ValueError(f"AlphaGenome no devolvió tracks para {modality} con las ontologías solicitadas")
            values, metadata, grid_start, resolution = normalizar_array(track)

            # Aplicar el filtro explícitamente antes de calcular perfiles y guardarlos.
            # No asumimos que predict_sequence haya limitado los tracks a ontology_terms.
            if ontology_terms:
                if "ontology_curie" not in getattr(metadata, "columns", []):
                    raise ValueError(
                        f"{modality}: la metadata de AlphaGenome no contiene ontology_curie; "
                        "se cancela para no guardar tracks sin filtrar."
                    )
                wanted_ontologies = {str(term).strip() for term in ontology_terms if str(term).strip()}
                ontology_mask = metadata["ontology_curie"].astype(str).str.strip().isin(wanted_ontologies).to_numpy()
                if not ontology_mask.any():
                    observed = sorted(set(metadata["ontology_curie"].astype(str).str.strip()))[:20]
                    raise ValueError(
                        f"{modality}: 0 tracks coinciden con ontology_terms; "
                        f"solicitadas={sorted(wanted_ontologies)}, ejemplos_recibidos={observed}. "
                        "No se guardará una caché sin filtrar."
                    )
                values = values[:, ontology_mask]
                metadata = metadata.loc[ontology_mask].reset_index(drop=True)

            if len(metadata) == 0 or values.shape[1] == 0:
                raise ValueError(
                    f"AlphaGenome devolvió cero tracks para {modality} tras filtrar las ontologías: "
                    + ",".join(ontology_terms)
                )
            signature = tuple(firmas_metadata(metadata))
            full_centers = grid_start + np.arange(values.shape[0], dtype=np.int64) * resolution + (resolution // 2)

            if kind == "baseline":
                mask = centros_union_eventos(full_centers, target_sites or [], radius)
                centers = full_centers[mask]
                local_values = values[mask, :].astype(np.float32, copy=True)
                if not len(centers):
                    raise ValueError(f"No hay bins de {modality} dentro de las ventanas de los candidatos")
                payload = {
                    "method": ALPHAGENOME_PREDICTION_METHOD_VERSION,
                    "identity_sha256": item["digest"], "kind": "baseline", "modality": modality,
                    "centers": centers, "values": local_values,
                    "metadata": metadata.copy(), "signature": signature,
                    "grid_start": grid_start, "resolution": resolution, "radius": radius,
                    "sequence_sha256": sequence_hash,
                }
            else:
                baseline_profile = item["baseline_profile"]
                if baseline_profile is None:
                    raise ValueError(f"Falta perfil REF para {modality}")
                if grid_start != int(baseline_profile["grid_start"]) or resolution != int(baseline_profile["resolution"]):
                    raise ValueError(
                        f"La cuadrícula de {modality} cambió entre REF y ALT: "
                        f"start {grid_start}/{baseline_profile['grid_start']}, "
                        f"resolution {resolution}/{baseline_profile['resolution']}"
                    )
                if signature != tuple(baseline_profile["signature"]):
                    raise ValueError(f"La metadata/orden de tracks de {modality} no coincide entre REF y ALT")
                centers = np.asarray(target_centers_by_modality[modality], dtype=np.int64)
                baseline_centers = np.asarray(baseline_profile["centers"], dtype=np.int64)
                baseline_values_all = np.asarray(baseline_profile["values"], dtype=np.float32)
                indexes = np.searchsorted(baseline_centers, centers)
                if len(indexes) and (
                    np.any(indexes >= len(baseline_centers))
                    or not np.array_equal(baseline_centers[indexes], centers)
                ):
                    raise ValueError("Los centros solicitados no son un subconjunto exacto de la cuadrícula REF")
                baseline_values = baseline_values_all[indexes, :]
                if baseline_values.shape != (len(centers), len(metadata)):
                    raise ValueError("La forma de la base REF seleccionada no coincide con los centros y tracks")
                aligned_alt = muestrear_alt_en_coordenadas_ref(
                    values, grid_start, resolution, centers, alignment_events or []
                )
                if aligned_alt.shape != baseline_values.shape:
                    raise ValueError("La forma ALT alineada no coincide con la referencia")
                np.subtract(baseline_values, aligned_alt, out=aligned_alt)
                payload = {
                    "method": ALPHAGENOME_PREDICTION_METHOD_VERSION,
                    "identity_sha256": item["digest"], "kind": "delta", "modality": modality,
                    "centers": centers.copy(), "delta": aligned_alt,
                    "metadata": metadata.copy(), "signature": signature,
                    "grid_start": grid_start, "resolution": resolution, "radius": radius,
                    "sequence_sha256": sequence_hash,
                }

            fd, temp_name = tempfile.mkstemp(prefix=cache_path.name + ".", suffix=".tmp", dir=str(cache_path.parent))
            try:
                with os.fdopen(fd, "wb") as raw_handle:
                    with gzip.GzipFile(fileobj=raw_handle, mode="wb", compresslevel=1) as handle:
                        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    raw_handle.flush()
                    os.fsync(raw_handle.fileno())
                os.replace(temp_name, cache_path)
            finally:
                try:
                    Path(temp_name).unlink(missing_ok=True)
                except OSError:
                    pass
            payloads[modality] = payload
            task.update({
                "status": "completed", "updated_at": now, "completed_at": now,
                "cache_file": str(cache_path.relative_to(out_dir)), "last_error": "", "blocked_by": [],
                "modality": modality, "kind": kind,
            })
            print(
                f"      [PERFIL GUARDADO] {modality}: {len(payload['centers'])} bins, "
                f"{len(payload['metadata'])} tracks",
                flush=True,
            )
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            task.update({"status": "failed", "updated_at": now, "last_error": message})
            task.setdefault("error_history", []).append(message)
            checkpoint.setdefault("failure_history", []).append({
                "timestamp": now, "task_id": item["task_id"], "attempt": task["attempts"], "error": message,
            })
            print(f"      [FALLO PERFIL {modality}] {message}; se intentarán las demás modalidades", flush=True)
        finally:
            try:
                del output
            except UnboundLocalError:
                pass
            try:
                del track, values
            except UnboundLocalError:
                pass
            try:
                import gc
                gc.collect()
            except Exception:
                pass
            guardar_json_atomico(checkpoint_path, checkpoint)

    return payloads

def track_factor_directo(metadata_row: dict[str, Any]) -> str:
    haystack = " ".join(limpiar(metadata_row.get(key, "")) for key in (
        "name", "target", "Assay title", "assay_title", "gene_symbol", "data_source"
    )).lower()
    if re.search(r"androgen receptor|(?<![a-z0-9])ar(?![a-z0-9])", haystack):
        return "AR"
    if re.search(r"\besr1\b|estrogen receptor alpha|oestrogen receptor alpha|er[- ]?alpha|erα", haystack):
        return "ESR1"
    if re.search(r"\besr2\b|estrogen receptor beta|oestrogen receptor beta|er[- ]?beta|erβ", haystack):
        return "ESR2"
    return ""


def calcular_metricas_track(
    modality: str, track_idx: int, track_id: str, meta: dict[str, str], baseline_values: Any,
    joint_delta: Any, sum_single_delta: Any, sum_single_abs: float,
    category: str, event_count: int, centers: Any, radius: int,
) -> dict[str, Any]:
    """Calcula magnitudes e interacción para un único track, sobre los mismos bins REF."""
    import numpy as np
    joint_delta = np.asarray(joint_delta, dtype=np.float32)
    sum_single_delta = np.asarray(sum_single_delta, dtype=np.float32)
    base = np.asarray(baseline_values, dtype=np.float32)
    if joint_delta.shape != sum_single_delta.shape or joint_delta.shape != base.shape:
        raise ValueError("Las matrices de efecto conjunto, suma individual y referencia no coinciden")
    if track_idx < 0 or track_idx >= joint_delta.shape[1]:
        raise IndexError(f"Índice de track fuera de rango: {track_idx}/{joint_delta.shape[1]}")
    joint_track = joint_delta[:, track_idx]
    single_track = sum_single_delta[:, track_idx]
    joint_signed = float(np.sum(joint_track, dtype=np.float64))
    joint_abs = float(np.sum(np.abs(joint_track), dtype=np.float64))
    singles_signed = float(np.sum(single_track, dtype=np.float64))
    interaction_signed = joint_signed - singles_signed
    interaction_abs = 0.0
    for start in range(0, len(joint_track), 8192):
        end = min(start + 8192, len(joint_track))
        residual_track = joint_track[start:end] - single_track[start:end]
        interaction_abs += float(np.sum(np.abs(residual_track), dtype=np.float64))
    baseline_track = base[:, track_idx]
    baseline_abs = float(np.sum(np.abs(baseline_track), dtype=np.float64))
    ancestral_track = baseline_track - joint_track
    ancestral_abs = float(np.sum(np.abs(ancestral_track), dtype=np.float64))
    n_bins = len(centers)
    return {
        "efecto_conjunto_firmado_humano_menos_ancestral": joint_signed,
        "efecto_conjunto_abs_L1": joint_abs,
        "efecto_conjunto_abs_medio_por_bin": joint_abs / n_bins if n_bins else 0.0,
        "efecto_conjunto_relativo_pct_vs_senal_humana": 100.0 * joint_signed / (baseline_abs + 1e-8),
        "suma_efectos_individuales_firmados": singles_signed,
        "suma_efectos_individuales_abs_L1": float(sum_single_abs),
        "no_aditividad_interaccion_firmada": interaction_signed,
        "no_aditividad_interaccion_abs_L1": interaction_abs,
        "interaccion_abs_sobre_suma_individual_abs": interaction_abs / sum_single_abs if sum_single_abs else "",
        "efecto_conjunto_abs_por_variante": joint_abs / event_count if event_count else "",
        "suma_individual_abs_por_variante": sum_single_abs / event_count if event_count else "",
        "interaccion_abs_por_variante": interaction_abs / event_count if event_count else "",
        "senal_humana_abs_integrada": baseline_abs,
        "senal_ancestral_abs_integrada": ancestral_abs,
        "track_idx": track_idx,
        "track_id": track_id,
    }


def ejecutar_alphagenome(
    candidates: list[dict[str, Any]], region_contexts: list[dict[str, Any]],
    out_dir: Path, key_file: Path, ontology_terms: list[str], reset_checkpoint: bool = False,
    modality_names: list[str] | None = None,
) -> None:
    import numpy as np
    try:
        from alphagenome.models import dna_client
    except ImportError:
        try:
            from alphagenome.models.v1 import dna_client
        except ImportError as exc:
            raise RuntimeError("Falta AlphaGenome en este entorno: " + str(exc)) from exc
    try:
        from alphagenome.data import genome
    except ImportError as exc:
        raise RuntimeError("No se pudo importar alphagenome.data.genome: " + str(exc)) from exc

    checkpoint_path = out_dir / "evaluacion_deleciones_alphagenome_checkpoint.json"
    cache_dir = out_dir / "alphagenome_cache"
    work_dir = out_dir / "alphagenome_work"
    if reset_checkpoint:
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
        if checkpoint_path.exists():
            checkpoint_path.unlink()
        if work_dir.exists():
            shutil.rmtree(work_dir)
        for name in (
            "evaluacion_deleciones_alphagenome_tareas.csv",
            "evaluacion_deleciones_alphagenome_historial_fallos.csv",
        ):
            path = out_dir / name
            if path.exists():
                path.unlink()

    checkpoint = cargar_checkpoint(checkpoint_path)
    checkpoint.setdefault("tasks", {})
    checkpoint.setdefault("failure_history", [])
    previous_method = checkpoint.get("prediction_method_version")
    if previous_method != ALPHAGENOME_PREDICTION_METHOD_VERSION:
        archivar_resultados_metodo_anterior(out_dir)
        checkpoint["prediction_method_version"] = ALPHAGENOME_PREDICTION_METHOD_VERSION
        checkpoint["joint_analysis_status"] = "joint_haplotypes_enabled_custom_coordinate_mapped_multi_indel_v5_local"
        for old_task_id, old_entry in checkpoint["tasks"].items():
            old_entry["status"] = "superseded_previous_method"
            old_entry["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
            old_entry["last_error"] = (
                f"Tarea conservada para historial; no se reutiliza con el método {ALPHAGENOME_PREDICTION_METHOD_VERSION}."
            )
        print(
            "Método actualizado: perfiles REF y ALT por secuencia, una solicitud RPC independiente por modalidad faltante, "
            "caché local por modalidad/ventana, indels ancestrales individuales y haplotipos conjuntos "
            "por categoría/región, y no-aditividad calculada sobre los mismos bins hg38. "
            "Los resultados antiguos se archivan; sus cachés se conservan separadas por identidad y método.",
            flush=True,
        )
    guardar_estado_checkpoint(out_dir, checkpoint)

    key = obtener_api_key(key_file)
    rpc_timeout_seconds = max(30, int(os.environ.get("ALPHAGENOME_RPC_TIMEOUT_SECONDS", "300")))
    instalar_timeout_rpc_alphagenome(rpc_timeout_seconds)
    model = dna_client.create(key, timeout=30)

    available_names = [entry[0] for entry in MODALIDADES]
    modality_priority = ["ATAC", "DNASE", "CHIP_TF", "CHIP_HISTONE", "RNA_SEQ", "CAGE", "PROCAP"]
    if modality_names is None:
        chosen_names = modality_priority.copy()
    else:
        chosen = {str(name).strip().upper() for name in modality_names if str(name).strip()}
        unknown = sorted(chosen - set(available_names))
        if unknown:
            raise ValueError(
                "Modalidades AlphaGenome no reconocidas: " + ", ".join(unknown)
                + ". Valores válidos: " + ", ".join(available_names)
            )
        chosen_names = [name for name in modality_priority if name in chosen]
        if not chosen_names:
            raise ValueError("Debes solicitar al menos una modalidad AlphaGenome con --modalities")
    all_specs = {
        name: (attr, getattr(dna_client.OutputType, name), int(radius))
        for name, _enum, attr, radius in MODALIDADES
    }
    modality_specs = {name: all_specs[name] for name in chosen_names}

    individual_rows: list[dict[str, Any]] = []
    joint_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    started_all = time.monotonic()
    durations: list[float] = []
    failed_modalities_this_run: set[str] = set()

    print(
        f"Candidatos seleccionados por C++ para AlphaGenome: {len(candidates)} | regiones: {len(region_contexts)}",
        flush=True,
    )
    print(
        "Flujo integrado: una RPC independiente por modalidad y secuencia; calcula efectos individuales, "
        "haplotipos multindel por categoría y no-aditividad. Si una modalidad agota el timeout, se omite el resto "
        "de esa modalidad durante este run y continúan las demás; las modalidades completas se guardan en caché. "
        "Los perfiles individuales se cachean por ventana local y se recargan uno por uno para no mantener todas las matrices en RAM.",
        flush=True,
    )
    print(
        "Alineamiento multindel usado en esta versión: predicción de la secuencia ALT reconstruida y remapeo de "
        "los centros de señal a hg38 mediante offsets acumulados, con interpolación en límites fraccionarios de bin. "
        "No se afirma que la API aplique indel stitching oficial a múltiples indels internos.",
        flush=True,
    )
    print(f"Timeout RPC por solicitud: {rpc_timeout_seconds} s", flush=True)
    print("Modalidades solicitadas: " + ",".join(chosen_names), flush=True)
    print("Ontologías solicitadas: " + ",".join(ontology_terms), flush=True)
    print(f"Checkpoint: {checkpoint_path}", flush=True)
    print(f"Caché: {cache_dir}", flush=True)

    def guardar_tablas_parciales() -> None:
        individual_fields = list(dict.fromkeys(k for row in individual_rows for k in row.keys())) or [
            "region_id", "variant_id", "categoria", "modalidad", "track_id",
            "efecto_individual_firmado_humano_menos_ancestral", "magnitud_abs_L1_individual",
        ]
        joint_fields = list(dict.fromkeys(k for row in joint_rows for k in row.keys())) or [
            "region_id", "categoria", "modalidad", "track_id", "estado", "nota_metodologica",
        ]
        region_fields = list(dict.fromkeys(k for row in region_rows for k in row.keys())) or [
            "region_id", "cromosoma", "n_eventos", "estado", "nota_metodologica",
        ]
        escribir_csv(out_dir / "evaluacion_deleciones_alphagenome_tracks.csv", individual_rows, individual_fields)
        escribir_csv(out_dir / "evaluacion_deleciones_impacto_conjunto_por_track.csv", joint_rows, joint_fields)
        escribir_csv(out_dir / "evaluacion_deleciones_haplotipos_regiones.csv", region_rows, region_fields)
        individual_comparison = comparar_efectos_individuales(individual_rows)
        escribir_csv(
            out_dir / "evaluacion_deleciones_comparacion_individual_AR_vs_ESR1.csv", individual_comparison,
            list(dict.fromkeys(k for row in individual_comparison for k in row.keys())) or [
                "modalidad", "track_id", "ganador_por_media_abs_L1",
            ],
        )
        joint_comparison = comparar_track_rows(joint_rows)
        escribir_csv(
            out_dir / "evaluacion_deleciones_comparacion_AR_vs_ESR1.csv", joint_comparison,
            list(dict.fromkeys(k for row in joint_comparison for k in row.keys())) or [
                "modalidad", "track_id", "ganador_AR_vs_ESR1_por_impacto_abs",
            ],
        )
        guardar_estado_checkpoint(out_dir, checkpoint)

    failed_regions = 0
    stop_run = False
    for region_no, context in enumerate(region_contexts, start=1):
        region = context.get("events", [])
        if not region:
            continue
        started_region = time.monotonic()
        chrom = normalizar_cromosoma(region[0]["cromosoma"])
        start0 = int(context["start0"])
        human_seq = str(context["human_seq"]).upper()
        if len(human_seq) not in {16384, 131072, 524288, 1048576}:
            raise ValueError(f"Longitud de contexto AlphaGenome no admitida: {len(human_seq)} bp")
        if set(human_seq) - set("ACGTN"):
            raise ValueError("La secuencia hg38 contiene símbolos no admitidos para AlphaGenome")
        interval = genome.Interval(chromosome=chrom, start=start0, end=start0 + len(human_seq))
        positions = [int(event["punto_union_interbase_hg38_0based"]) for event in region]
        region_identity = (
            f"{chrom}_{min(positions)}_{max(positions)}_"
            f"{hashlib.sha1('|'.join(sorted(str(x['variant_id']) for x in region)).encode()).hexdigest()[:8]}"
        )

        category_events: dict[str, list[dict[str, Any]]] = {cat: [] for cat in CATEGORIAS_HORMONALES}
        category_conflicts: dict[str, list[int]] = {}
        for category in CATEGORIAS_HORMONALES:
            raw = [event for event in region if category in event.get("categorias_hormonales", [])]
            category_events[category], category_conflicts[category] = eventos_insercion_unicos(raw)

        print(
            f"\n[REGIÓN {region_no}/{len(region_contexts)}] {chrom}:{min(positions)}-{max(positions)} | "
            f"eventos={len(region)} | categorías={sum(bool(x) for x in category_events.values())}",
            flush=True,
        )
        region_errors: list[str] = []
        for category in CATEGORIAS_HORMONALES:
            if category_conflicts.get(category):
                print(
                    f"  [AVISO {category}] se excluyen posiciones con alelos ancestrales incompatibles: "
                    f"{category_conflicts[category][:12]}", flush=True,
                )

        if not any(category_events.values()):
            region_rows.append({
                "region_id": region_no, "cromosoma": chrom, "n_eventos": len(region),
                "estado": "sin_eventos_compatibles",
                "detalle_error": "No hay secuencias ancestrales compatibles dentro de las categorías hormonales.",
            })
            guardar_tablas_parciales()
            continue

        try:
            baseline_profiles = perfiles_locales_lote_con_checkpoint(
                model, interval, human_seq, modality_specs, ontology_terms,
                out_dir, cache_dir, checkpoint,
                f"v5:{region_identity}:baseline", kind="baseline", target_sites=region,
                failed_modalities_this_run=failed_modalities_this_run,
            )
        except Exception as exc:
            message = f"baseline: {type(exc).__name__}: {exc}"
            region_errors.append(message)
            print(f"  [FALLO BASELINE] {message}", flush=True)
            baseline_profiles = {}
            if es_error_rpc_que_debe_detener_el_run(exc):
                stop_run = True

        missing_baseline_modalities = sorted(set(modality_specs) - set(baseline_profiles))
        if missing_baseline_modalities:
            region_errors.append("baseline sin modalidades: " + ",".join(missing_baseline_modalities))
            print(
                "  [AVISO BASELINE] modalidades no disponibles en esta región: "
                + ",".join(missing_baseline_modalities), flush=True,
            )

        if not baseline_profiles:
            state = "interrumpido_rpc_reanudable" if stop_run else "fallo_baseline"
            region_rows.append({
                "region_id": region_no, "cromosoma": chrom, "n_eventos": len(region),
                "inicio_interbase_0based": min(positions), "fin_interbase_0based": max(positions),
                "estado": state, "detalle_error": " | ".join(region_errors),
            })
            guardar_tablas_parciales()
            if stop_run:
                failed_regions += 1
                break
            continue

        all_centers_by_modality = {
            modality: np.asarray(profile["centers"], dtype=np.int64)
            for modality, profile in baseline_profiles.items()
        }
        all_events_by_id: dict[str, dict[str, Any]] = {}
        event_categories: dict[str, set[str]] = defaultdict(set)
        for category, sites in category_events.items():
            for event in sites:
                variant_id = str(event.get("variant_id", ""))
                if variant_id:
                    all_events_by_id.setdefault(variant_id, event)
                    event_categories[variant_id].add(category)
        ordered_events = sorted(
            all_events_by_id.values(),
            key=lambda event: (int(event["punto_union_interbase_hg38_0based"]), str(event.get("variant_id", ""))),
        )
        single_modalities_by_variant: dict[str, set[str]] = {}

        print(
            f"  [REF GUARDADO] modalidades={','.join(baseline_profiles)}; "
            f"bins locales por modalidad={ {k: len(v['centers']) for k, v in baseline_profiles.items()} }",
            flush=True,
        )

        # Un ALT individual por variante y región; cada modalidad faltante usa una RPC independiente.
        for event_idx, event in enumerate(ordered_events, start=1):
            variant_id = str(event.get("variant_id", ""))
            try:
                position = int(event["punto_union_interbase_hg38_0based"])
                alternate_sequence, single_events, _total_inserted = construir_secuencia_alternativa(
                    human_seq, start0, [event]
                )
                single_centers_by_modality: dict[str, Any] = {}
                for modality, baseline in baseline_profiles.items():
                    base_centers = np.asarray(baseline["centers"], dtype=np.int64)
                    radius = int(modality_specs[modality][2])
                    local_mask = np.abs(base_centers - position) <= radius
                    if not local_mask.any():
                        region_errors.append(f"single {variant_id}/{modality}: no hay bins en la ventana local")
                        continue
                    single_centers_by_modality[modality] = base_centers[local_mask]
                delta_profiles = perfiles_locales_lote_con_checkpoint(
                    model, interval, alternate_sequence,
                    {name: modality_specs[name] for name in single_centers_by_modality}, ontology_terms,
                    out_dir, cache_dir, checkpoint,
                    f"v5:{region_identity}:single:{variant_id}", kind="delta",
                    target_centers_by_modality=single_centers_by_modality,
                    alignment_events=single_events,
                    baseline_profiles={name: baseline_profiles[name] for name in single_centers_by_modality},
                    failed_modalities_this_run=failed_modalities_this_run,
                )
                single_modalities_by_variant[variant_id] = set(delta_profiles)
                missing_modes = sorted(set(baseline_profiles) - set(delta_profiles))
                if missing_modes:
                    region_errors.append(f"single {variant_id}: faltan modalidades {','.join(missing_modes)}")

                for category in sorted(event_categories.get(variant_id, set())):
                    for modality, delta_payload in delta_profiles.items():
                        baseline = baseline_profiles[modality]
                        centers = np.asarray(delta_payload["centers"], dtype=np.int64)
                        base_centers = np.asarray(baseline["centers"], dtype=np.int64)
                        base_indexes = np.searchsorted(base_centers, centers)
                        if len(base_indexes) and (
                            np.any(base_indexes >= len(base_centers))
                            or not np.array_equal(base_centers[base_indexes], centers)
                        ):
                            raise ValueError(f"Los centros locales ALT/REF no coinciden para {variant_id}/{modality}")
                        if not len(centers):
                            region_errors.append(f"single {variant_id}/{modality}: no hay bins en la ventana local")
                            continue
                        base_local = np.asarray(baseline["values"][base_indexes, :], dtype=np.float32)
                        delta_local = np.asarray(delta_payload["delta"], dtype=np.float32)
                        if delta_local.shape != base_local.shape:
                            raise ValueError(f"Las formas del delta y REF no coinciden para {variant_id}/{modality}")
                        ancestral_local = base_local - delta_local
                        for track_idx in range(delta_local.shape[1]):
                            meta = fila_metadata(baseline["metadata"], track_idx)
                            signature = baseline["signature"][track_idx]
                            track_id = json.dumps([modality, list(signature)], ensure_ascii=False, separators=(",", ":"))
                            signed = float(np.sum(delta_local[:, track_idx], dtype=np.float64))
                            abs_l1 = float(np.sum(np.abs(delta_local[:, track_idx]), dtype=np.float64))
                            direct_factor = track_factor_directo(meta)
                            direct_direction = ""
                            if category == "AR_loss" and direct_factor == "AR":
                                direct_direction = "concordante" if signed < 0 else "discordante_o_nulo"
                            elif category == "ESR1_gain" and direct_factor == "ESR1":
                                direct_direction = "concordante" if signed > 0 else "discordante_o_nulo"
                            individual_rows.append({
                                "region_id": region_no, "cromosoma": chrom,
                                "inicio_region_interbase_0based": min(positions), "fin_region_interbase_0based": max(positions),
                                "variant_id": variant_id, "gen_locus": event.get("gen_locus", ""),
                                "punto_union_interbase_hg38_0based": position,
                                "longitud_ancestral_bp": len(limpiar(event.get("secuencia_ancestral_inferida", ""))),
                                "secuencia_ancestral_inferida": event.get("secuencia_ancestral_inferida", ""),
                                "categoria": category, "categorias_hormonales": category,
                                "modalidad": modality, "track_id": track_id,
                                "track_name": meta.get("name", ""), "strand": meta.get("strand", ""),
                                "ontology_curie": meta.get("ontology_curie", ""),
                                "biosample_name": meta.get("biosample_name", ""), "biosample_type": meta.get("biosample_type", ""),
                                "gtex_tissue": meta.get("gtex_tissue", ""), "target": meta.get("target", ""),
                                "histone_mark": meta.get("histone_mark", ""),
                                "assay_title": meta.get("Assay title", "") or meta.get("assay_title", ""),
                                "data_source": meta.get("data_source", ""), "gene_id": meta.get("gene_id", ""),
                                "gene_symbol": meta.get("gene_symbol", ""), "factor_directo_detectado": direct_factor,
                                "radio_local_bp": int(modality_specs[modality][2]), "n_bins_locales": int(len(centers)),
                                "senal_humana_abs_integrada": float(np.sum(np.abs(base_local[:, track_idx]), dtype=np.float64)),
                                "senal_ancestral_abs_integrada": float(np.sum(np.abs(ancestral_local[:, track_idx]), dtype=np.float64)),
                                "efecto_individual_firmado_humano_menos_ancestral": signed,
                                "magnitud_abs_L1_individual": abs_l1,
                                "magnitud_absoluta_score": abs_l1,
                                "efecto_individual_abs_medio_por_bin": abs_l1 / int(len(centers)),
                                "direccion_concordante_track_factor_directo": direct_direction,
                                "metodo": "predict_sequence de REF/ALT; señal ALT alineada a coordenadas hg38 por offsets acumulados de inserciones ancestrales.",
                            })
                if event_idx % 25 == 0:
                    print(f"  [INDIVIDUALES {event_idx}/{len(ordered_events)}] {region_identity}", flush=True)
                # Keep per-variant profiles on disk; only scalar CSV rows and status remain in memory.
                delta_profiles = None
                delta_payload = None
                base_local = None
                delta_local = None
                ancestral_local = None
                alternate_sequence = None
            except Exception as exc:
                message = f"single {variant_id}: {type(exc).__name__}: {exc}"
                region_errors.append(message)
                print(f"  [FALLO INDIVIDUAL] {message}", flush=True)
                if es_error_rpc_que_debe_detener_el_run(exc):
                    stop_run = True
                    break

        # Haplotipos compuestos y no-aditividad en las mismas coordenadas/track que la suma individual.
        if not stop_run:
            for category in CATEGORIAS_HORMONALES:
                sites = category_events[category]
                if not sites:
                    continue
                print(f"  [HAPLOTIPO {category}] variantes={len(sites)}", flush=True)
                cat_data: dict[str, dict[str, Any]] = {}
                for modality, baseline in baseline_profiles.items():
                    all_centers = np.asarray(baseline["centers"], dtype=np.int64)
                    radius = int(modality_specs[modality][2])
                    category_mask = centros_union_eventos(all_centers, sites, radius)
                    centers = all_centers[category_mask]
                    base_values = np.asarray(baseline["values"][category_mask, :], dtype=np.float32)
                    if not len(centers):
                        region_errors.append(f"joint {category}/{modality}: no hay bins locales")
                        continue
                    sum_single_delta = np.zeros_like(base_values, dtype=np.float32)
                    sum_single_abs = np.zeros(base_values.shape[1], dtype=np.float64)
                    missing_events = [
                        str(site.get("variant_id", "")) for site in sites
                        if modality not in single_modalities_by_variant.get(str(site.get("variant_id", "")), set())
                    ]
                    if not missing_events:
                        for site in sites:
                            variant_id = str(site.get("variant_id", ""))
                            alt_sequence, single_events, _n_inserted = construir_secuencia_alternativa(
                                human_seq, start0, [site]
                            )
                            local_mask = np.abs(centers - int(site["punto_union_interbase_hg38_0based"])) <= radius
                            single_centers = centers[local_mask]
                            one_payloads = perfiles_locales_lote_con_checkpoint(
                                model, interval, alt_sequence, {modality: modality_specs[modality]}, ontology_terms,
                                out_dir, cache_dir, checkpoint,
                                f"v5:{region_identity}:single:{variant_id}", kind="delta",
                                target_centers_by_modality={modality: single_centers},
                                alignment_events=single_events,
                                baseline_profiles={modality: baseline_profiles[modality]},
                                failed_modalities_this_run=failed_modalities_this_run,
                            )
                            payload = one_payloads.get(modality)
                            if payload is None:
                                missing_events.append(variant_id)
                                continue
                            payload_centers = np.asarray(payload["centers"], dtype=np.int64)
                            indexes = np.searchsorted(centers, payload_centers)
                            if len(indexes) and (
                                np.any(indexes >= len(centers))
                                or not np.array_equal(centers[indexes], payload_centers)
                            ):
                                raise ValueError(f"Los bins locales de {variant_id}/{modality} no están dentro de la ventana conjunta")
                            delta_single_local = np.asarray(payload["delta"], dtype=np.float32)
                            sum_single_delta[indexes, :] += delta_single_local
                            sum_single_abs += np.sum(np.abs(delta_single_local), axis=0, dtype=np.float64)
                            # Explicitly drop the decompressed cache object before the next variant.
                            payload = None
                            one_payloads = None
                            alt_sequence = None
                            delta_single_local = None
                    if missing_events:
                        cat_data[modality] = {
                            "state": "no_calculado_faltan_individuales", "missing_events": missing_events,
                            "centers": centers, "base_values": base_values, "radius": radius,
                        }
                        continue
                    base_cat = {
                        "centers": centers, "values": base_values,
                        "metadata": baseline["metadata"], "signature": baseline["signature"],
                        "grid_start": baseline["grid_start"], "resolution": baseline["resolution"],
                        "sequence_sha256": baseline["sequence_sha256"],
                    }
                    cat_data[modality] = {
                        "state": "ready", "centers": centers, "category_mask": category_mask,
                        "base_values": base_values, "sum_single_delta": sum_single_delta,
                        "sum_single_abs": sum_single_abs, "base_cat": base_cat, "radius": radius,
                    }

                valid_modes = [m for m, item in cat_data.items() if item["state"] == "ready"]
                joint_payloads: dict[str, dict[str, Any]] = {}
                joint_error = ""
                if len(sites) == 1:
                    # Identidad exacta: para un solo indel el efecto conjunto es el individual; no se hace RPC adicional.
                    for modality in valid_modes:
                        joint_payloads[modality] = {
                            "delta": cat_data[modality]["sum_single_delta"],
                            "centers": cat_data[modality]["centers"],
                            "state": "single_variant_identity",
                        }
                elif valid_modes:
                    try:
                        joint_sequence, joint_events, _total_inserted = construir_secuencia_alternativa(
                            human_seq, start0, sites
                        )
                        joint_payloads = perfiles_locales_lote_con_checkpoint(
                            model, interval, joint_sequence,
                            {m: modality_specs[m] for m in valid_modes}, ontology_terms,
                            out_dir, cache_dir, checkpoint,
                            f"v5:{region_identity}:joint:{category}", kind="delta",
                            target_centers_by_modality={m: cat_data[m]["centers"] for m in valid_modes},
                            alignment_events=joint_events,
                            baseline_profiles={m: cat_data[m]["base_cat"] for m in valid_modes},
                            failed_modalities_this_run=failed_modalities_this_run,
                        )
                        lost = sorted(set(valid_modes) - set(joint_payloads))
                        if lost:
                            joint_error = f"AlphaGenome no devolvió perfiles para {','.join(lost)}"
                    except Exception as exc:
                        joint_error = f"{type(exc).__name__}: {exc}"
                        region_errors.append(f"joint {category}: {joint_error}")
                        print(f"  [FALLO HAPLOTIPO {category}] {joint_error}", flush=True)
                        if es_error_rpc_que_debe_detener_el_run(exc):
                            stop_run = True

                for modality, item in cat_data.items():
                    if item["state"] != "ready":
                        joint_rows.append({
                            "region_id": region_no, "cromosoma": chrom, "categoria": category,
                            "modalidad": modality, "n_variantes_conjuntas": len(sites),
                            "estado": item["state"],
                            "eventos_sin_prediccion": ";".join(item.get("missing_events", [])),
                            "nota_metodologica": "No se calcula interacción cuando faltan deltas individuales.",
                        })
                        continue
                    if modality not in joint_payloads:
                        joint_rows.append({
                            "region_id": region_no, "cromosoma": chrom, "categoria": category,
                            "modalidad": modality, "n_variantes_conjuntas": len(sites),
                            "estado": "fallo_prediccion_conjunta", "detalle_error": joint_error or "sin perfil conjunto",
                            "nota_metodologica": "No se calcula interacción al faltar el perfil ALT conjunto.",
                        })
                        continue
                    baseline = baseline_profiles[modality]
                    base_values = item["base_values"]
                    centers = item["centers"]
                    joint_delta = np.asarray(joint_payloads[modality]["delta"], dtype=np.float32)
                    sum_single_delta = item["sum_single_delta"]
                    if joint_delta.shape != sum_single_delta.shape:
                        raise ValueError(f"La delta conjunta no coincide con la suma individual en {category}/{modality}")
                    for track_idx in range(joint_delta.shape[1]):
                        meta = fila_metadata(baseline["metadata"], track_idx)
                        signature = baseline["signature"][track_idx]
                        track_id = json.dumps([modality, list(signature)], ensure_ascii=False, separators=(",", ":"))
                        stats = calcular_metricas_track(
                            modality, track_idx, track_id, meta, base_values, joint_delta,
                            sum_single_delta, float(item["sum_single_abs"][track_idx]),
                            category, len(sites), centers, int(item["radius"]),
                        )
                        direct_factor = track_factor_directo(meta)
                        direct_direction = ""
                        signed_joint = float(stats["efecto_conjunto_firmado_humano_menos_ancestral"])
                        if category == "AR_loss" and direct_factor == "AR":
                            direct_direction = "concordante" if signed_joint < 0 else "discordante_o_nulo"
                        elif category == "ESR1_gain" and direct_factor == "ESR1":
                            direct_direction = "concordante" if signed_joint > 0 else "discordante_o_nulo"
                        joint_rows.append({
                            "region_id": region_no, "cromosoma": chrom,
                            "inicio_region_interbase_0based": min(positions), "fin_region_interbase_0based": max(positions),
                            "categoria": category, "modalidad": modality, "track_id": track_id,
                            "track_name": meta.get("name", ""), "strand": meta.get("strand", ""),
                            "ontology_curie": meta.get("ontology_curie", ""),
                            "biosample_name": meta.get("biosample_name", ""), "biosample_type": meta.get("biosample_type", ""),
                            "gtex_tissue": meta.get("gtex_tissue", ""), "target": meta.get("target", ""),
                            "histone_mark": meta.get("histone_mark", ""),
                            "assay_title": meta.get("Assay title", "") or meta.get("assay_title", ""),
                            "data_source": meta.get("data_source", ""), "gene_id": meta.get("gene_id", ""),
                            "gene_symbol": meta.get("gene_symbol", ""), "factor_directo_detectado": direct_factor,
                            "n_variantes_conjuntas": len(sites), "n_bins_evaluados": len(centers),
                            "radio_ventana_bp": item["radius"], **stats,
                            "direccion_concordante_en_track_factor_directo": direct_direction,
                            "estado": "single_variant_identity" if len(sites) == 1 else "completed",
                            "mapeo_indel": "ALT bins remapeados a centros hg38 por offsets acumulados; interpolación lineal en fracciones de bin.",
                            "nota_metodologica": (
                                "La interacción es delta conjunta menos suma de deltas individuales, evaluados en los mismos bins de REF. "
                                "El remapeo multindel es propio de este script; no equivale a indel stitching oficial de un indel único."
                            ),
                        })
                if stop_run:
                    break
            if stop_run:
                pass

        elapsed_region = time.monotonic() - started_region
        state = "completado_con_errores" if region_errors else "completado_individuales_y_conjuntos"
        if stop_run:
            state = "interrumpido_rpc_reanudable"
        region_rows.append({
            "region_id": region_no, "cromosoma": chrom, "n_eventos": len(region),
            "inicio_interbase_0based": min(positions), "fin_interbase_0based": max(positions),
            "ancho_contexto_bp": len(human_seq),
            "n_AR_loss": len(category_events["AR_loss"]), "n_ESR1_gain": len(category_events["ESR1_gain"]),
            "n_ERE_gain_shared": len(category_events["ERE_gain_shared"]), "n_ESR2_change": len(category_events["ESR2_change"]),
            "estado": state, "detalle_error": " | ".join(region_errors),
            "nota_metodologica": (
                "Predicción REF humana y ALT ancestral. Deltas locales por modalidad; haplotipos de todos los indels "
                "compatibles dentro de cada categoría/región; no-aditividad = delta conjunto - suma delta individual "
                "en la unión de ventanas hg38; cada delta individual se incorpora solo dentro de su propia ventana por modalidad. "
                "La asignación ALT→REF usa offsets acumulados e interpolación en fracciones de bin; la alineación multindel es "
                "una implementación propia que requiere validación independiente."
            ),
            "duracion_segundos": round(elapsed_region, 2),
        })
        guardar_tablas_parciales()
        durations.append(time.monotonic() - started_region)
        if stop_run:
            failed_regions += 1
            break
        average = sum(durations) / len(durations)
        remaining = len(region_contexts) - region_no
        print(f"  [REGIÓN {state}] duración={durations[-1]:.1f}s; ETA regiones restantes≈{average * remaining:.1f}s", flush=True)

    comparison = comparar_efectos_individuales(individual_rows)
    escribir_csv(
        out_dir / "evaluacion_deleciones_comparacion_individual_AR_vs_ESR1.csv", comparison,
        list(dict.fromkeys(k for row in comparison for k in row.keys())) or ["modalidad", "track_id", "ganador_por_media_abs_L1"],
    )
    joint_comparison = comparar_track_rows(joint_rows)
    escribir_csv(
        out_dir / "evaluacion_deleciones_comparacion_AR_vs_ESR1.csv", joint_comparison,
        list(dict.fromkeys(k for row in joint_comparison for k in row.keys())) or ["modalidad", "track_id", "ganador_AR_vs_ESR1_por_impacto_abs"],
    )
    event_rows = []
    for candidate in candidates:
        variant_id = str(candidate.get("variant_id", ""))
        matching_tasks = [
            task for task_id, task in checkpoint.get("tasks", {}).items()
            if task_id.startswith("v5:") and f":single:{variant_id}:" in task_id
        ]
        event_rows.append({
            "variant_id": variant_id, "gen_locus": candidate.get("gen_locus", ""),
            "cromosoma": candidate["cromosoma"],
            "position_alphagenome_1based": candidate.get("position_alphagenome_1based", ""),
            "punto_union_interbase_hg38_0based": candidate.get("punto_union_interbase_hg38_0based", ""),
            "longitud_secuencia_ancestral_bp": candidate.get("longitud_secuencia_ancestral_bp", ""),
            "secuencia_ancestral_inferida": candidate.get("secuencia_ancestral_inferida", ""),
            "categorias_hormonales": ";".join(candidate.get("categorias_hormonales", [])),
            "estado": ";".join(sorted({str(task.get("status", "not_started")) for task in matching_tasks})) if matching_tasks else "not_started",
            "modalidades_solicitadas": ";".join(chosen_names),
            "metodo": ALPHAGENOME_PREDICTION_METHOD_VERSION,
        })
    escribir_csv(
        out_dir / "evaluacion_deleciones_alphagenome_eventos.csv", event_rows,
        list(dict.fromkeys(k for row in event_rows for k in row.keys())) or ["variant_id", "estado"],
    )
    guardar_estado_checkpoint(out_dir, checkpoint)
    failed_tasks = sum(
        task.get("status") in {"failed", "running", "blocked_by_failed_individuals"}
        for task_id, task in checkpoint.get("tasks", {}).items() if task_id.startswith("v5:")
    )
    elapsed = time.monotonic() - started_all
    print(f"\nAlphaGenome terminó en {elapsed/60:.1f} min.", flush=True)
    print(f"Tareas v5 fallidas/en curso: {failed_tasks}", flush=True)
    print(f"Regiones detenidas por RPC: {failed_regions}", flush=True)
    print(f"Estado reanudable: {checkpoint_path}", flush=True)
    print(f"Individuales: {out_dir / 'evaluacion_deleciones_alphagenome_tracks.csv'}", flush=True)
    print(f"Conjuntos/no-aditividad: {out_dir / 'evaluacion_deleciones_impacto_conjunto_por_track.csv'}", flush=True)
    print(f"Comparación AR vs ESR1: {out_dir / 'evaluacion_deleciones_comparacion_AR_vs_ESR1.csv'}", flush=True)
    if stop_run:
        print("El run se detuvo por un error de solicitud; al repetir el mismo comando se reutilizan perfiles cacheados.", flush=True)

def leer_csv_generico(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def cargar_filtro_previo_y_contextos(out_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    """Carga la selección PWM previa sin volver a ejecutar MAF ni el filtro C++.

    El CSV de candidatos guarda las variantes validadas; el CSV de motivos guarda
    las categorías ya calculadas. Se mantiene la agrupación regional de TODOS los
    candidatos originales y solo se vuelve a recuperar hg38 para los contextos
    AlphaGenome de las regiones que contienen eventos seleccionados.
    """
    candidate_path = out_dir / "evaluacion_deleciones_candidatos.csv"
    motif_path = out_dir / "evaluacion_deleciones_motivos_hormonales.csv"
    if not candidate_path.exists() or not motif_path.exists():
        raise FileNotFoundError(
            "Para --run-alphagenome sin repetir el filtro C++, primero deben existir "
            f"{candidate_path.name} y {motif_path.name}. Ejecuta py evaluar_deleciones.py "
            "una vez para generar los resultados PWM."
        )

    validation_rows = leer_csv_generico(candidate_path)
    motif_rows = leer_csv_generico(motif_path)
    summary_by_id = {
        limpiar(row.get("variant_id", "")): row
        for row in motif_rows if limpiar(row.get("variant_id", ""))
    }
    all_candidates: list[dict[str, Any]] = []
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in validation_rows:
        variant_id = limpiar(row.get("variant_id", ""))
        if not variant_id or variant_id in seen:
            continue
        if limpiar(row.get("estado_validacion", "")) != "lista_para_validar_referencia":
            continue
        summary = summary_by_id.get(variant_id)
        if not summary:
            continue
        categories = [
            item for item in limpiar(summary.get("categorias_hormonales", "")).split(";")
            if item in CATEGORIAS_HORMONALES
        ]
        chrom = normalizar_cromosoma(row.get("cromosoma", ""))
        point = entero(row.get("punto_union_interbase_hg38_0based"), "punto_union_interbase_hg38_0based")
        ancestral = limpiar(row.get("secuencia_ancestral_inferida", "")).upper()
        if not ancestral or set(ancestral) - BASES:
            continue
        candidate: dict[str, Any] = dict(row)
        candidate.update({
            "variant_id": variant_id,
            "gen_locus": limpiar(row.get("gen_locus", "")),
            "cromosoma": chrom,
            "punto_union_interbase_hg38_0based": point,
            "position_alphagenome_1based": entero(
                row.get("position_alphagenome_1based") or point,
                "position_alphagenome_1based",
            ),
            "longitud_secuencia_ancestral_bp": len(ancestral),
            "secuencia_ancestral_inferida": ancestral,
            "categorias_hormonales": categories,
        })
        all_candidates.append(candidate)
        if categories:
            selected.append(candidate)
        seen.add(variant_id)

    all_motif_ids = {
        limpiar(row.get("variant_id", ""))
        for row in motif_rows if limpiar(row.get("variant_id", ""))
    }
    missing = all_motif_ids - seen
    if missing:
        raise RuntimeError(
            f"Los resultados PWM contienen {len(missing)} candidatos que no se pudieron "
            "reconstruir desde evaluacion_deleciones_candidatos.csv. No se repetirá el filtro automáticamente; "
            "revisa que ambos CSV pertenezcan a la misma ejecución."
        )
    if not selected:
        return [], [], len(motif_rows)

    contexts: list[dict[str, Any]] = []
    regions = agrupar_regiones(all_candidates)
    regions_with_hits = [
        region for region in regions if any(event.get("categorias_hormonales") for event in region)
    ]
    print(
        f"Reutilizando filtro PWM previo: {len(motif_rows)} candidatos evaluados; "
        f"{len(selected)} seleccionados. No se volverán a escanear los MAF ni a ejecutar C++. ",
        flush=True,
    )
    print(
        f"Recuperando secuencias hg38 para {len(regions_with_hits)} regiones AlphaGenome…",
        flush=True,
    )
    for region_no, region in enumerate(regions_with_hits, start=1):
        chrom = region[0]["cromosoma"]
        positions = [int(item["punto_union_interbase_hg38_0based"]) for item in region]
        insertion_by_position: dict[int, int] = {}
        for event in region:
            event_pos = int(event["punto_union_interbase_hg38_0based"])
            insertion_by_position[event_pos] = max(
                insertion_by_position.get(event_pos, 0), len(event["secuencia_ancestral_inferida"])
            )
        effective_max = max(positions) + sum(
            length for pos, length in insertion_by_position.items() if pos < max(positions)
        )
        midpoint = (min(positions) + effective_max) // 2
        try:
            start0, human_seq = obtener_ventana_referencia(chrom, midpoint)
        except Exception as exc:
            raise RuntimeError(
                f"No se pudo recuperar la secuencia hg38 de la región {region_no}/{len(regions_with_hits)} "
                f"({chrom}:{min(positions)}-{max(positions)}): {exc}"
            ) from exc
        passing_region = [event for event in region if event.get("categorias_hormonales")]
        contexts.append({"events": passing_region, "start0": start0, "human_seq": human_seq})
        print(
            f"  [CONTEXTO {region_no}/{len(regions_with_hits)}] {chrom}:{min(positions)}-{max(positions)} listo",
            flush=True,
        )
    return selected, contexts, len(motif_rows)


def ejecutar_alphagenome_desde_filtro_previo(args: argparse.Namespace) -> int:
    """Modo AlphaGenome rápido: reutiliza los CSV PWM ya generados y omite el filtro C++."""
    terms = [term.strip() for term in args.ontology_terms.split(",") if term.strip()]
    if not terms:
        raise ValueError("No quedan ontology_terms válidos después de procesar --ontology-terms")
    selected, screened_contexts, n_total = cargar_filtro_previo_y_contextos(args.out_dir)
    if not selected:
        crear_salidas_alphagenome_vacias(
            args.out_dir,
            "el CSV de motivos previo no contiene candidatos seleccionados; no se llamó a AlphaGenome.",
        )
        return 0
    if args.max_events is not None:
        selected = selected[:args.max_events]
    selected_ids = {event["variant_id"] for event in selected}
    region_contexts = []
    for context in screened_contexts:
        kept = [event for event in context["events"] if event["variant_id"] in selected_ids]
        if kept:
            region_contexts.append({
                "events": kept,
                "start0": context["start0"],
                "human_seq": context["human_seq"],
            })
    print(
        f"AlphaGenome: filtro previo reutilizado ({n_total} evaluados; "
        f"{len(selected)} seleccionados para este run).",
        flush=True,
    )
    ejecutar_alphagenome(
        candidates=selected,
        region_contexts=region_contexts,
        out_dir=args.out_dir,
        key_file=args.api_key_file,
        ontology_terms=terms,
        reset_checkpoint=args.reset_checkpoint,
        modality_names=[x.strip().upper() for x in args.modalities.split(",") if x.strip()],
    )
    return 0


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Filtra motivos hormonales con C++ y evalúa efectos individuales y haplotipos de indels ancestrales con AlphaGenome, alineando perfiles ALT a coordenadas hg38."
    )
    parser.add_argument("--input", type=Path, default=here / "deleciones_todos_loci.csv", help="CSV previo de Genoma1.py")
    parser.add_argument("--out-dir", type=Path, default=here, help="Carpeta para resultados nuevos")
    parser.add_argument("--api-key-file", type=Path, default=here / "api_key.txt", help="Archivo local con la API key de AlphaGenome")
    parser.add_argument("--min-consensus-support", type=int, default=2, help="Número mínimo de primates con la misma base por columna")
    parser.add_argument("--maf-dir", type=Path, default=None, help="Carpeta con los MAF originales; por defecto, la carpeta del CSV de entrada")
    parser.add_argument("--small-indel-min-bp", type=int, default=MIN_INDEL_MAF_DEFAULT, help="Longitud mínima para recuperar huecos humanos pequeños desde MAF (por defecto 1)")
    parser.add_argument("--small-indel-max-bp", type=int, default=MAX_INDEL_MAF_DEFAULT, help="Longitud máxima para recuperar huecos humanos pequeños desde MAF (por defecto 9)")
    parser.add_argument("--no-small-indel-scan", action="store_true", help="Desactiva la búsqueda suplementaria de huecos de 1–9 pb en MAF")
    parser.add_argument("--max-events", type=int, default=None, help="Limita candidatos que ya pasaron el filtro PWM C++ y se envían a AlphaGenome")
    parser.add_argument("--run-alphagenome", action="store_true", help="Ejecuta AlphaGenome reutilizando los CSV del filtro PWM existentes; si faltan, los genera")
    parser.add_argument("--rebuild-motif-filter", action="store_true", help="Fuerza reconstrucción de candidatos y ejecución de C++ incluso con --run-alphagenome")
    parser.add_argument("--reset-checkpoint", action="store_true", help="Borra el checkpoint y la caché persistente de AlphaGenome para empezar desde cero")
    parser.add_argument("--ontology-terms", default=",".join(ONTOLOGIAS_MAMARIAS_DEFAULT), help="CURIEs separados por coma para filtrar tracks; por defecto, tejidos/células mamarias")
    parser.add_argument(
        "--modalities", default=",".join(entry[0] for entry in MODALIDADES),
        help="Modalidades AlphaGenome separadas por coma: RNA_SEQ,ATAC,DNASE,CAGE,PROCAP,CHIP_TF,CHIP_HISTONE; cada modalidad se solicita por separado y se cachea independientemente",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.min_consensus_support < 2:
        raise ValueError("--min-consensus-support debe ser al menos 2")
    if args.max_events is not None and args.max_events < 1:
        raise ValueError("--max-events debe ser positivo")
    if args.small_indel_min_bp < 1 or args.small_indel_max_bp < args.small_indel_min_bp:
        raise ValueError("El intervalo de indels pequeños debe ser válido: min >= 1 y max >= min")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Si el usuario ya ejecutó el filtro, `--run-alphagenome` reutiliza esos CSV.
    # Así no vuelve a recorrer todos los MAF ni repite el filtro C++ de varios minutos.
    prior_candidates = args.out_dir / "evaluacion_deleciones_candidatos.csv"
    prior_motifs = args.out_dir / "evaluacion_deleciones_motivos_hormonales.csv"
    if args.run_alphagenome and not args.rebuild_motif_filter and prior_candidates.exists() and prior_motifs.exists():
        return ejecutar_alphagenome_desde_filtro_previo(args)

    rows = leer_csv(args.input)
    candidates, validation = preparar_candidatos(rows, args.min_consensus_support)
    for candidate in candidates:
        if not limpiar(candidate.get("fuente_candidato", "")):
            candidate["fuente_candidato"] = "CSV_EXISTENTE_Genoma1"
    for row in validation:
        if not limpiar(row.get("fuente_candidato", "")):
            row["fuente_candidato"] = "CSV_EXISTENTE_Genoma1"

    maf_dir = (args.maf_dir or args.input.resolve().parent).resolve()
    scan_rows: list[dict[str, Any]] = []
    maf_summary: list[dict[str, Any]] = []
    supplemental_validation: list[dict[str, Any]] = []
    supplemental_added: list[dict[str, Any]] = []
    supplement_status: dict[str, str] = {}
    if args.no_small_indel_scan:
        maf_summary.append({
            "gen_locus": "", "archivo_maf": "", "ruta_maf": str(maf_dir), "estado": "BUSQUEDA_SUPLEMENTARIA_DESACTIVADA",
            "bloques": 0, "bloques_utiles": 0, "huecos_1_9": 0, "candidatos_generados": 0,
            "rechazados_por_calidad": 0, "detalle": "Se usó --no-small-indel-scan",
        })
    else:
        seed_rows, scan_rows, maf_summary = detectar_indels_pequenos_maf(
            maf_dir, args.small_indel_min_bp, args.small_indel_max_bp
        )
        if seed_rows:
            supplemental, supplemental_validation = preparar_candidatos(seed_rows, args.min_consensus_support)
            supplemental_added, supplemental_validation, supplement_status = combinar_candidatos_suplementarios(
                candidates, supplemental, supplemental_validation
            )
            for candidate in supplemental_added:
                candidate["fuente_candidato"] = "MAF_SUPLEMENTARIO_1_9BP"
            candidates.extend(supplemental_added)
            validation.extend(supplemental_validation)
        else:
            # También se registra el estado de reconstrucción aunque no haya eventos adecuados.
            _, supplemental_validation = preparar_candidatos(seed_rows, args.min_consensus_support)
            validation.extend(supplemental_validation)

    # Enlaza cada hueco MAF escaneado con el estado final de deduplicación/reconstrucción.
    for scan in scan_rows:
        scan["estado_final"] = scan.get("estado", "")
        if scan.get("estado") == "LISTO_PARA_VALIDACION":
            match = next((row for row in supplemental_validation
                          if limpiar(row.get("archivo_maf", "")) == limpiar(scan.get("archivo_maf", ""))
                          and limpiar(row.get("bloque_maf", "")) == limpiar(scan.get("bloque_maf", ""))
                          and limpiar(row.get("columna_inicio_bloque_1based", "")) == limpiar(scan.get("columna_inicio_bloque_1based", ""))), None)
            if match is not None:
                scan["variant_id"] = match.get("variant_id", "")
                scan["estado_final"] = match.get("estado_validacion", scan["estado"])
                scan["detalle_final"] = match.get("motivo_validacion", "")
                if scan.get("variant_id") in supplement_status:
                    scan["estado_final"] = supplement_status[scan["variant_id"]]
                    scan["detalle_final"] = supplement_status[scan["variant_id"]]

    validation_path = args.out_dir / "evaluacion_deleciones_candidatos.csv"
    validation_fields = list(dict.fromkeys(CAMPOS_CANDIDATO + ["fila_entrada"]))
    escribir_csv(validation_path, validation, validation_fields)
    scan_fields = list(dict.fromkeys(k for row in scan_rows for k in row.keys())) or [
        "gen_locus", "archivo_maf", "bloque_maf", "longitud_indel_bp", "estado", "detalle"
    ]
    escribir_csv(args.out_dir / "evaluacion_indels_pequenos_maf.csv", scan_rows, scan_fields)
    maf_summary_fields = list(dict.fromkeys(k for row in maf_summary for k in row.keys())) or [
        "gen_locus", "archivo_maf", "estado", "detalle"
    ]
    escribir_csv(args.out_dir / "evaluacion_indels_pequenos_maf_resumen.csv", maf_summary, maf_summary_fields)

    n_valid = sum(row.get("estado_validacion") == "lista_para_validar_referencia" for row in validation)
    n_reject = sum(str(row.get("estado_validacion", "")).startswith("rechazada") for row in validation)
    n_small_found = sum(row.get("estado") == "LISTO_PARA_VALIDACION" for row in scan_rows)
    n_small_accepted = len(supplemental_added)
    n_small_failed = sum(row.get("estado_final", "").startswith(("rechazada", "conflicto", "duplicado")) for row in scan_rows)
    print(f"Filas de entrada: {len(rows)} | candidatos reconstruidos del CSV: {len([x for x in validation if x.get('fuente_candidato') == 'CSV_EXISTENTE_Genoma1' and x.get('estado_validacion') == 'lista_para_validar_referencia'])} | rechazados: {n_reject}")
    if not args.no_small_indel_scan:
        print(f"Búsqueda MAF suplementaria ({args.small_indel_min_bp}-{args.small_indel_max_bp} pb): {n_small_found} huecos candidatos de calidad suficiente | {n_small_accepted} añadidos | {n_small_failed} duplicados/conflictos/rechazados")
        print(f"Informe de indels pequeños: {args.out_dir / 'evaluacion_indels_pequenos_maf.csv'}")
        print(f"Resumen de archivos MAF: {args.out_dir / 'evaluacion_indels_pequenos_maf_resumen.csv'}")
    print(f"Total de candidatos listos para el filtro PWM: {len(candidates)}")
    print(f"Validación guardada en: {validation_path}")

    selected, screened_contexts = prefiltrar_por_motivos_cpp(candidates, args.out_dir)

    if not args.run_alphagenome:
        print("Filtro PWM C++ completado; no se llamó a AlphaGenome. Para ejecutar predicciones añade --run-alphagenome.")
        return 0

    terms = [term.strip() for term in args.ontology_terms.split(",") if term.strip()]
    if not terms:
        raise ValueError("No quedan ontology_terms válidos después de procesar --ontology-terms")

    if not selected:
        crear_salidas_alphagenome_vacias(
            args.out_dir,
            "ningún candidato presentó un cambio de motivo hormonal por encima de los umbrales; "
            "no se inicializó ni se llamó a AlphaGenome.",
        )
        return 0

    if args.max_events is not None:
        selected = selected[:args.max_events]
    selected_ids = {event["variant_id"] for event in selected}
    region_contexts = []
    for context in screened_contexts:
        kept = [event for event in context["events"] if event["variant_id"] in selected_ids]
        if kept:
            region_contexts.append({
                "events": kept,
                "start0": context["start0"],
                "human_seq": context["human_seq"],
            })

    excluded_count = len(candidates) - len({event["variant_id"] for context in screened_contexts for event in context["events"]})
    print(
        f"\nFiltro previo terminado. Candidatos originales: {len(candidates)}; "
        f"seleccionados por motivos: {len({event['variant_id'] for context in screened_contexts for event in context['events']})}; "
        f"excluidos antes de AlphaGenome: {excluded_count}; "
        f"seleccionados para este run: {len(selected)}.",
        flush=True,
    )
    if len(selected) < 1:
        crear_salidas_alphagenome_vacias(args.out_dir, "el límite de eventos no dejó candidatos para analizar.")
        return 0

    ejecutar_alphagenome(
        candidates=selected,
        region_contexts=region_contexts,
        out_dir=args.out_dir,
        key_file=args.api_key_file,
        ontology_terms=terms,
        reset_checkpoint=args.reset_checkpoint,
        modality_names=[x.strip().upper() for x in args.modalities.split(",") if x.strip()],
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrumpido por el usuario. Las predicciones completadas permanecen en la caché persistente; las tareas incompletas se reintentarán al ejecutar de nuevo.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as error:
        print(f"ERROR: {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1)
