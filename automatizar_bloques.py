import csv
import hashlib
import json
import math
import os
import re
import ssl
import sys
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

import certifi
import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client

CARPETA = Path(__file__).resolve().parent
ARCHIVO_ENTRADA = CARPETA / "variantes_todos_loci.csv"
ARCHIVO_TRACKS = CARPETA / "resultados_impacto_conjunto_por_track.csv"
ARCHIVO_INDIVIDUAL = CARPETA / "resultados_impacto_individual_por_track.csv"
ARCHIVO_COMPARACION = CARPETA / "resultados_impacto_conjunto_comparacion.csv"
ARCHIVO_MODALIDADES = CARPETA / "resultados_impacto_conjunto_resumen_modalidad.csv"
ARCHIVO_REGIONES = CARPETA / "resultados_impacto_conjunto_regiones.csv"

UMBRAL_MOTIVO = 0.70
DELTA_MOTIVO = 0.05
DNA = {"A", "C", "G", "T"}
CATEGORIAS = ("AR_loss", "ESR1_gain")
LONGITUD_SECUENCIA = 2**20
MAX_SPAN_HAPLOTIPO = 850_000
RADIO_POR_MODALIDAD = {
    "RNA_SEQ": 100_000,
    "CAGE": 50_000,
    "PROCAP": 50_000,
    "ATAC": 5_000,
    "DNASE": 5_000,
    "CHIP_TF": 5_000,
    "CHIP_HISTONE": 5_000,
}
MODALIDADES = (
    ("RNA_SEQ", dna_client.OutputType.RNA_SEQ, "rna_seq"),
    ("ATAC", dna_client.OutputType.ATAC, "atac"),
    ("DNASE", dna_client.OutputType.DNASE, "dnase"),
    ("CAGE", dna_client.OutputType.CAGE, "cage"),
    ("PROCAP", dna_client.OutputType.PROCAP, "procap"),
    ("CHIP_TF", dna_client.OutputType.CHIP_TF, "chip_tf"),
    ("CHIP_HISTONE", dna_client.OutputType.CHIP_HISTONE, "chip_histone"),
)
SALIDA_SOLICITADA = [entry[1] for entry in MODALIDADES]
CAMPOS_TRACK = (
    "name", "strand", "ontology_curie", "biosample_name", "biosample_type",
    "gtex_tissue", "target", "histone_mark", "Assay title", "assay_title",
    "data_source", "gene_id", "gene_symbol"
)
PALABRAS_TEJIDOS = re.compile(
    r"breast|mammary|adipose|ovary|ovarian|uterus|uterine|testis|testicular|"
    r"pituitary|prostate|placenta|mcf[- ]?7|t47d|mda[- ]?mb",
    re.IGNORECASE,
)
SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
HTTP_HEADERS = {"Accept": "application/json", "User-Agent": "AlphaGenomeJointHaplotypeAnalysis/2.0"}


def clean(value):
    if value is None:
        return ""
    try:
        if isinstance(value, float) and math.isnan(value):
            return ""
    except TypeError:
        pass
    result = str(value).strip()
    return "" if result.lower() in {"nan", "none", "<na>"} else result


def get_api_key():
    path = CARPETA / "api_key.txt"
    if path.exists():
        key = path.read_text(encoding="utf-8").strip()
        if key:
            return key
    return os.environ.get("ALPHAGENOME_API_KEY", "").strip()


def request_json(url, timeout=90):
    request = Request(url, headers=HTTP_HEADERS)
    with urlopen(request, timeout=timeout, context=SSL_CONTEXT) as response:
        if response.status != 200:
            raise RuntimeError(f"HTTP {response.status} al consultar el servicio de secuencia.")
        return json.loads(response.read().decode("utf-8"))


def normalizar_cromosoma(value):
    chrom = clean(value)
    if not chrom:
        raise ValueError("Cromosoma vacío.")
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.upper()
    if chrom in {"M", "MT"}:
        return "chrM"
    if chrom in {"X", "Y"}:
        return "chr" + chrom
    if chrom.isdigit() and 1 <= int(chrom) <= 22:
        return "chr" + str(int(chrom))
    raise ValueError(f"Cromosoma no válido: {value!r}")


def dna_base(value, field):
    base = clean(value).upper()
    if base not in DNA:
        raise ValueError(f"{field} no es una base A/C/G/T: {value!r}")
    return base


def number(row, field):
    try:
        value = float(clean(row.get(field, "")))
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def create_candidates():
    required = {
        "coordenada_hg38_1based", "base_ancestral_inferida", "base_humana_maf", "cromosoma",
        "AR_PWM_ancestral", "AR_PWM_delta_ancestral_menos_humano",
        "ESR1_PWM_humano", "ESR1_PWM_delta_ancestral_menos_humano", "ESR1_motif_id"
    }
    candidates = []
    omitted = 0
    seen = set()
    with open(ARCHIVO_ENTRADA, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Faltan columnas obligatorias: " + ", ".join(sorted(missing)))
        for row_number, row in enumerate(reader, start=2):
            ar_anc = number(row, "AR_PWM_ancestral")
            ar_delta = number(row, "AR_PWM_delta_ancestral_menos_humano")
            er_human = number(row, "ESR1_PWM_humano")
            er_delta = number(row, "ESR1_PWM_delta_ancestral_menos_humano")
            er_motif = clean(row.get("ESR1_motif_id"))
            categories = []
            if ar_anc is not None and ar_delta is not None and ar_anc >= UMBRAL_MOTIVO and ar_delta >= DELTA_MOTIVO:
                categories.append("AR_loss")
            if er_motif and er_human is not None and er_delta is not None and er_human >= UMBRAL_MOTIVO and er_delta <= -DELTA_MOTIVO:
                categories.append("ESR1_gain")
            if not categories:
                continue
            try:
                pos = int(clean(row.get("coordenada_hg38_1based")))
                if pos < 1:
                    raise ValueError("Coordenada no positiva.")
                chrom = normalizar_cromosoma(row.get("cromosoma"))
                human = dna_base(row.get("base_humana_maf"), "base_humana_maf")
                ancestral = dna_base(row.get("base_ancestral_inferida"), "base_ancestral_inferida")
                if human == ancestral:
                    raise ValueError("Alelo humano igual al ancestral.")
            except (ValueError, TypeError):
                omitted += 1
                continue
            for category in categories:
                key = (category, chrom, pos, human, ancestral)
                if key in seen:
                    continue
                seen.add(key)
                candidates.append({
                    "fila_csv": row_number,
                    "cromosoma": chrom,
                    "posicion": pos,
                    "gen_locus": clean(row.get("gen_locus")),
                    "categoria": category,
                    "humano": human,
                    "ancestral": ancestral,
                })

    by_position = defaultdict(set)
    for item in candidates:
        by_position[(item["cromosoma"], item["posicion"])].add((item["humano"], item["ancestral"]))
    conflict_positions = {key for key, alleles in by_position.items() if len(alleles) > 1}
    if conflict_positions:
        before = len(candidates)
        candidates = [x for x in candidates if (x["cromosoma"], x["posicion"]) not in conflict_positions]
        print(f"[AVISO] Se excluyeron {before - len(candidates)} entradas en {len(conflict_positions)} coordenadas con alelos conflictivos.")
    print(f"Candidatos para predicción conjunta: {len(candidates)}")
    if omitted:
        print(f"Entradas omitidas por alelos/coordenadas inválidos: {omitted}")
    for category in CATEGORIAS:
        print(f"  {category}: {sum(x['categoria'] == category for x in candidates)} variantes candidatas")
    return candidates


def fetch_sequence(chrom, start0, end0):
    url = (
        "https://api.genome.ucsc.edu/getData/sequence?"
        f"genome=hg38;chrom={chrom};start={start0};end={end0}"
    )
    data = request_json(url)
    sequence = clean(data.get("dna")).upper()
    expected = end0 - start0
    if len(sequence) != expected:
        raise RuntimeError(
            f"UCSC devolvió {len(sequence):,} bases para {chrom}:{start0}-{end0}; se esperaban {expected:,}."
        )
    bad = set(sequence) - DNA - {"N"}
    if bad:
        raise RuntimeError(f"La secuencia contiene caracteres inesperados: {sorted(bad)}")
    return sequence


def partition_regions(candidates):
    by_chrom = defaultdict(list)
    for item in candidates:
        by_chrom[item["cromosoma"]].append(item)
    regions = []
    for chrom in sorted(by_chrom):
        items = sorted(by_chrom[chrom], key=lambda x: x["posicion"])
        current = []
        first_pos = None
        for item in items:
            pos = item["posicion"]
            if current and pos - first_pos > MAX_SPAN_HAPLOTIPO:
                regions.append(current)
                current = []
                first_pos = None
            if not current:
                first_pos = pos
            current.append(item)
        if current:
            regions.append(current)
    return regions


def output_ontology_terms(model):
    fallback = [
        "UBERON:0000310", "UBERON:0001347", "UBERON:0000992",
        "UBERON:0000995", "UBERON:0000473", "UBERON:0000007"
    ]
    terms = {"UBERON:0000310"}
    try:
        metadata = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
        for modality, _, attr in MODALIDADES:
            frame = getattr(metadata, attr, None)
            if frame is None or not hasattr(frame, "iterrows"):
                continue
            for _, row in frame.iterrows():
                values = [clean(v) for v in row.tolist()]
                haystack = " ".join(values)
                if PALABRAS_TEJIDOS.search(haystack):
                    term = clean(row.get("ontology_curie", ""))
                    if term and ":" in term:
                        terms.add(term)
        if len(terms) == 1:
            terms.update(fallback)
    except Exception as error:
        print(f"[AVISO] No se pudo inspeccionar la metadata de tejidos: {type(error).__name__}: {error}")
        terms.update(fallback)
    result = sorted(terms)
    print(f"Ontologías de predicción seleccionadas: {len(result)}")
    return result


def output_blocks(output):
    blocks = {}
    for modality, _, attr in MODALIDADES:
        td = getattr(output, attr, None)
        if td is None or getattr(td, "values", None) is None:
            continue
        array = np.asarray(td.values)
        metadata = td.metadata.copy()
        if array.ndim == 1:
            if len(metadata) == 1 and len(array) > 1:
                array = array.reshape(-1, 1)
            elif len(array) == len(metadata):
                array = array.reshape(1, -1)
            else:
                raise ValueError(f"Forma inesperada en {modality}: {array.shape}; tracks={len(metadata)}")
        if array.ndim != 2:
            raise ValueError(f"Se esperaba una matriz de dos dimensiones para {modality}, recibida {array.shape}.")
        if array.shape[1] != len(metadata):
            raise ValueError(f"No coinciden las columnas de {modality} con su metadata.")
        if array.size == 0 or len(metadata) == 0:
            continue
        interval = getattr(td, "interval", None)
        if interval is None:
            raise ValueError(f"AlphaGenome no incluyó el intervalo de salida en {modality}.")
        resolution = int(getattr(td, "resolution", 1))
        blocks[modality] = {
            "values": np.nan_to_num(array, nan=0.0, posinf=0.0, neginf=0.0),
            "metadata": metadata,
            "signature": track_signature(metadata),
            "start": int(interval.start),
            "resolution": resolution,
        }
    return blocks


def track_signature(metadata):
    result = []
    for _, row in metadata.iterrows():
        item = {field: clean(row.get(field, "")) for field in CAMPOS_TRACK if field in metadata.columns}
        if not item.get("name"):
            item["name"] = clean(row.get("name", ""))
        if not item.get("strand"):
            item["strand"] = clean(row.get("strand", "."))
        result.append(json.dumps(item, sort_keys=True, ensure_ascii=False))
    return tuple(result)


def site_mask(n_bins, start0, resolution, positions, radius):
    centers = start0 + np.arange(n_bins, dtype=np.int64) * resolution + (resolution // 2)
    mask = np.zeros(n_bins, dtype=bool)
    for pos in positions:
        center0 = int(pos) - 1
        mask |= (centers >= center0 - radius) & (centers < center0 + radius + 1)
    if not mask.any() and n_bins:
        for pos in positions:
            idx = int(round(((int(pos) - 1) - start0) / max(1, resolution)))
            if 0 <= idx < n_bins:
                mask[idx] = True
    return mask


def track_fields(metadata_row):
    out = {field: clean(metadata_row.get(field, "")) for field in CAMPOS_TRACK}
    searchable = " ".join(out.values()).lower()
    factor = ""
    if re.search(r"androgen receptor|(?<![a-z0-9])ar(?![a-z0-9])", searchable):
        factor = "AR"
    elif re.search(r"\besr1\b|estrogen receptor alpha|oestrogen receptor alpha|er[- ]alpha|erα", searchable):
        factor = "ESR1"
    out["factor_directo_detectado"] = factor
    return out


def add_prediction_batch(model, sequence_jobs, interval, ontology_terms, callback, batch_size=3, workers=2):
    for begin in range(0, len(sequence_jobs), batch_size):
        batch = sequence_jobs[begin:begin + batch_size]
        seqs = [job["sequence"] for job in batch]
        intervals = [interval] * len(batch)
        outputs = model.predict_sequences(
            sequences=seqs,
            intervals=intervals,
            requested_outputs=SALIDA_SOLICITADA,
            ontology_terms=ontology_terms,
            progress_bar=False,
            max_workers=min(workers, len(batch)),
        )
        if len(outputs) != len(batch):
            raise RuntimeError(f"AlphaGenome devolvió {len(outputs)} resultados para {len(batch)} secuencias.")
        for job, output in zip(batch, outputs):
            callback(job, output)
        completed = min(begin + len(batch), len(sequence_jobs))
        print(f"    Predicciones: {completed}/{len(sequence_jobs)}", end="\r", flush=True)
    if sequence_jobs:
        print()


def write_csv(path, rows):
    if not rows:
        with open(path, "w", newline="", encoding="utf-8-sig") as handle:
            handle.write("sin_resultados\n")
        print(f"[AVISO] {path.name}: no hubo filas; se dejó constancia en el archivo.")
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize_tracks(track_rows):
    grouped = {}
    for row in track_rows:
        key = (row["categoria"], row["modalidad"], row["track_id"])
        if key not in grouped:
            grouped[key] = {
                "categoria": row["categoria"], "modalidad": row["modalidad"],
                "track_id": row["track_id"], "metadata": row,
                "n_regiones": 0, "n_variantes": 0,
                "sum_abs_joint": 0.0, "sum_signed_joint": 0.0,
                "sum_abs_individual": 0.0, "sum_signed_individual": 0.0,
                "sum_abs_interaction": 0.0, "sum_signed_interaction": 0.0,
                "sum_abs_human_signal": 0.0, "sum_abs_ancestral_signal": 0.0,
            }
        g = grouped[key]
        g["n_regiones"] += 1
        g["n_variantes"] += int(row["n_variantes_conjuntas"])
        for out_key, in_key in (
            ("sum_abs_joint", "efecto_conjunto_abs_L1"),
            ("sum_signed_joint", "efecto_conjunto_firmado_humano_menos_ancestral"),
            ("sum_abs_individual", "suma_efectos_individuales_abs_L1"),
            ("sum_signed_individual", "suma_efectos_individuales_firmados"),
            ("sum_abs_interaction", "no_aditividad_interaccion_abs_L1"),
            ("sum_signed_interaction", "no_aditividad_interaccion_firmada"),
            ("sum_abs_human_signal", "senal_humana_abs_integrada"),
            ("sum_abs_ancestral_signal", "senal_ancestral_abs_integrada"),
        ):
            try:
                g[out_key] += float(row[in_key])
            except (ValueError, TypeError, KeyError):
                pass
    rows = []
    for g in grouped.values():
        metadata = g["metadata"]
        n = g["n_variantes"]
        rows.append({
            "categoria": g["categoria"], "modalidad": g["modalidad"], "track_id": g["track_id"],
            "factor_directo_detectado": metadata.get("factor_directo_detectado", ""),
            "track_name": metadata.get("name", ""), "strand": metadata.get("strand", ""),
            "ontology_curie": metadata.get("ontology_curie", ""),
            "biosample_name": metadata.get("biosample_name", ""),
            "biosample_type": metadata.get("biosample_type", ""),
            "gtex_tissue": metadata.get("gtex_tissue", ""), "target": metadata.get("target", ""),
            "n_regiones": g["n_regiones"], "n_variantes_en_regiones": n,
            "suma_abs_efecto_conjunto": g["sum_abs_joint"],
            "suma_firmada_efecto_conjunto_humano_menos_ancestral": g["sum_signed_joint"],
            "suma_abs_efectos_individuales": g["sum_abs_individual"],
            "suma_firmada_efectos_individuales": g["sum_signed_individual"],
            "suma_abs_no_aditividad_interaccion": g["sum_abs_interaction"],
            "suma_firmada_no_aditividad_interaccion": g["sum_signed_interaction"],
            "efecto_conjunto_abs_por_variante": g["sum_abs_joint"] / n if n else "",
            "efectos_individuales_abs_por_variante": g["sum_abs_individual"] / n if n else "",
            "interaccion_abs_por_variante": g["sum_abs_interaction"] / n if n else "",
            "interaccion_sobre_suma_individual_abs": (
                g["sum_abs_interaction"] / g["sum_abs_individual"] if g["sum_abs_individual"] else ""
            ),
        })
    return rows


def make_comparison(summary_rows):
    indexed = {(row["categoria"], row["modalidad"], row["track_id"]): row for row in summary_rows}
    keys = sorted({(row["modalidad"], row["track_id"]) for row in summary_rows})
    out = []
    for modality, track_id in keys:
        ar = indexed.get(("AR_loss", modality, track_id))
        es = indexed.get(("ESR1_gain", modality, track_id))
        if ar is None or es is None:
            continue
        ar_total = float(ar["suma_abs_efecto_conjunto"])
        es_total = float(es["suma_abs_efecto_conjunto"])
        ar_per_variant = ar.get("efecto_conjunto_abs_por_variante", "")
        es_per_variant = es.get("efecto_conjunto_abs_por_variante", "")
        out.append({
            "modalidad": modality, "track_id": track_id,
            "track_name": ar.get("track_name", ""), "strand": ar.get("strand", ""),
            "ontology_curie": ar.get("ontology_curie", ""),
            "biosample_name": ar.get("biosample_name", ""),
            "biosample_type": ar.get("biosample_type", ""),
            "gtex_tissue": ar.get("gtex_tissue", ""), "target": ar.get("target", ""),
            "n_regiones_AR_loss": ar.get("n_regiones", 0), "n_regiones_ESR1_gain": es.get("n_regiones", 0),
            "n_variantes_AR_loss": ar.get("n_variantes_en_regiones", 0),
            "n_variantes_ESR1_gain": es.get("n_variantes_en_regiones", 0),
            "impacto_conjunto_abs_AR_loss": ar_total,
            "impacto_conjunto_abs_ESR1_gain": es_total,
            "diferencia_abs_AR_menos_ESR1": ar_total - es_total,
            "ratio_abs_AR_sobre_ESR1": ar_total / es_total if es_total else "",
            "impacto_abs_AR_por_variante": ar_per_variant,
            "impacto_abs_ESR1_por_variante": es_per_variant,
            "ganador_por_impacto_conjunto": "AR_loss" if ar_total > es_total else ("ESR1_gain" if es_total > ar_total else "EMPATE"),
            "interaccion_abs_AR_loss": ar.get("suma_abs_no_aditividad_interaccion", ""),
            "interaccion_abs_ESR1_gain": es.get("suma_abs_no_aditividad_interaccion", ""),
        })
    return out


def summarize_modalities(summary_rows):
    output = []
    for category in CATEGORIAS:
        for modality in sorted({r["modalidad"] for r in summary_rows}):
            values = [r for r in summary_rows if r["categoria"] == category and r["modalidad"] == modality]
            if not values:
                continue
            cumulative = [float(r["suma_abs_efecto_conjunto"]) for r in values]
            per_variant = [float(r["efecto_conjunto_abs_por_variante"]) for r in values if r["efecto_conjunto_abs_por_variante"] != ""]
            interaction = [float(r["suma_abs_no_aditividad_interaccion"]) for r in values]
            direct_count = sum(r.get("factor_directo_detectado") == ("AR" if category == "AR_loss" else "ESR1") for r in values)
            output.append({
                "categoria": category, "modalidad": modality, "n_tracks_con_resultados": len(values),
                "n_tracks_factor_directo_esperado": direct_count,
                "mediana_impacto_conjunto_abs_por_track": float(np.median(cumulative)),
                "media_impacto_conjunto_abs_por_track": float(np.mean(cumulative)),
                "mediana_impacto_abs_por_variante_y_track": float(np.median(per_variant)) if per_variant else "",
                "mediana_no_aditividad_abs_por_track": float(np.median(interaction)),
            })
    return output


def process_region(model, region_items, region_number, ontology_terms, track_rows, individual_rows, region_rows):
    chrom = region_items[0]["cromosoma"]
    sites = {}
    for item in region_items:
        key = (item["posicion"], item["humano"], item["ancestral"])
        sites[key] = {"posicion": item["posicion"], "humano": item["humano"], "ancestral": item["ancestral"]}
    positions = sorted({key[0] for key in sites})
    if not positions:
        return
    by_position = defaultdict(set)
    for pos, human, ancestral in sites:
        by_position[pos].add((human, ancestral))
    conflicting = {pos for pos, alleles in by_position.items() if len(alleles) > 1}
    if conflicting:
        region_items = [x for x in region_items if x["posicion"] not in conflicting]
        sites = {(x["posicion"], x["humano"], x["ancestral"]): x for x in region_items}
        positions = sorted({x["posicion"] for x in region_items})
        if not positions:
            print(f"[AVISO] Región {region_number}: todas las posiciones tenían alelos incompatibles.")
            return

    center = (positions[0] + positions[-1]) // 2
    start0 = max(0, center - LONGITUD_SECUENCIA // 2)
    end0 = start0 + LONGITUD_SECUENCIA
    raw_sequence = fetch_sequence(chrom, start0, end0)
    human_chars = list(raw_sequence)
    mismatch_count = 0
    site_offset = {}
    for pos, human, ancestral in sorted(sites):
        offset = pos - 1 - start0
        if offset < 0 or offset >= len(human_chars):
            raise ValueError(f"La variante {chrom}:{pos} quedó fuera del intervalo de contexto.")
        if raw_sequence[offset] != human:
            mismatch_count += 1
        human_chars[offset] = human
        site_offset[(pos, human, ancestral)] = offset
    human_sequence = "".join(human_chars)
    interval = genome.Interval(chromosome=chrom, start=start0, end=end0)

    by_category = {cat: [] for cat in CATEGORIAS}
    cats_by_variant = defaultdict(set)
    for item in region_items:
        variant_key = (item["posicion"], item["humano"], item["ancestral"])
        by_category[item["categoria"]].append(variant_key)
        cats_by_variant[variant_key].add(item["categoria"])
    for cat in CATEGORIAS:
        by_category[cat] = sorted(set(by_category[cat]))

    print(
        f"\n[REGIÓN {region_number}] {chrom}:{positions[0]}-{positions[-1]} | "
        f"sitios={len(sites)} | AR_loss={len(by_category['AR_loss'])} | ESR1_gain={len(by_category['ESR1_gain'])}"
    )
    if mismatch_count:
        print(f"  [INFO] {mismatch_count} sitios: la base de GRCh38 no era la base humana del CSV; se simuló explícitamente el alelo humano del CSV.")

    baseline_output = model.predict_sequence(
        sequence=human_sequence,
        interval=interval,
        requested_outputs=SALIDA_SOLICITADA,
        ontology_terms=ontology_terms,
    )
    baseline_blocks = output_blocks(baseline_output)
    state = {cat: {} for cat in CATEGORIAS}
    category_positions = {cat: [key[0] for key in by_category[cat]] for cat in CATEGORIAS}
    for cat in CATEGORIAS:
        if not by_category[cat]:
            continue
        for modality, block in baseline_blocks.items():
            radius = RADIO_POR_MODALIDAD[modality]
            mask = site_mask(block["values"].shape[0], block["start"], block["resolution"], category_positions[cat], radius)
            if not mask.any():
                continue
            ref_masked = block["values"][mask, :].astype(np.float32, copy=True)
            state[cat][modality] = {
                "mask": mask,
                "ref": ref_masked,
                "sum_single_delta": np.zeros_like(ref_masked, dtype=np.float32),
                "sum_single_abs": np.zeros_like(ref_masked, dtype=np.float32),
                "metadata": block["metadata"].copy(),
                "signature": block["signature"],
                "n_bins": int(mask.sum()),
                "radius": radius,
            }
    del baseline_output, baseline_blocks

    # Prediction of each individual variant on the same human genetic background.
    single_jobs = []
    for variant_key in sorted(cats_by_variant):
        seq_chars = list(human_sequence)
        seq_chars[site_offset[variant_key]] = variant_key[2]
        seq = "".join(seq_chars)
        single_jobs.append({
            "sequence": seq,
            "variant_key": variant_key,
            "categories": sorted(cats_by_variant[variant_key]),
        })

    def single_callback(job, output):
        blocks = output_blocks(output)
        pos, human, ancestral = job["variant_key"]
        variant_id = f"{chrom}:{pos}:{human}>{ancestral}"
        for cat in job["categories"]:
            for modality, group_state in state[cat].items():
                block = blocks.get(modality)
                if block is None or block["signature"] != group_state["signature"]:
                    continue
                alt_masked = block["values"][group_state["mask"], :].astype(np.float32, copy=False)
                delta = group_state["ref"] - alt_masked
                group_state["sum_single_delta"] += delta
                group_state["sum_single_abs"] += np.abs(delta)
                meta = group_state["metadata"]
                for idx in range(delta.shape[1]):
                    fields = track_fields(meta.iloc[idx].to_dict())
                    d = delta[:, idx]
                    individual_rows.append({
                        "region_id": region_number, "cromosoma": chrom,
                        "inicio_region_hg38_1based": positions[0], "fin_region_hg38_1based": positions[-1],
                        "categoria": cat, "variante": variant_id, "modalidad": modality,
                        "track_id": group_state["signature"][idx], **fields,
                        "n_variantes_del_grupo_en_region": len(by_category[cat]),
                        "n_bins_evaluados": group_state["n_bins"],
                        "radio_ventana_bp": group_state["radius"],
                        "efecto_individual_firmado_humano_menos_ancestral": float(np.sum(d)),
                        "efecto_individual_abs_L1": float(np.sum(np.abs(d))),
                        "efecto_individual_abs_medio_por_bin": float(np.mean(np.abs(d))) if len(d) else 0.0,
                    })
        del blocks

    print(f"  Efectos individuales condicionados al fondo humano: {len(single_jobs)} variantes únicas")
    add_prediction_batch(model, single_jobs, interval, ontology_terms, single_callback, batch_size=3, workers=2)

    # Prediction of the two joint haplotypes: only the evaluated group reverts to ancestral.
    combined_map = {}
    for cat in CATEGORIAS:
        if not by_category[cat]:
            continue
        chars = list(human_sequence)
        for variant_key in by_category[cat]:
            chars[site_offset[variant_key]] = variant_key[2]
        seq = "".join(chars)
        digest = hashlib.sha1(seq.encode("ascii", errors="ignore")).hexdigest()
        if digest not in combined_map:
            combined_map[digest] = {"sequence": seq, "categories": []}
        combined_map[digest]["categories"].append(cat)
    combined_jobs = list(combined_map.values())

    def combined_callback(job, output):
        blocks = output_blocks(output)
        for cat in job["categories"]:
            for modality, group_state in state[cat].items():
                block = blocks.get(modality)
                if block is None or block["signature"] != group_state["signature"]:
                    continue
                alt_masked = block["values"][group_state["mask"], :].astype(np.float32, copy=False)
                joint_delta = group_state["ref"] - alt_masked
                additive_delta = group_state["sum_single_delta"]
                residual = joint_delta - additive_delta
                meta = group_state["metadata"]
                for idx in range(joint_delta.shape[1]):
                    fields = track_fields(meta.iloc[idx].to_dict())
                    jd = joint_delta[:, idx]
                    sd = additive_delta[:, idx]
                    sa = group_state["sum_single_abs"][:, idx]
                    inter = residual[:, idx]
                    ref_signal = group_state["ref"][:, idx]
                    alt_signal = alt_masked[:, idx]
                    joint_signed = float(np.sum(jd))
                    joint_abs = float(np.sum(np.abs(jd)))
                    single_signed = float(np.sum(sd))
                    single_abs = float(np.sum(sa))
                    interaction_abs = float(np.sum(np.abs(inter)))
                    direct_expected = ""
                    if cat == "AR_loss" and fields["factor_directo_detectado"] == "AR":
                        direct_expected = joint_signed < 0
                    elif cat == "ESR1_gain" and fields["factor_directo_detectado"] == "ESR1":
                        direct_expected = joint_signed > 0
                    track_rows.append({
                        "region_id": region_number, "cromosoma": chrom,
                        "inicio_region_hg38_1based": positions[0], "fin_region_hg38_1based": positions[-1],
                        "categoria": cat, "modalidad": modality,
                        "track_id": group_state["signature"][idx], **fields,
                        "n_variantes_conjuntas": len(by_category[cat]),
                        "n_bins_evaluados": group_state["n_bins"],
                        "radio_ventana_bp": group_state["radius"],
                        "senal_humana_abs_integrada": float(np.sum(np.abs(ref_signal))),
                        "senal_ancestral_abs_integrada": float(np.sum(np.abs(alt_signal))),
                        "efecto_conjunto_firmado_humano_menos_ancestral": joint_signed,
                        "efecto_conjunto_abs_L1": joint_abs,
                        "efecto_conjunto_abs_medio_por_bin": float(np.mean(np.abs(jd))) if len(jd) else 0.0,
                        "efecto_conjunto_relativo_pct_vs_senal_humana": 100.0 * joint_signed / (float(np.sum(np.abs(ref_signal))) + 1e-8),
                        "suma_efectos_individuales_firmados": single_signed,
                        "suma_efectos_individuales_abs_L1": single_abs,
                        "no_aditividad_interaccion_firmada": float(np.sum(inter)),
                        "no_aditividad_interaccion_abs_L1": interaction_abs,
                        "interaccion_abs_sobre_suma_individual_abs": interaction_abs / single_abs if single_abs else "",
                        "efecto_conjunto_abs_por_variante": joint_abs / len(by_category[cat]) if by_category[cat] else "",
                        "suma_individual_abs_por_variante": single_abs / len(by_category[cat]) if by_category[cat] else "",
                        "interaccion_abs_por_variante": interaction_abs / len(by_category[cat]) if by_category[cat] else "",
                        "direccion_concordante_en_track_factor_directo": direct_expected,
                    })
        del blocks

    print(f"  Haplotípos conjuntos a predecir: {len(combined_jobs)}")
    add_prediction_batch(model, combined_jobs, interval, ontology_terms, combined_callback, batch_size=2, workers=2)
    region_rows.append({
        "region_id": region_number, "cromosoma": chrom,
        "inicio_region_hg38_1based": positions[0], "fin_region_hg38_1based": positions[-1],
        "n_sitios_unicos": len(sites), "n_variantes_AR_loss": len(by_category["AR_loss"]),
        "n_variantes_ESR1_gain": len(by_category["ESR1_gain"]),
        "n_posiciones_GRCh38_distintas_del_alelo_humano_csv": mismatch_count,
        "ancho_contexto_modelo_bp": LONGITUD_SECUENCIA,
        "metodo": "predict_sequence: haplotipo humano versus reversión ancestral del grupo; interacción = conjunto - suma de individuales",
    })


def main():
    if not ARCHIVO_ENTRADA.exists():
        print(f"[ERROR] No se encontró {ARCHIVO_ENTRADA.name} en {CARPETA}")
        return 2
    key = get_api_key()
    if not key:
        print("[ERROR] Falta la API key en api_key.txt o en ALPHAGENOME_API_KEY.")
        return 2
    try:
        model = dna_client.create(key)
        candidates = create_candidates()
        if not candidates:
            print("[ERROR] Ninguna variante pasó los filtros de motivos configurados.")
            return 2
        ontology_terms = output_ontology_terms(model)
    except Exception as error:
        print(f"[ERROR] Inicialización fallida: {type(error).__name__}: {error}")
        return 2

    regions = partition_regions(candidates)
    print(f"Regiones/haplotipos por ventana: {len(regions)}; máximo span de variantes: {MAX_SPAN_HAPLOTIPO:,} bp.")
    print("Predicción: AlphaGenome DNA model, no AVI global de Atlas.")
    print("Métrica clave: efecto conjunto en tracks comparables; interacción = efecto conjunto menos suma de efectos individuales.")

    all_track_rows = []
    all_individual_rows = []
    region_rows = []
    for region_number, region in enumerate(regions, start=1):
        try:
            process_region(model, region, region_number, ontology_terms, all_track_rows, all_individual_rows, region_rows)
        except Exception as error:
            print(f"[ERROR] Región {region_number} falló: {type(error).__name__}: {error}")
            region_rows.append({
                "region_id": region_number,
                "cromosoma": region[0]["cromosoma"] if region else "",
                "inicio_region_hg38_1based": min(x["posicion"] for x in region) if region else "",
                "fin_region_hg38_1based": max(x["posicion"] for x in region) if region else "",
                "estado": "ERROR",
                "detalle_error": f"{type(error).__name__}: {error}",
            })

    summary_rows = summarize_tracks(all_track_rows)
    comparisons = make_comparison(summary_rows)
    modality_summary = summarize_modalities(summary_rows)
    write_csv(ARCHIVO_TRACKS, all_track_rows)
    write_csv(ARCHIVO_INDIVIDUAL, all_individual_rows)
    write_csv(ARCHIVO_COMPARACION, comparisons)
    write_csv(ARCHIVO_MODALIDADES, modality_summary)
    write_csv(ARCHIVO_REGIONES, region_rows)

    print("\n=== ANÁLISIS CONJUNTO TERMINADO ===")
    print(f"Regiones procesadas: {sum(1 for x in region_rows if x.get('estado') != 'ERROR')} / {len(regions)}")
    print(f"Predicciones conjuntas por track: {len(all_track_rows)} filas")
    print(f"Predicciones individuales por track: {len(all_individual_rows)} filas")
    print(f"Tracks comparables entre AR_loss y ESR1_gain: {len(comparisons)}")
    print(f"1. RESULTADO PRINCIPAL: {ARCHIVO_COMPARACION.name}")
    print(f"2. Efectos conjuntos y no aditividad: {ARCHIVO_TRACKS.name}")
    print(f"3. Aportes de variantes individuales: {ARCHIVO_INDIVIDUAL.name}")
    print(f"4. Resumen por modalidad: {ARCHIVO_MODALIDADES.name}")
    print(f"5. Control de regiones y errores: {ARCHIVO_REGIONES.name}")

    direct_rows = [r for r in all_track_rows if r.get("factor_directo_detectado") in {"AR", "ESR1"}]
    direct_matches = [
        r for r in direct_rows
        if (r["categoria"] == "AR_loss" and r["factor_directo_detectado"] == "AR")
        or (r["categoria"] == "ESR1_gain" and r["factor_directo_detectado"] == "ESR1")
    ]
    print(f"Tracks con factor directo AR/ESR1 identificable: {len(direct_rows)}; concordantes con dirección esperada en el grupo/factor: {sum(r.get('direccion_concordante_en_track_factor_directo') is True for r in direct_matches)}")
    if not comparisons:
        print("[AVISO] No hay tracks idénticos compartidos entre las dos categorías para una comparación directa. Revise los CSV por categoría; no se fabricó un ganador global.")
    else:
        print("\n=== COMPARACIÓN DIRECTA POR MODALIDAD (solo tracks idénticos) ===")
        for modality in sorted({row["modalidad"] for row in comparisons}):
            subset = [row for row in comparisons if row["modalidad"] == modality]
            ratios = [float(row["ratio_abs_AR_sobre_ESR1"]) for row in subset if row["ratio_abs_AR_sobre_ESR1"] != ""]
            ar_wins = sum(row["ganador_por_impacto_conjunto"] == "AR_loss" for row in subset)
            er_wins = sum(row["ganador_por_impacto_conjunto"] == "ESR1_gain" for row in subset)
            ties = sum(row["ganador_por_impacto_conjunto"] == "EMPATE" for row in subset)
            median_ratio = float(np.median(ratios)) if ratios else float("nan")
            print(
                f"{modality}: mediana ratio AR_loss/ESR1_gain={median_ratio:.4g}; "
                f"AR_loss mayor en {ar_wins} tracks, ESR1_gain mayor en {er_wins}, empates={ties}."
            )
        print("Interprete cada modalidad por separado; no sume unidades de modalidades distintas en un único número.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
