from collections import Counter, defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import re
import shutil
import ssl
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import certifi

CARPETA = Path(__file__).resolve().parent
BASES = set("ACGT")
ESPECIES_ANCESTRALES = ("chimpance", "gorila", "orangutan")
ESPECIES_SALIDA = ("humano", "chimpance", "gorila", "orangutan", "bonobo")
MIN_DELECION = 10
MIN_E_NO_ALINEADO = 1000
FLANCO_MOTIVOS = 40

# --- UMBRALES EQUILIBRADOS ---
UMBRAL_MOTIVO = 0.70
DELTA_MOTIVO = 0.05
MAX_DISTANCIA_INTERACCION = 30

SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
COMPLEMENTO = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn-", "TGCAYRKMVHDBNtgcayrkmvhdbn-"
)
ARE_FIJAS = {
    0: "A",
    1: "G",
    2: "A",
    3: "A",
    4: "C",
    5: "A",
    9: "T",
    10: "G",
    11: "T",
    12: "T",
    13: "C",
    14: "T",
}
ERE_CONSENSO = "AGGTCANNN TGACCT".replace(" ", "")
UMBRAL_CONSENSO_ERE = 0.65
DELTA_CONSENSO_ERE = 0.05
UMBRAL_CONSENSO_ARE = 0.65
DELTA_CONSENSO_ARE = 0.05

URL_ENHANCERS = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed"
URL_JASPAR_MEME = "https://jaspar2026.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
FILER_BASE = "https://filer2.niagads.org/Annotationtracks/EpiMap_enhancers/ChromHMM/bed9/hg38/formatted_output_{sample}_{source}_Enh_hg38.bed.gz"

GENES = {
    "AR": {
        "chrom": "chrX",
        "files": (
            "alineamiento_ar_ampliado.maf",
            "AR_region_distal.maf",
            "alineamiento_ar.maf",
        ),
        "tokens": ("ar", "androgenreceptor", "receptorandrogenico"),
    },
    "ESR1": {
        "chrom": "chr6",
        "files": ("alineamiento_esr1.maf", "ESR1_region.maf", "esr1.maf"),
        "tokens": (
            "esr1",
            "eralfa",
            "eralpha",
            "estrogenreceptoralpha",
        ),
    },
    "ESR2": {
        "chrom": "chr14",
        "files": ("alineamiento_esr2.maf", "ESR2_region.maf", "esr2.maf"),
        "tokens": (
            "esr2",
            "erbeta",
            "erbeta",
            "estrogenreceptorbeta",
        ),
    },
    "CYP19A1": {
        "chrom": "chr15",
        "files": ("alineamiento_cyp19a1.maf", "CYP19A1_region.maf", "cyp19a1.maf"),
        "tokens": ("cyp19a1", "aromatase", "aromatasa"),
    },
}

EPIMAP_CELL_TRACKS = {
    "luminal_epithelial_mammary": {
        "label": "Mammary luminal epithelial cell (adult; EpiMap BSS01185)",
        "sample": "BSS01185",
        "sources": ["Roadmap_2015", "Roadmap_2021", "ENCODE"],
        "primary_mammary": True,
    },
    "mammary_epithelial_50y": {
        "label": "Mammary gland epithelial cell (adult, 50 y; EpiMap BSS01213)",
        "sample": "BSS01213",
        "sources": ["ENCODE", "ENCODE_2012", "Roadmap_2015"],
        "primary_mammary": True,
    },
    "mammary_epithelial_18y": {
        "label": "Mammary gland epithelial cell (adult, 18 y; EpiMap BSS01211)",
        "sample": "BSS01211",
        "sources": ["ENCODE", "ENCODE_2012", "Roadmap_2015"],
        "primary_mammary": True,
    },
    "mammary_epithelial_other": {
        "label": "Mammary gland epithelial cell (EpiMap BSS01209)",
        "sample": "BSS01209",
        "sources": ["ENCODE", "ENCODE_2012", "Roadmap_2015"],
        "primary_mammary": True,
    },
    "mammary_myoepithelial_33y": {
        "label": "Mammary myoepithelial cell (adult, 33 y; EpiMap BSS01340)",
        "sample": "BSS01340",
        "sources": ["Roadmap_2015", "Roadmap_2021", "ENCODE"],
        "primary_mammary": True,
    },
    "mammary_myoepithelial_36y": {
        "label": "Mammary myoepithelial cell (adult, 36 y; EpiMap BSS01341)",
        "sample": "BSS01341",
        "sources": ["Roadmap_2015", "Roadmap_2021", "ENCODE"],
        "primary_mammary": True,
    },
    "mammary_fibroblast": {
        "label": "Mammary/breast fibroblast (EpiMap BSS00343)",
        "sample": "BSS00343",
        "sources": ["Roadmap_2015", "ENCODE", "Roadmap_2021"],
        "primary_mammary": True,
    },
    "adipocyte_general": {
        "label": (
            "Adipocyte differentiated in vitro; not breast-specific (EpiMap"
            " BSS00039)"
        ),
        "sample": "BSS00039",
        "sources": ["Roadmap_2015", "Roadmap_2021", "ENCODE"],
        "primary_mammary": False,
    },
    "adipose_tissue_general": {
        "label": (
            "Adult adipose tissue; not breast-specific (EpiMap BSS01665)"
        ),
        "sample": "BSS01665",
        "sources": ["Roadmap_2015", "Roadmap_2021", "ENCODE"],
        "primary_mammary": False,
    },
}
OPTIONAL_PROGENITOR_BED = CARPETA / "mammary_progenitor_enhancers.bed"


def descargar(url, destino, comprimido=False):
  destino = Path(destino)
  if destino.exists() and destino.stat().st_size > 1000:
    if not comprimido:
      return destino
    try:
      with gzip.open(destino, "rb") as prueba:
        prueba.read(1)
      return destino
    except Exception:
      destino.unlink()
  temporal = destino.with_suffix(destino.suffix + ".part")
  print(f"Descargando {destino.name}...")
  try:
    with urlopen(url, context=SSL_CONTEXT, timeout=300) as respuesta, open(
        temporal, "wb"
    ) as salida:
      shutil.copyfileobj(respuesta, salida)
    if comprimido:
      with gzip.open(temporal, "rb") as prueba:
        prueba.read(1)
    temporal.replace(destino)
  except Exception:
    if temporal.exists():
      temporal.unlink()
    raise
  return destino


def identificar_especie(src):
  prefijo = src.split(".")[0].lower()
  texto = src.lower()
  if (
      prefijo in {"hg38", "human", "homo_sapiens", "hg19"}
      or "homo_sapiens" in texto
  ):
    return "humano"
  if (
      prefijo.startswith("pantro5")
      or "pan_troglodytes" in texto
      or prefijo in {"chimp", "chimpanzee"}
  ):
    return "chimpance"
  if (
      prefijo.startswith("gorgor5")
      or "gorilla_gorilla" in texto
      or prefijo == "gorilla"
  ):
    return "gorila"
  if (
      prefijo.startswith("ponabe2")
      or "pongo_abelii" in texto
      or prefijo == "orangutan"
  ):
    return "orangutan"
  if (
      prefijo.startswith("panpan2")
      or "pan_paniscus" in texto
      or prefijo == "bonobo"
  ):
    return "bonobo"
  return None


def cromosoma_src(src):
  partes = src.split(".", 1)
  if len(partes) < 2:
    return ""
  return partes[1].split(":", 1)[0]


def complementarseq(seq):
  return seq.translate(COMPLEMENTO)


def reverse_complement(seq):
  return seq.translate(COMPLEMENTO)[::-1]


def resolver_maf(gene):
  config = GENES[gene]
  for nombre in config["files"]:
    ruta = CARPETA / nombre
    if ruta.exists():
      return ruta
  
  maf_files = list(CARPETA.glob("*.maf"))
  for ruta in maf_files:
    stem_lower = ruta.stem.lower()
    if gene.lower() in stem_lower:
      return ruta
      
  tokens = set(config["tokens"])
  for ruta in maf_files:
    normalizado = re.sub(r"[^a-z0-9]+", "_", ruta.stem.lower()).strip("_")
    partes = set(normalizado.split("_"))
    if partes.intersection(tokens):
      return ruta

  if maf_files:
    mapa_orden = {"AR": 0, "ESR1": 1, "ESR2": 2, "CYP19A1": 3}
    idx = mapa_orden.get(gene, 0)
    if idx < len(maf_files):
      return maf_files[idx]
    return maf_files[0]

  raise FileNotFoundError(
      f"No encontré ningún archivo .maf válido para el gen {gene}."
  )


def iter_maf_bloques(ruta):
  bloque = []
  with open(ruta, encoding="utf-8", errors="replace") as archivo:
    for linea in archivo:
      linea = linea.rstrip("\r\n")
      if not linea.strip():
        if bloque:
          yield bloque
          bloque = []
      elif not linea.startswith("#") and not linea.startswith("track "):
        bloque.append(linea)
  if bloque:
    yield bloque


def parsear_bloque(lineas, cromosoma_esperado):
  filas_s, filas_e = [], []
  for linea in lineas:
    campos = linea.split()
    if campos and campos[0] == "s" and len(campos) >= 7:
      filas_s.append({
          "src": campos[1],
          "chrom": cromosoma_src(campos[1]),
          "start": int(campos[2]),
          "size": int(campos[3]),
          "strand": campos[4],
          "src_size": int(campos[5]),
          "text": campos[6].upper(),
          "especie": identificar_especie(campos[1]),
      })
    elif campos and campos[0] == "e" and len(campos) >= 7:
      filas_e.append({
          "src": campos[1],
          "chrom": cromosoma_src(campos[1]),
          "start": int(campos[2]),
          "size": int(campos[3]),
          "strand": campos[4],
          "src_size": int(campos[5]),
          "status": campos[6],
          "especie": identificar_especie(campos[1]),
      })
  if not filas_s:
    return None
  longitud = len(filas_s[0]["text"])
  if any(len(f["text"]) != longitud for f in filas_s):
    return None
  humanos = [
      f
      for f in filas_s
      if f["especie"] == "humano"
      and f["src"].split(".")[0].lower() == "hg38"
      and f["chrom"].lower() == cromosoma_esperado.lower()
  ]
  if len(humanos) != 1:
    return None
  humano = humanos[0]
  filas = {}
  for fila in filas_s:
    especie = fila["especie"]
    if especie is None:
      continue
    if especie not in filas or fila["size"] > filas[especie]["size"]:
      filas[especie] = fila
  coords = [None] * longitud
  contador = 0
  for i, base in enumerate(humano["text"]):
    if base == "-":
      continue
    if humano["strand"] == "+":
      coords[i] = humano["start"] + contador
    else:
      coords[i] = humano["src_size"] - humano["start"] - 1 - contador
    contador += 1
  for fila in filas.values():
    fila["text_hg38"] = (
        complementarseq(fila["text"])
        if humano["strand"] == "-"
        else fila["text"]
    )
  coords_validas = [c for c in coords if c is not None]
  if not coords_validas:
    return None
  return {
      "humano": humano,
      "filas": filas,
      "coords": coords,
      "e_rows": filas_e,
      "inicio": min(coords_validas),
      "fin": max(coords_validas) + 1,
      "columnas": longitud,
  }


def guardar_csv(ruta, filas, campos=None):
  ruta = Path(ruta)
  if campos is None:
    campos = (
        list(dict.fromkeys(k for fila in filas for k in fila)) if filas else []
    )
  if not campos:
    campos = ["sin_registros"]
  with open(ruta, "w", newline="", encoding="utf-8-sig") as archivo:
    escritor = csv.DictWriter(archivo, fieldnames=campos, extrasaction="ignore")
    escritor.writeheader()
    if filas:
      escritor.writerows(filas)


def extraer_variantes_y_deleciones(gene, ruta_maf):
  meta = GENES[gene]
  sustituciones, deleciones, segmentos_e = {}, [], []
  fuentes = Counter()
  n_bloques, validos = 0, 0
  limites_maf = []
  claves_e = set()
  claves_delecion = set()
  for lineas in iter_maf_bloques(ruta_maf):
    n_bloques += 1
    for linea in lineas:
      campos = linea.split()
      if campos and campos[0] == "s" and len(campos) >= 7:
        fuentes[campos[1].split(".")[0]] += 1
    bloque = parsear_bloque(lineas, meta["chrom"])
    if bloque is None:
      continue
    validos += 1
    limites_maf.append((bloque["inicio"], bloque["fin"]))
    filas, coords, hrow = (
        bloque["filas"],
        bloque["coords"],
        bloque["filas"]["humano"],
    )
    hseq = hrow["text_hg38"]
    longitud = bloque["columnas"]
    if all(especie in filas for especie in ESPECIES_ANCESTRALES):
      for i, coord in enumerate(coords):
        if coord is None or hseq[i] not in BASES:
          continue
        simios = {
            especie: filas[especie]["text_hg38"][i]
            for especie in ESPECIES_ANCESTRALES
        }
        if (
            not all(base in BASES for base in simios.values())
            or len(set(simios.values())) != 1
        ):
          continue
        ancestral, humano = next(iter(simios.values())), hseq[i]
        if humano == ancestral:
          continue
        bonobo = filas.get("bonobo", {}).get("text_hg38", "")
        bonobo = bonobo[i] if len(bonobo) > i and bonobo[i] in BASES else ""
        coord1 = coord + 1
        dato = {
            "gen_locus": gene,
            "cromosoma": meta["chrom"],
            "coordenada_hg38_1based": coord1,
            "base_ancestral_inferida": ancestral,
            "base_humana_maf": humano,
            "base_chimpance": simios["chimpance"],
            "base_gorila": simios["gorila"],
            "base_orangutan": simios["orangutan"],
            "base_bonobo": bonobo,
            "bonobo_apoya_base_ancestral": bonobo == ancestral
            if bonobo
            else "",
            "bloque_maf": n_bloques,
            "columna_bloque_1based": i + 1,
            "estado_ancestral": (
                "consenso_chimpance_gorila_orangutan; no es reconstruccion"
                " filogenetica formal"
            ),
        }
        if coord1 not in sustituciones:
          sustituciones[coord1] = dato
        else:
          anterior = sustituciones[coord1]
          anterior["repeticiones_en_MAF"] = (
              int(anterior.get("repeticiones_en_MAF", 1)) + 1
          )
          if (
              anterior["base_ancestral_inferida"] != ancestral
              or anterior["base_humana_maf"] != humano
          ):
            anterior["conflicto_entre_bloques"] = True
    i = 0
    while i < longitud:
      soporte = sum(
          especie in filas and filas[especie]["text_hg38"][i] in BASES
          for especie in ESPECIES_ANCESTRALES
      )
      if hrow["text"][i] != "-" or soporte < 2:
        i += 1
        continue
      inicio_col = i
      while i < longitud:
        soporte_i = sum(
            especie in filas and filas[especie]["text_hg38"][i] in BASES
            for especie in ESPECIES_ANCESTRALES
        )
        if hrow["text"][i] != "-" or soporte_i < 2:
          break
        i += 1
      fin_col = i
      if fin_col - inicio_col < MIN_DELECION:
        continue
      izquierda = next(
          (j for j in range(inicio_col - 1, -1, -1) if coords[j] is not None),
          None,
      )
      derecha = next(
          (j for j in range(fin_col, longitud) if coords[j] is not None), None
      )
      if izquierda is None or derecha is None:
        continue
      ci, cd = coords[izquierda], coords[derecha]
      columnas_soporte = [
          sum(
              especie in filas and filas[especie]["text_hg38"][j] in BASES
              for especie in ESPECIES_ANCESTRALES
          )
          for j in range(inicio_col, fin_col)
      ]
      punto_union = max(ci, cd)
      clave_del = (punto_union, fin_col - inicio_col)
      if clave_del in claves_delecion:
        continue
      claves_delecion.add(clave_del)

      def secuencia_especie(esp):
        if esp not in filas:
          return ""
        s = filas[esp]["text_hg38"][inicio_col:fin_col]
        return s[::-1] if hrow["strand"] == "-" else s

      deleciones.append({
          "gen_locus": gene,
          "cromosoma": meta["chrom"],
          "bloque_maf": n_bloques,
          "columna_inicio_bloque_1based": inicio_col + 1,
          "columna_fin_bloque_1based": fin_col,
          "longitud_columnas": fin_col - inicio_col,
          "flanco_hg38_menor_1based": min(ci, cd) + 1,
          "flanco_hg38_mayor_1based": max(ci, cd) + 1,
          "flancos_humanos_consecutivos": abs(ci - cd) == 1,
          "punto_union_interbase_hg38_0based": punto_union,
          "soporte_minimo_simios": min(columnas_soporte),
          "secuencia_chimpance": secuencia_especie("chimpance"),
          "secuencia_gorila": secuencia_especie("gorila"),
          "secuencia_orangutan": secuencia_especie("orangutan"),
          "secuencia_bonobo": secuencia_especie("bonobo"),
          "estado": (
              "candidato_por_hueco_MAF; validar ortologia y alineamiento"
          ),
      })
    if bloque["fin"] > bloque["inicio"]:
      for erow in bloque["e_rows"]:
        esp = erow["especie"]
        if esp not in ESPECIES_SALIDA or erow["size"] < MIN_E_NO_ALINEADO:
          continue
        e_start = (
            erow["start"]
            if erow["strand"] == "+"
            else erow["src_size"] - erow["start"] - erow["size"]
        )
        e_end = e_start + erow["size"]
        key = (
            erow["src"],
            e_start,
            e_end,
            erow["status"],
            bloque["inicio"],
            bloque["fin"],
        )
        if key in claves_e:
          continue
        claves_e.add(key)
        segmentos_e.append({
            "gen_locus": gene,
            "especie": esp,
            "src": erow["src"],
            "inicio_especie_1based": e_start + 1,
            "fin_especie_1based": e_end,
            "longitud": erow["size"],
            "estado_MAF": erow["status"],
            "bloque_humano_inicio_hg38_1based": bloque["inicio"] + 1,
            "bloque_humano_fin_hg38_1based": bloque["fin"],
        })
  variantes = [sustituciones[k] for k in sorted(sustituciones)]
  if not limites_maf:
    raise RuntimeError(
        f"No encontré bloques MAF con la referencia humana {meta['chrom']} en"
        f" {ruta_maf.name}. Verifica el ensamblaje y los nombres src del MAF."
    )
  region = (
      max(0, min(x[0] for x in limites_maf) - 1000),
      max(x[1] for x in limites_maf) + 1000,
  )
  print(
      f"{gene}: {ruta_maf.name} | bloques leídos {n_bloques:,} | bloques útiles"
      f" {validos:,} | sustituciones {len(variantes):,} | huecos candidatos"
      f" {len(deleciones):,} | región aproximada"
      f" {meta['chrom']}:{region[0]}-{region[1]}"
  )
  return variantes, deleciones, segmentos_e, region, dict(fuentes)


def inferir_clase_ccre(raw, rgb):
  texto = raw.lower()
  for token in (
      "dels",
      "pels",
      "pls",
      "ca-ctcf",
      "ca-h3k4me3",
      "ca-tf",
      "ctcf",
  ):
    if token in texto:
      return token.upper()
  return {
      "255,0,0": "PLS",
      "255,167,0": "pELS",
      "255,205,0": "dELS",
      "0,176,240": "CA-CTCF",
      "166,166,166": "CA",
      "192,192,192": "CA",
  }.get(rgb, "clase_no_identificada")


def leer_bed_region(
    ruta,
    cromosoma,
    inicio,
    fin,
    origen,
    comprimido=False,
    clase_forzada=None,
):
  intervalos = []
  opener = gzip.open if comprimido else open
  modo = "rt" if comprimido else "r"
  with opener(ruta, modo, encoding="utf-8", errors="replace") as archivo:
    for linea in archivo:
      if not linea.strip() or linea.startswith(("#", "track", "browser")):
        continue
      campos = linea.rstrip("\r\n").split("\t")
      if len(campos) < 3:
        campos = linea.split()
      if len(campos) < 3 or campos[0] != cromosoma:
        continue
      try:
        ini, fin_bed = int(campos[1]), int(campos[2])
      except ValueError:
        continue
      if fin_bed <= inicio or ini >= fin:
        continue
      raw = " ".join(campos[3:])
      rgb = campos[8] if len(campos) > 8 else ""
      clase = clase_forzada or inferir_clase_ccre(raw, rgb)
      intervalos.append({
          "inicio": ini,
          "fin": fin_bed,
          "id": campos[3]
          if len(campos) > 3
          else f"{cromosoma}:{ini}-{fin_bed}",
          "clase": clase,
          "rgb": rgb,
          "origen": origen,
          "raw": raw,
      })
  intervalos.sort(key=lambda x: x["inicio"])
  return intervalos


def descargar_anotaciones(gene, cromosoma, inicio, fin):
  archivo_ccre = CARPETA / "GRCh38_cCREs_ELS_V4.bed"
  descargar(URL_ENHANCERS, archivo_ccre)
  enhancers = leer_bed_region(
      archivo_ccre, cromosoma, inicio, fin, "SCREEN Registry V4 ELS cCRE"
  )
  tracks, estado = {}, []
  for clave, config in EPIMAP_CELL_TRACKS.items():
    ruta_ok, fuente, errores = None, "", []
    for source in config["sources"]:
      ruta = CARPETA / f"EpiMap_{config['sample']}_{source}_Enh_hg38.bed.gz"
      url = FILER_BASE.format(sample=config["sample"], source=source)
      try:
        descargar(url, ruta, comprimido=True)
        ruta_ok, fuente = ruta, source
        break
      except Exception as error:
        errores.append(f"{source}: {error}")
    if ruta_ok is None:
      estado.append({
          "gen_locus": gene,
          "grupo_celular": clave,
          "muestra": config["sample"],
          "etiqueta": config["label"],
          "estado": "no_disponible",
          "fuente": "",
          "archivo": "",
          "n_intervalos_region": 0,
          "detalle": " | ".join(errores[-2:]),
      })
      print(f"[SIN DATOS] {gene} / {clave}")
      continue
    intervalos = leer_bed_region(
        ruta_ok,
        cromosoma,
        inicio,
        fin,
        f"EpiMap: {config['label']}",
        comprimido=True,
        clase_forzada="EpiMap_enhancer_chromatin_state",
    )
    tracks[clave] = {
        "label": config["label"],
        "sample": config["sample"],
        "source": fuente,
        "primary_mammary": config["primary_mammary"],
        "intervalos": intervalos,
        "file": ruta_ok.name,
    }
    estado.append({
        "gen_locus": gene,
        "grupo_celular": clave,
        "muestra": config["sample"],
        "etiqueta": config["label"],
        "estado": "descargado",
        "fuente": fuente,
        "archivo": ruta_ok.name,
        "n_intervalos_region": len(intervalos),
        "detalle": "ChromHMM enhancers; hg38",
    })
  if OPTIONAL_PROGENITOR_BED.exists():
    custom = leer_bed_region(
        OPTIONAL_PROGENITOR_BED,
        cromosoma,
        inicio,
        fin,
        "BED local de usuario",
        clase_forzada="progenitor_mamario_enhancer_usuario",
    )
    tracks["mammary_progenitor_custom"] = {
        "label": "Mammary stem/progenitor (custom BED supplied by user)",
        "sample": "custom",
        "source": "local BED",
        "primary_mammary": True,
        "intervalos": custom,
        "file": OPTIONAL_PROGENITOR_BED.name,
    }
    estado.append({
        "gen_locus": gene,
        "grupo_celular": "mammary_progenitor_custom",
        "muestra": "custom",
        "etiqueta": "Progenitor mamario; BED local",
        "estado": "descargado",
        "fuente": "local BED",
        "archivo": OPTIONAL_PROGENITOR_BED.name,
        "n_intervalos_region": len(custom),
        "detalle": "Validez dependiente del BED suministrado",
    })
  guardar_csv(CARPETA / f"estado_anotaciones_{gene}.csv", estado)
  return enhancers, tracks, estado


def construir_indice(intervalos):
  starts, maximos, max_end = [], [], -1
  for item in intervalos:
    starts.append(item["inicio"])
    max_end = max(max_end, item["fin"])
    maximos.append(max_end)
  return starts, maximos


def solapamientos(posicion, intervalos, indice):
  from bisect import bisect_right

  starts, maximos = indice
  i = bisect_right(starts, posicion) - 1
  hits = []
  while i >= 0 and maximos[i] > posicion:
    x = intervalos[i]
    if x["inicio"] <= posicion < x["fin"]:
      hits.append(x)
    i -= 1
  return list(reversed(hits))


def distancia_intervalo(posicion, intervalo):
  if intervalo["inicio"] <= posicion < intervalo["fin"]:
    return 0
  if posicion < intervalo["inicio"]:
    return intervalo["inicio"] - posicion
  return posicion - intervalo["fin"] + 1


def descargar_secuencia_hg38(gene, cromosoma, inicio, fin):
  ruta = CARPETA / f"hg38_{gene}_{cromosoma}_{inicio}_{fin}.fasta"
  if ruta.exists():
    lineas = ruta.read_text(encoding="utf-8").splitlines()
    seq = "".join(x.strip() for x in lineas if not x.startswith(">")).upper()
    if len(seq) == fin - inicio and set(seq) <= BASES | {"N"}:
      return seq
  partes = []
  tamano = 200_000
  for a in range(inicio, fin, tamano):
    b = min(fin, a + tamano)
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={cromosoma};start={a};end={b}"
    with urlopen(url, context=SSL_CONTEXT, timeout=180) as respuesta:
      datos = json.load(respuesta)
    if "dna" not in datos:
      raise RuntimeError(
          f"UCSC no devolvió ADN para {cromosoma}:{a}-{b}: {datos}"
      )
    partes.append(datos["dna"].upper())
    print(
        f"Secuencia hg38 {gene}:"
        f" {min(b - inicio, fin - inicio):,}/{fin - inicio:,} bp",
        end="\r",
    )
  seq = "".join(partes)
  print()
  if len(seq) != fin - inicio:
    raise RuntimeError(
        f"Longitud inesperada de hg38 para {gene}: {len(seq)}; esperada"
        f" {fin-inicio}"
    )
  ruta.write_text(
      f">hg38_{gene}_{cromosoma}_{inicio}_{fin}\n"
      + "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80))
      + "\n",
      encoding="utf-8",
  )
  return seq


def cargar_motivos_jaspar():
  archivo = CARPETA / "JASPAR2026_CORE_vertebrates_nonredundant.meme.txt"
  descargar(URL_JASPAR_MEME, archivo)
  lineas = archivo.read_text(encoding="utf-8", errors="replace").splitlines()
  motivos, i = [], 0
  while i < len(lineas):
    if not lineas[i].startswith("MOTIF "):
      i += 1
      continue
    partes = lineas[i].split(maxsplit=2)
    motif_id = partes[1] if len(partes) > 1 else "unknown"
    nombre = partes[2] if len(partes) > 2 else motif_id
    j = i + 1
    while j < len(lineas) and not lineas[j].startswith(
        "letter-probability matrix:"
    ):
      if lineas[j].startswith("MOTIF "):
        break
      j += 1
    if j >= len(lineas) or not lineas[j].startswith(
        "letter-probability matrix:"
    ):
      i += 1
      continue
    match = re.search(r"\bw\s*=\s*(\d+)", lineas[j])
    if not match:
      i = j + 1
      continue
    width = int(match.group(1))
    if width < 5 or width > 30:
      i = j + 1 + width
      continue
    matrix, ok = [], True
    for fila in lineas[j + 1 : j + 1 + width]:
      try:
        values = [float(x) for x in fila.split()]
      except ValueError:
        ok = False
        break
      if len(values) != 4:
        ok = False
        break
      matrix.append(values)
    if not ok or len(matrix) != width:
      i = j + 1
      continue
    log_matrix, min_score, max_score = [], 0.0, 0.0
    for values in matrix:
      total = sum(values)
      probs = [(max(0.0, v) + 0.001) / (total + 0.004) for v in values]
      row = {base: math.log2(probs[k] / 0.25) for k, base in enumerate("ACGT")}
      log_matrix.append(row)
      min_score += min(row.values())
      max_score += max(row.values())
    if max_score > min_score:
      normalized_name = re.sub(
          r"[^a-z0-9]+", "", (nombre + " " + motif_id).lower()
      )
      base_id = motif_id.split(".")[0]
      is_ar = (
          base_id == "MA0007"
          or bool(re.search(r"(^|[^a-z])ar([^a-z]|$)", nombre.lower()))
          or "androgenreceptor" in normalized_name
      )
      is_esr1 = (
          base_id == "MA0112"
          or "esr1" in normalized_name
          or "estrogenreceptoralpha" in normalized_name
          or "estrogenreceptoralfa" in normalized_name
      )
      is_esr2 = (
          "esr2" in normalized_name
          or "estrogenreceptorbeta" in normalized_name
          or "estrogenreceptorbeta" in normalized_name
      )
      motivos.append({
          "id": motif_id,
          "name": nombre,
          "width": width,
          "matrix": log_matrix,
          "min_score": min_score,
          "max_score": max_score,
          "is_ar": is_ar,
          "is_esr1": is_esr1,
          "is_esr2": is_esr2,
          "is_ere_family": is_esr1 or is_esr2,
      })
    i = j + 1 + width
  if not motivos:
    raise RuntimeError("No se pudieron leer matrices JASPAR del archivo MEME.")
  print(
      f"Matrices JASPAR cargadas: {len(motivos):,}; AR:"
      f" {sum(m['is_ar'] for m in motivos)}; ESR1:"
      f" {sum(m['is_esr1'] for m in motivos)}; ESR2 específicas:"
      f" {sum(m['is_esr2'] for m in motivos)}"
  )
  return motivos


def puntuacion_pwm(ventana, motivo):
  if len(ventana) != motivo["width"] or not all(base in BASES for base in ventana):
    return None
  score = sum(motivo["matrix"][i][base] for i, base in enumerate(ventana))
  return (score - motivo["min_score"]) / (
      motivo["max_score"] - motivo["min_score"]
  )


def mejor_par_pwm_centrada(seq_anc, seq_hum, centro, motivo):
  width = motivo["width"]
  lo, hi = max(0, centro - width + 1), min(centro, len(seq_anc) - width)
  best = None
  for start in range(lo, hi + 1):
    va, vh = seq_anc[start : start + width], seq_hum[start : start + width]
    for orientacion in ("+", "-"):
      aa, hh = (
          (reverse_complement(va), reverse_complement(vh))
          if orientacion == "-"
          else (va, vh)
      )
      sa, sh = puntuacion_pwm(aa, motivo), puntuacion_pwm(hh, motivo)
      if sa is None or sh is None:
        continue
      if best is None or (abs(sa - sh), max(sa, sh)) > (
          abs(best["score_ancestral"] - best["score_humano"]),
          max(best["score_ancestral"], best["score_humano"]),
      ):
        best = {
            "score_ancestral": sa,
            "score_humano": sh,
            "start": start,
            "orientacion": orientacion,
            "secuencia_ancestral": aa,
            "secuencia_humana": hh,
        }
  return best


def consenso_score(ventana, tipo):
  if len(ventana) != 15 or not all(base in BASES for base in ventana):
    return None
  if tipo == "ARE":
    return (
        sum(ventana[p] == b for p, b in ARE_FIJAS.items()) / len(ARE_FIJAS)
    )
  posiciones = [(i, base) for i, base in enumerate(ERE_CONSENSO) if base != "N"]
  return sum(ventana[i] == base for i, base in posiciones) / len(posiciones)


def mejor_par_consenso(seq_anc, seq_hum, centro, tipo):
  width = 15
  lo, hi = max(0, centro - width + 1), min(centro, len(seq_anc) - width)
  best = None
  for start in range(lo, hi + 1):
    va, vh = seq_anc[start : start + width], seq_hum[start : start + width]
    for orientacion in ("+", "-"):
      aa, hh = (
          (reverse_complement(va), reverse_complement(vh))
          if orientacion == "-"
          else (va, vh)
      )
      sa, sh = consenso_score(aa, tipo), consenso_score(hh, tipo)
      if sa is None or sh is None:
        continue
      if best is None or (abs(sa - sh), max(sa, sh)) > (
          abs(best["score_ancestral"] - best["score_humano"]),
          max(best["score_ancestral"], best["score_humano"]),
      ):
        best = {
            "score_ancestral": sa,
            "score_humano": sh,
            "start": start,
            "orientacion": orientacion,
            "secuencia_ancestral": aa,
            "secuencia_humana": hh,
        }
  return best


def clasificar_factor(nombre):
  n = nombre.upper().replace(" ", "")
  if n in {"AR", "ANDROGENRECEPTOR"}:
    return "AR; contexto funcional no inferible del motivo"
  if "ESR1" in n or "ESTROGENRECEPTORALPHA" in n:
    return "ESR1/ERalpha; motivo candidato"
  if "ESR2" in n or "ESTROGENRECEPTORBETA" in n:
    return "ESR2/ERbeta; motivo candidato"
  if any(x in n for x in ("FOXA1", "FOXA2", "FOXA3")):
    return "FOXA; factor pionero/candidato"
  if "GATA3" in n or "GATA2" in n:
    return "GATA; linaje epitelial/candidato"
  if "TFAP2C" in n or "TFAP2A" in n:
    return "TFAP2; factor epitelial/candidato"
  if "GRHL2" in n or "GRHL3" in n:
    return "GRHL; factor epitelial/candidato"
  if any(x in n for x in ("JUN", "FOS", "AP1")):
    return "AP-1; candidato"
  if "TEAD" in n:
    return "TEAD; candidato"
  if "STAT5" in n or "ELF5" in n or "PGR" in n:
    return "desarrollo mamario/hormonal; candidato"
  return "otro TF; función a verificar"


def puntuar_variante(variante, genoma, inicio_region, motivos):
  coord1 = int(variante["coordenada_hg38_1based"])
  pos0, idx = coord1 - 1, coord1 - 1 - inicio_region
  if idx < 0 or idx >= len(genoma):
    return {"estado_modelado": "fuera_de_secuencia_hg38"}, []
  humano, ancestral = genoma[idx], variante["base_ancestral_inferida"]
  if humano not in BASES or ancestral not in BASES:
    return {"estado_modelado": "base_no_ACGT"}, []
  if humano != variante["base_humana_maf"]:
    return {"estado_modelado": "discordancia_MAF_hg38", "base_hg38": humano}, []
  a, b = max(inicio_region, pos0 - FLANCO_MOTIVOS), min(
      inicio_region + len(genoma), pos0 + FLANCO_MOTIVOS + 1
  )
  human_seq = genoma[a - inicio_region : b - inicio_region]
  center = pos0 - a
  ancestral_seq = (
      human_seq[:center] + ancestral + human_seq[center + 1 :]
  )
  resultados_motivo = []
  for motivo in motivos:
    pair = mejor_par_pwm_centrada(ancestral_seq, human_seq, center, motivo)
    if pair is None:
      continue
    sa, sh = pair["score_ancestral"], pair["score_humano"]
    delta = sa - sh
    if sa >= UMBRAL_MOTIVO and delta >= DELTA_MOTIVO:
      efecto = "posible_debilitamiento_motivo"
    elif sh >= UMBRAL_MOTIVO and delta <= -DELTA_MOTIVO:
      efecto = "posible_ganancia_motivo_humano"
    else:
      continue
    resultados_motivo.append({
        "gen_locus": variante["gen_locus"],
        "cromosoma": variante["cromosoma"],
        "coordenada_hg38_1based": coord1,
        "base_ancestral": ancestral,
        "base_humana": humano,
        "motif_id": motivo["id"],
        "TF_motif": motivo["name"],
        "es_motivo_AR": motivo["is_ar"],
        "es_motivo_ESR1": motivo["is_esr1"],
        "es_motivo_ESR2_especifico": motivo["is_esr2"],
        "motivo_familia_ERE": motivo["is_ere_family"],
        "grupo_funcional": clasificar_factor(motivo["name"]),
        "score_ancestral": round(sa, 4),
        "score_humano": round(sh, 4),
        "delta_ancestral_menos_humano": round(delta, 4),
        "efecto_predicho": efecto,
        "ventana_inicio_hg38_1based": a + pair["start"] + 1,
        "orientacion_motivo": pair["orientacion"],
        "secuencia_ancestral": pair["secuencia_ancestral"],
        "secuencia_humana": pair["secuencia_humana"],
    })
  are = mejor_par_consenso(ancestral_seq, human_seq, center, "ARE")
  ere = mejor_par_consenso(ancestral_seq, human_seq, center, "ERE")

  def registrar_consenso(pair, tipo):
    if pair is None:
      return
    score_anc = pair["score_ancestral"]
    score_hum = pair["score_humano"]
    delta = score_anc - score_hum
    if tipo == "ERE":
      es_perdida = (
          score_anc >= UMBRAL_CONSENSO_ERE and delta >= DELTA_CONSENSO_ERE
      )
      es_ganancia = (
          score_hum >= UMBRAL_CONSENSO_ERE and delta <= -DELTA_CONSENSO_ERE
      )
      motif_id = "ERE_CONSENSUS_SHARED_ER_ALPHA_BETA"
      nombre = "Shared ERalpha/ERbeta estrogen response element consensus"
      es_ar = False
      es_esr1 = False
      es_esr2 = False
      familia_ere = True
    else:
      es_perdida = (
          score_anc >= UMBRAL_CONSENSO_ARE and delta >= DELTA_CONSENSO_ARE
      )
      es_ganancia = (
          score_hum >= UMBRAL_CONSENSO_ARE and delta <= -DELTA_CONSENSO_ARE
      )
      motif_id = "ARE_CONSENSUS_AR"
      nombre = "Androgen response element consensus"
      es_ar = True
      es_esr1 = False
      es_esr2 = False
      familia_ere = False
    if not es_perdida and not es_ganancia:
      return
    resultados_motivo.append({
        "gen_locus": variante["gen_locus"],
        "tipo_evento": "SNP_candidato",
        "cromosoma": variante["cromosoma"],
        "coordenada_hg38_1based": coord1,
        "base_ancestral": ancestral,
        "base_humana": humano,
        "motif_id": motif_id,
        "TF_motif": nombre,
        "es_motivo_AR": es_ar,
        "es_motivo_ESR1": es_esr1,
        "es_motivo_ESR2_especifico": es_esr2,
        "motivo_familia_ERE": familia_ere,
        "ERE_compartido_ERalpha_ERbeta": tipo == "ERE",
        "grupo_funcional": (
            "ERE compartido por ERalpha y ERbeta"
            if tipo == "ERE"
            else "ARE; unión candidata de AR"
        ),
        "score_ancestral": round(score_anc, 4),
        "score_humano": round(score_hum, 4),
        "delta_ancestral_menos_humano": round(delta, 4),
        "efecto_predicho": (
            "posible_debilitamiento_motivo"
            if es_perdida
            else "posible_ganancia_motivo_humano"
        ),
        "ventana_inicio_hg38_1based": a + pair["start"] + 1,
        "orientacion_motivo": pair["orientacion"],
        "secuencia_ancestral": pair["secuencia_ancestral"],
        "secuencia_humana": pair["secuencia_humana"],
        "nota": (
            "Consenso de secuencia exploratorio; no demuestra unión in vivo,"
            " activación transcripcional ni dirección del efecto. ERE puede"
            " ser reconocido por ERalpha y ERbeta."
        ),
    })

  registrar_consenso(are, "ARE")
  registrar_consenso(ere, "ERE")
  focus_pairs = {}
  for m in motivos:
    if m["is_ar"] or m["is_esr1"] or m["is_esr2"]:
      p = mejor_par_pwm_centrada(ancestral_seq, human_seq, center, m)
      if p is not None:
        focus_pairs[m["id"]] = (p, m)

  def best_flag(flag):
    choices = [(p, m) for p, m in focus_pairs.values() if m[flag]]
    if not choices:
      return {}
    p, m = max(choices, key=lambda x: x[0]["score_ancestral"])
    return {
        "id": m["id"],
        "name": m["name"],
        "anc": p["score_ancestral"],
        "hum": p["score_humano"],
        "delta": p["score_ancestral"] - p["score_humano"],
        "start": a + p["start"] + 1,
    }

  ar, esr1, esr2 = best_flag("is_ar"), best_flag("is_esr1"), best_flag("is_esr2")
  res = {
      "estado_modelado": "comparado_con_hg38; efecto PWM teórico",
      "base_hg38": humano,
      "ARE_consenso_ancestral": round(are["score_ancestral"], 4)
      if are
      else "",
      "ARE_consenso_humano": round(are["score_humano"], 4) if are else "",
      "ARE_delta_ancestral_menos_humano": round(
          are["score_ancestral"] - are["score_humano"], 4
      )
      if are
      else "",
      "ERE_consenso_ancestral": round(ere["score_ancestral"], 4)
      if ere
      else "",
      "ERE_consenso_humano": round(ere["score_humano"], 4) if ere else "",
      "ERE_delta_ancestral_menos_humano": round(
          ere["score_ancestral"] - ere["score_humano"], 4
      )
      if ere
      else "",
      "ERE_consenso_inicio_hg38_1based": a + ere["start"] + 1 if ere else "",
      "ERE_consenso_orientacion": ere["orientacion"] if ere else "",
      "ERE_consenso_secuencia_ancestral": ere["secuencia_ancestral"]
      if ere
      else "",
      "ERE_consenso_secuencia_humana": ere["secuencia_humana"] if ere else "",
      "AR_motif_id": ar.get("id", ""),
      "AR_PWM_ancestral": round(ar["anc"], 4) if ar else "",
      "AR_PWM_humano": round(ar["hum"], 4) if ar else "",
      "AR_PWM_delta_ancestral_menos_humano": round(ar["delta"], 4)
      if ar
      else "",
      "ESR1_motif_id": esr1.get("id", ""),
      "ESR1_PWM_ancestral": round(esr1["anc"], 4) if esr1 else "",
      "ESR1_PWM_humano": round(esr1["hum"], 4) if esr1 else "",
      "ESR1_PWM_delta_ancestral_menos_humano": round(esr1["delta"], 4)
      if esr1
      else "",
      "ESR2_motif_id": esr2.get("id", ""),
      "ESR2_PWM_ancestral": round(esr2["anc"], 4) if esr2 else "",
      "ESR2_PWM_humano": round(esr2["hum"], 4) if esr2 else "",
      "ESR2_PWM_delta_ancestral_menos_humano": round(esr2["delta"], 4)
      if esr2
      else "",
      "n_perdidas_motivos": sum(
          x["efecto_predicho"] == "posible_debilitamiento_motivo"
          for x in resultados_motivo
      ),
      "n_ganancias_motivos": sum(
          x["efecto_predicho"] == "posible_ganancia_motivo_humano"
          for x in resultados_motivo
      ),
      "perdidas_motivos_top": "; ".join(
          f"{x['TF_motif']} {x['delta_ancestral_menos_humano']:+.3f}"
          for x in sorted(
              (
                  x
                  for x in resultados_motivo
                  if x["efecto_predicho"] == "posible_debilitamiento_motivo"
              ),
              key=lambda y: y["delta_ancestral_menos_humano"],
              reverse=True,
          )[:5]
      ),
      "ganancias_motivos_top": "; ".join(
          f"{x['TF_motif']} {x['delta_ancestral_menos_humano']:+.3f}"
          for x in sorted(
              (
                  x
                  for x in resultados_motivo
                  if x["efecto_predicho"] == "posible_ganancia_motivo_humano"
              ),
              key=lambda y: y["delta_ancestral_menos_humano"],
          )[:5]
      ),
  }
  return res, resultados_motivo


def escanear_interacciones(
    gene, variantes, genoma, inicio_region, motivos, posiciones_regulatorias
):
  filas = []
  if not variantes:
    return filas
  focales = [m for m in motivos if m["is_ar"] or m["is_esr1"] or m["is_esr2"]]
  ordenadas = sorted(variantes, key=lambda x: int(x["coordenada_hg38_1based"]))
  for i, v1 in enumerate(ordenadas):
    p1 = int(v1["coordenada_hg38_1based"]) - 1
    for v2 in ordenadas[i + 1 :]:
      p2 = int(v2["coordenada_hg38_1based"]) - 1
      distancia = p2 - p1
      if distancia > MAX_DISTANCIA_INTERACCION:
        break
      if not (
          posiciones_regulatorias.get(p1, False)
          or posiciones_regulatorias.get(p2, False)
      ):
        continue
      lo = max(inicio_region, p1 - FLANCO_MOTIVOS)
      hi = min(inicio_region + len(genoma), p2 + FLANCO_MOTIVOS + 1)
      human_seq = genoma[lo - inicio_region : hi - inicio_region]
      idx1, idx2 = p1 - lo, p2 - lo
      anc1, anc2 = (
          v1["base_ancestral_inferida"],
          v2["base_ancestral_inferida"],
      )
      if human_seq[idx1] not in BASES or human_seq[idx2] not in BASES:
        continue
      if (
          human_seq[idx1] != v1["base_humana_maf"]
          or human_seq[idx2] != v2["base_humana_maf"]
      ):
        continue
      seq1 = human_seq[:idx1] + anc1 + human_seq[idx1 + 1 :]
      seq2 = human_seq[:idx2] + anc2 + human_seq[idx2 + 1 :]
      seq12 = seq1[:idx2] + anc2 + seq1[idx2 + 1 :]
      for m in focales:
        width = m["width"]
        if distancia >= width:
          continue
        starts = range(
            max(0, idx2 - width + 1), min(idx1, len(human_seq) - width) + 1
        )
        mejor = None
        for start in starts:
          for orient in ("+", "-"):
            s_h, s_1, s_2, s_12 = (
                human_seq[start : start + width],
                seq1[start : start + width],
                seq2[start : start + width],
                seq12[start : start + width],
            )
            if orient == "-":
              s_h, s_1, s_2, s_12 = map(
                  reverse_complement, (s_h, s_1, s_2, s_12)
              )
            sh, s1, s2, s12 = [
                puntuacion_pwm(s, m) for s in (s_h, s_1, s_2, s_12)
            ]
            if any(x is None for x in (sh, s1, s2, s12)):
              continue
            if mejor is None or sh > mejor["humano"]:
              mejor = {
                  "humano": sh,
                  "variante1_ancestral": s1,
                  "variante2_ancestral": s2,
                  "ambas_ancestrales": s12,
                  "orientacion": orient,
                  "start": start,
              }
        if mejor is None:
          continue
        delta1 = mejor["variante1_ancestral"] - mejor["humano"]
        delta2 = mejor["variante2_ancestral"] - mejor["humano"]
        delta_both = mejor["ambas_ancestrales"] - mejor["humano"]
        residual = delta_both - delta1 - delta2
        if abs(residual) < 0.05:
          continue
        filas.append({
            "gen_locus": gene,
            "cromosoma": GENES[gene]["chrom"],
            "coordenada_variante_1": p1 + 1,
            "cambio_1": (
                f"{v1['base_ancestral_inferida']}>{v1['base_humana_maf']}"
            ),
            "coordenada_variante_2": p2 + 1,
            "cambio_2": (
                f"{v2['base_ancestral_inferida']}>{v2['base_humana_maf']}"
            ),
            "distancia_bp": distancia,
            "motif_id": m["id"],
            "TF_motif": m["name"],
            "es_AR": m["is_ar"],
            "es_ESR1": m["is_esr1"],
            "es_ESR2": m["is_esr2"],
            "PWM_humano": round(mejor["humano"], 4),
            "delta_variante1_ancestral": round(delta1, 4),
            "delta_variante2_ancestral": round(delta2, 4),
            "delta_ambas_ancestrales": round(delta_both, 4),
            "residuo_no_aditivo_exploratorio": round(residual, 4),
            "orientacion": mejor["orientacion"],
            "ventana_inicio_hg38_1based": lo + mejor["start"] + 1,
            "nota": (
                "Proxy PWM; no demuestra epistasis biológica ni actividad"
                " reguladora"
            ),
        })
      width_ere = 15
      starts_ere = range(
          max(0, idx2 - width_ere + 1), min(idx1, len(human_seq) - width_ere) + 1
      )
      best_ere = None
      for start_ere in starts_ere:
        for orient_ere in ("+", "-"):
          seq_h = human_seq[start_ere : start_ere + width_ere]
          seq_1 = seq1[start_ere : start_ere + width_ere]
          seq_2 = seq2[start_ere : start_ere + width_ere]
          seq_12 = seq12[start_ere : start_ere + width_ere]
          if orient_ere == "-":
            seq_h, seq_1, seq_2, seq_12 = map(
                reverse_complement, (seq_h, seq_1, seq_2, seq_12)
            )
          score_h, score_1, score_2, score_12 = [
              consenso_score(x, "ERE") for x in (seq_h, seq_1, seq_2, seq_12)
          ]
          if any(x is None for x in (score_h, score_1, score_2, score_12)):
            continue
          if best_ere is None or score_h > best_ere["humano"]:
            best_ere = {
                "humano": score_h,
                "variante1_ancestral": score_1,
                "variante2_ancestral": score_2,
                "ambas_ancestrales": score_12,
                "orientacion": orient_ere,
                "start": start_ere,
            }
      if best_ere is not None:
        d1 = best_ere["variante1_ancestral"] - best_ere["humano"]
        d2 = best_ere["variante2_ancestral"] - best_ere["humano"]
        d12 = best_ere["ambas_ancestrales"] - best_ere["humano"]
        nonadditive = d12 - d1 - d2
        if abs(nonadditive) >= 0.05:
          filas.append({
              "gen_locus": gene,
              "cromosoma": GENES[gene]["chrom"],
              "coordenada_variante_1": p1 + 1,
              "cambio_1": (
                  f"{v1['base_ancestral_inferida']}>{v1['base_humana_maf']}"
              ),
              "coordenada_variante_2": p2 + 1,
              "cambio_2": (
                  f"{v2['base_ancestral_inferida']}>{v2['base_humana_maf']}"
              ),
              "distancia_bp": distancia,
              "motif_id": "ERE_CONSENSUS_SHARED_ER_ALPHA_BETA",
              "TF_motif": (
                  "Shared ERalpha/ERbeta estrogen response element consensus"
              ),
              "es_AR": False,
              "es_ESR1": False,
              "es_ESR2": False,
              "ERE_compartido_ERalpha_ERbeta": True,
              "PWM_humano": round(best_ere["humano"], 4),
              "delta_variante1_ancestral": round(d1, 4),
              "delta_variante2_ancestral": round(d2, 4),
              "delta_ambas_ancestrales": round(d12, 4),
              "residuo_no_aditivo_exploratorio": round(nonadditive, 4),
              "orientacion": best_ere["orientacion"],
              "ventana_inicio_hg38_1based": lo + best_ere["start"] + 1,
              "nota": (
                  "Interacción no aditiva de secuencia en consenso ERE"
                  " compartido; es un proxy, no prueba epistasis biológica ni"
                  " activación de ERalpha/ERbeta."
              ),
          })
  return filas


def inferir_secuencia_ancestral_delecion(delecion):
  secuencias = [
      str(delecion.get(f"secuencia_{esp}", "")).upper()
      for esp in ("chimpance", "gorila", "orangutan", "bonobo")
  ]
  secuencias = [s for s in secuencias if s]
  if not secuencias or len({len(s) for s in secuencias}) != 1:
    return "", 0.0
  consenso, soportes = [], []
  for columna in zip(*secuencias):
    votos = Counter(base for base in columna if base in BASES)
    if not votos:
      return "", 0.0
    base, cantidad = votos.most_common(1)[0]
    if cantidad < 2:
      return "", 0.0
    consenso.append(base)
    soportes.append(cantidad / len(secuencias))
  return "".join(consenso), min(soportes) if soportes else 0.0


def mejor_puntaje_evento(
    secuencia,
    inicio_evento,
    fin_evento,
    motivo=None,
    tipo_consenso=None,
    frontera=False,
):
  width = motivo["width"] if motivo else 15
  if frontera:
    lo = max(0, inicio_evento - width + 1)
    hi = min(inicio_evento - 1, len(secuencia) - width)
  else:
    lo = max(0, inicio_evento - width + 1)
    hi = min(fin_evento - 1, len(secuencia) - width)
  mejor = None
  for start in range(lo, hi + 1):
    ventana = secuencia[start : start + width]
    for orientacion in ("+", "-"):
      orientada = (
          reverse_complement(ventana) if orientacion == "-" else ventana
      )
      score = (
          puntuacion_pwm(orientada, motivo)
          if motivo
          else consenso_score(orientada, tipo_consenso)
      )
      if score is not None and (mejor is None or score > mejor["score"]):
        mejor = {
            "score": score,
            "start": start,
            "orientacion": orientacion,
            "secuencia": orientada,
        }
  return mejor


def evaluar_motivos_delecion(
    gene,
    delecion,
    genoma,
    inicio_region,
    motivos,
    enhancers,
    enhancer_idx,
    tracks,
    track_idx,
):
  consenso, confianza = inferir_secuencia_ancestral_delecion(delecion)
  posicion = int(delecion["punto_union_interbase_hg38_0based"])
  if not consenso:
    return {
        "secuencia_ancestral_reconstruida": "",
        "confianza_consenso_minima": confianza,
        "estado_motivos_delecion": "secuencia_ancestral_no_reconstruible",
    }, []
  centro = posicion - inicio_region
  if centro < 0 or centro > len(genoma):
    return {
        "secuencia_ancestral_reconstruida": consenso,
        "confianza_consenso_minima": confianza,
        "estado_motivos_delecion": "fuera_de_region_hg38",
    }, []
  a, b = max(0, centro - FLANCO_MOTIVOS), min(len(genoma), centro + FLANCO_MOTIVOS)
  human_seq = genoma[a:b]
  centro_local = centro - a
  ancestral_seq = (
      human_seq[:centro_local] + consenso + human_seq[centro_local:]
  )
  hits_g = solapamientos(posicion, enhancers, enhancer_idx)
  hits_t = {
      k: solapamientos(posicion, value["intervalos"], track_idx[k])
      for k, value in tracks.items()
  }
  primary = any(
      tracks[k]["primary_mammary"] and hits_t.get(k) for k in tracks
  )
  adipose = any(
      hits_t.get(k) for k in ("adipocyte_general", "adipose_tissue_general")
  )
  regulatory = bool(hits_g or any(hits_t.values()))
  scan = (
      motivos
      if regulatory
      else [m for m in motivos if m["is_ar"] or m["is_esr1"] or m["is_esr2"]]
  )
  effects = []
  for motivo in scan:
    sa = mejor_puntaje_evento(
        ancestral_seq, centro_local, centro_local + len(consenso), motivo=motivo
    )
    sh = mejor_puntaje_evento(
        human_seq, centro_local, centro_local, motivo=motivo, frontera=True
    )
    if sa is None or sh is None:
      continue
    delta = sa["score"] - sh["score"]
    if sa["score"] >= UMBRAL_MOTIVO and delta >= DELTA_MOTIVO:
      predicho = "posible_debilitamiento_motivo_por_delecion"
    elif sh["score"] >= UMBRAL_MOTIVO and delta <= -DELTA_MOTIVO:
      predicho = "posible_ganancia_motivo_en_junction"
    else:
      continue
    effects.append({
        "gen_locus": gene,
        "tipo_evento": "delecion_candidata",
        "cromosoma": GENES[gene]["chrom"],
        "coordenada_evento_hg38_interbase_0based": posicion,
        "longitud_delecion_candidata": len(consenso),
        "motif_id": motivo["id"],
        "TF_motif": motivo["name"],
        "es_motivo_AR": motivo["is_ar"],
        "es_motivo_ESR1": motivo["is_esr1"],
        "es_motivo_ESR2_especifico": motivo["is_esr2"],
        "motivo_familia_ERE": motivo["is_ere_family"],
        "grupo_funcional": clasificar_factor(motivo["name"]),
        "score_ancestral": round(sa["score"], 4),
        "score_humano": round(sh["score"], 4),
        "delta_ancestral_menos_humano": round(delta, 4),
        "efecto_predicho": predicho,
        "ventana_inicio_ancestral_hg38_aproximada_1based": (
            inicio_region + a + sa["start"] + 1
        ),
        "orientacion_ancestral": sa["orientacion"],
        "orientacion_humana": sh["orientacion"],
        "secuencia_ancestral": sa["secuencia"],
        "secuencia_humana": sh["secuencia"],
        "confianza_consenso_minima": round(confianza, 3),
        "solapa_enhancer_global": bool(hits_g),
        "solapa_mamario_primario": primary,
        "solapa_adiposo_general": adipose,
        "tracks_solapados": ";".join(
            k for k, val in hits_t.items() if val
        ),
        "nota": (
            "Comparación exploratoria entre secuencia insertada inferida en"
            " simios y secuencia humana que une los flancos; no confirma"
            " actividad funcional."
        ),
    })
  are_a = mejor_puntaje_evento(
      ancestral_seq,
      centro_local,
      centro_local + len(consenso),
      tipo_consenso="ARE",
  )
  are_h = mejor_puntaje_evento(
      human_seq, centro_local, centro_local, tipo_consenso="ARE", frontera=True
  )
  ere_a = mejor_puntaje_evento(
      ancestral_seq,
      centro_local,
      centro_local + len(consenso),
      tipo_consenso="ERE",
  )
  ere_h = mejor_puntaje_evento(
      human_seq, centro_local, centro_local, tipo_consenso="ERE", frontera=True
  )
  if ere_a is not None and ere_h is not None:
    ere_delta = ere_a["score"] - ere_h["score"]
    if ere_a["score"] >= UMBRAL_CONSENSO_ERE and ere_delta >= DELTA_CONSENSO_ERE:
      ere_effect = "posible_debilitamiento_motivo_por_delecion"
    elif (
        ere_h["score"] >= UMBRAL_CONSENSO_ERE and ere_delta <= -DELTA_CONSENSO_ERE
    ):
      ere_effect = "posible_ganancia_motivo_en_junction"
    else:
      ere_effect = ""
    if ere_effect:
      effects.append({
          "gen_locus": gene,
          "tipo_evento": "delecion_candidata",
          "cromosoma": GENES[gene]["chrom"],
          "coordenada_evento_hg38_interbase_0based": posicion,
          "longitud_delecion_candidata": len(consenso),
          "motif_id": "ERE_CONSENSUS_SHARED_ER_ALPHA_BETA",
          "TF_motif": (
              "Shared ERalpha/ERbeta estrogen response element consensus"
          ),
          "es_motivo_AR": False,
          "es_motivo_ESR1": False,
          "es_motivo_ESR2_especifico": False,
          "motivo_familia_ERE": True,
          "ERE_compartido_ERalpha_ERbeta": True,
          "grupo_funcional": "ERE compartido por ERalpha y ERbeta",
          "score_ancestral": round(ere_a["score"], 4),
          "score_humano": round(ere_h["score"], 4),
          "delta_ancestral_menos_humano": round(ere_delta, 4),
          "efecto_predicho": ere_effect,
          "ventana_inicio_ancestral_hg38_aproximada_1based": (
              inicio_region + a + ere_a["start"] + 1
          ),
          "orientacion_ancestral": ere_a["orientacion"],
          "orientacion_humana": ere_h["orientacion"],
          "secuencia_ancestral": ere_a["secuencia"],
          "secuencia_humana": ere_h["secuencia"],
          "confianza_consenso_minima": round(confianza, 3),
          "solapa_enhancer_global": bool(hits_g),
          "solapa_mamario_primario": primary,
          "solapa_adiposo_general": adipose,
          "tracks_solapados": ";".join(
              k for k, val in hits_t.items() if val
          ),
          "nota": (
              "Comparación de consenso ERE alrededor de una unión de deleción"
              " candidata; no demuestra unión in vivo ni activación"
              " transcripcional."
          ),
      })
  losses = [
      e
      for e in effects
      if e["efecto_predicho"] == "posible_debilitamiento_motivo_por_delecion"
  ]
  gains = [
      e
      for e in effects
      if e["efecto_predicho"] == "posible_ganancia_motivo_en_junction"
  ]
  resumen = {
      "secuencia_ancestral_reconstruida": consenso,
      "confianza_consenso_minima": round(confianza, 3),
      "estado_motivos_delecion": "comparacion_PWM_exploratoria",
      "solapa_enhancer_global_en_junction": bool(hits_g),
      "solapa_anotacion_mamaria_en_junction": primary,
      "solapa_adiposo_general_en_junction": adipose,
      "tracks_solapados_en_junction": ";".join(
          k for k, val in hits_t.items() if val
      ),
      "ARE_ancestral_max": round(are_a["score"], 4) if are_a else "",
      "ARE_humano_junction_max": round(are_h["score"], 4) if are_h else "",
      "ARE_delta_ancestral_menos_humano": round(
          are_a["score"] - are_h["score"], 4
      )
      if are_a and are_h
      else "",
      "ERE_ancestral_max": round(ere_a["score"], 4) if ere_a else "",
      "ERE_humano_junction_max": round(ere_h["score"], 4) if ere_h else "",
      "ERE_delta_ancestral_menos_humano": round(
          ere_a["score"] - ere_h["score"], 4
      )
      if ere_a and ere_h
      else "",
      "n_perdidas_motivos": len(losses),
      "n_ganancias_motivos": len(gains),
      "perdidas_motivos_top": "; ".join(
          f"{e['TF_motif']} {e['delta_ancestral_menos_humano']:+.3f}"
          for e in sorted(
              losses, key=lambda x: x["delta_ancestral_menos_humano"], reverse=True
          )[:5]
      ),
      "ganancias_motivos_top": "; ".join(
          f"{e['TF_motif']} {e['delta_ancestral_menos_humano']:+.3f}"
          for e in sorted(
              gains, key=lambda x: x["delta_ancestral_menos_humano"]
          )[:5]
      ),
  }
  return resumen, effects


def anotar_un_locus(
    gene, variantes, deleciones, segmentos_e, region, motivos, roles_modelo
):
  meta = GENES[gene]
  chrom, inicio, fin = meta["chrom"], region[0], region[1]
  enhancers, tracks, estado = descargar_anotaciones(gene, chrom, inicio, fin)
  enh_idx = construir_indice(enhancers)
  track_idx = {k: construir_indice(v["intervalos"]) for k, v in tracks.items()}
  genoma = descargar_secuencia_hg38(gene, chrom, inicio, fin)
  resultados, efectos, prioridades = [], [], Counter()
  primary = {k for k, v in tracks.items() if v["primary_mammary"]}
  posiciones_regulatorias = {}
  for v in variantes:
    p = int(v["coordenada_hg38_1based"]) - 1
    hits_g = solapamientos(p, enhancers, enh_idx)
    hits_t = {
        k: solapamientos(p, d["intervalos"], track_idx[k])
        for k, d in tracks.items()
    }
    overlaps_cells = any(hits_t.values())
    hit_mammary = any(hits_t.get(k) for k in primary)
    hit_adipose = any(
        hits_t.get(k) for k in ("adipocyte_general", "adipose_tissue_general")
    )
    regulatory = bool(hits_g or overlaps_cells)
    posiciones_regulatorias[p] = regulatory
    row = dict(v)
    row.update({
        "base_hg38": genoma[p - inicio] if inicio <= p < fin else "N",
        "solapa_enhancer_global": bool(hits_g),
        "enhancers_globales": ";".join(x["id"] for x in hits_g),
        "overlap_cualquier_anotacion_celular": overlaps_cells,
        "overlap_enhancer_mamario_primario": hit_mammary,
        "overlap_adiposo_general_no_mamario": hit_adipose,
        "tracks_mamarios_solapados": ";".join(
            k for k in primary if hits_t.get(k)
        ),
        "tracks_celulares_solapados": ";".join(
            k for k, hit in hits_t.items() if hit
        ),
        "estado_modelado": "no_evaluado",
        "motif_scan_applied": False,
        "interpretacion_limite": (
            "Un solapamiento no prueba el gen diana ni la dirección funcional"
            " del elemento."
        ),
    })
    for k, data in tracks.items():
      hit = hits_t[k]
      row[f"overlap_{k}"] = bool(hit)
      row[f"ids_{k}"] = ";".join(x["id"] for x in hit)
      row[f"source_{k}"] = data["source"] if hit else ""
    seq_human = genoma[p - inicio] if 0 <= p - inicio < len(genoma) else "N"
    if seq_human != v["base_humana_maf"]:
      row["validacion_base"] = "discordancia_MAF_hg38"
      resumen, efectos_v = {"estado_modelado": "discordancia_MAF_hg38"}, []
    else:
      row["validacion_base"] = "coincide"
      scan_motifs = (
          motivos
          if regulatory
          else [m for m in motivos if m["is_ar"] or m["is_esr1"] or m["is_esr2"]]
      )
      resumen, efectos_v = puntuar_variante(v, genoma, inicio, scan_motifs)
      row.update(resumen)
      row["motif_scan_applied"] = True
      row["motif_scan_scope"] = (
          "todos_los_motivos_JASPAR_en_anotacion"
          if regulatory
          else "solo_AR_ESR1_ESR2_fuera_de_anotacion"
      )
    for e in efectos_v:
      e2 = dict(e)
      e2.update({
          "solapa_enhancer_global": bool(hits_g),
          "solapa_mamario_primario": hit_mammary,
          "solapa_adiposo_general": hit_adipose,
          "tracks_solapados": ";".join(k for k, hit in hits_t.items() if hit),
      })
      efectos.append(e2)
    ar_drop = False
    try:
      ar_drop = (
          float(row.get("AR_PWM_ancestral", "")) >= UMBRAL_MOTIVO
          and float(row.get("AR_PWM_delta_ancestral_menos_humano", ""))
          >= DELTA_MOTIVO
      )
    except (ValueError, TypeError):
      pass
    are_drop = False
    try:
      are_drop = (
          float(row.get("ARE_consenso_ancestral", "")) >= UMBRAL_CONSENSO_ARE
          and float(row.get("ARE_delta_ancestral_menos_humano", ""))
          >= DELTA_CONSENSO_ARE
      )
    except (ValueError, TypeError):
      pass
    ere_drop = False
    ere_gain = False
    try:
      ere_anc = float(row.get("ERE_consenso_ancestral", ""))
      ere_hum = float(row.get("ERE_consenso_humano", ""))
      ere_delta = float(row.get("ERE_delta_ancestral_menos_humano", ""))
      ere_drop = (
          ere_anc >= UMBRAL_CONSENSO_ERE and ere_delta >= DELTA_CONSENSO_ERE
      )
      ere_gain = (
          ere_hum >= UMBRAL_CONSENSO_ERE and ere_delta <= -DELTA_CONSENSO_ERE
      )
    except (ValueError, TypeError):
      pass
    esr1_gain = False
    try:
      esr1_gain = (
          float(row.get("ESR1_PWM_humano", "")) >= UMBRAL_MOTIVO
          and float(row.get("ESR1_PWM_delta_ancestral_menos_humano", ""))
          <= -DELTA_MOTIVO
      )
    except (ValueError, TypeError):
      pass
    n_losses_other = sum(
        x["efecto_predicho"] == "posible_debilitamiento_motivo"
        and not x["es_motivo_AR"]
        and x.get("motif_id")
        not in {"ARE_CONSENSUS_AR", "ERE_CONSENSUS_SHARED_ER_ALPHA_BETA"}
        for x in efectos_v
    )
    if hits_g and hit_mammary and (ar_drop or are_drop or ere_gain):
      priority = "P1_global_y_mamario_con_posible_perdida_AR_ARE_o_ganancia_ERE"
    elif hit_mammary and (ar_drop or are_drop or ere_gain):
      priority = "P2_mamario_con_posible_perdida_AR_ARE_o_ganancia_ERE"
    elif hits_g and hit_mammary:
      priority = "P3_global_y_mamario_sin_direccion_molecular"
    elif hits_g and (ar_drop or are_drop or ere_gain):
      priority = "P4_global_con_posible_perdida_AR_ARE_o_ganancia_ERE"
    elif hit_mammary and n_losses_other:
      priority = "P5_mamario_con_perdida_otro_motivo"
    elif hits_g and n_losses_other:
      priority = "P6_global_con_perdida_otro_motivo"
    elif regulatory:
      priority = "P7_solapamiento_regulador_sin_efecto_priorizado"
    elif hit_adipose:
      priority = "P8_solo_proxy_adiposo_general"
    elif ar_drop or are_drop or ere_drop or ere_gain:
      priority = "P9_cambio_AR_ARE_ERE_sin_anotacion_disponible"
    else:
      priority = "sin_solapamiento_o_cambio_priorizado"
    row["prioridad"] = priority
    ar_delta_value = float(
        row.get("AR_PWM_delta_ancestral_menos_humano") or 0
    )
    are_delta_value = float(
        row.get("ARE_delta_ancestral_menos_humano") or 0
    )
    row["candidato_perdida_AR_PWM"] = bool(ar_drop)
    row["candidato_perdida_ARE_consenso"] = bool(are_drop)
    row["candidato_perdida_AR_o_ARE"] = bool(ar_drop or are_drop)
    row["proxy_perdida_AR_PWM_sin_filtro_anotacion"] = (
        round(max(0, ar_delta_value), 4) if ar_drop else 0.0
    )
    row["proxy_perdida_ARE_sin_filtro_anotacion"] = (
        round(max(0, are_delta_value), 4) if are_drop else 0.0
    )
    row["proxy_modelo_perdida_frenos_AR"] = (
        round(
            max(
                ar_delta_value if ar_drop else 0.0,
                are_delta_value if are_drop else 0.0,
            ),
            4,
        )
        if (ar_drop or are_drop) and regulatory
        else 0.0
    )
    row["proxy_perdida_frenos_sin_filtro_anotacion"] = (
        round(
            max(
                ar_delta_value if ar_drop else 0.0,
                are_delta_value if are_drop else 0.0,
            ),
            4,
        )
        if (ar_drop or are_drop)
        else 0.0
    )
    row["proxy_perdida_frenos_fuera_anotacion"] = (
        row["proxy_perdida_frenos_sin_filtro_anotacion"]
        if not regulatory
        else 0.0
    )
    row["proxy_modelo_ganancia_aceleradores_ESR1"] = (
        round(
            max(
                0,
                -float(
                    row.get("ESR1_PWM_delta_ancestral_menos_humano") or 0
                ),
            ),
            4,
        )
        if esr1_gain and regulatory
        else 0.0
    )
    ere_delta_value = float(
        row.get("ERE_delta_ancestral_menos_humano") or 0
    )
    row["candidato_perdida_ERE_compartido"] = bool(ere_drop)
    row["candidato_ganancia_ERE_compartido"] = bool(ere_gain)
    row["ERE_familia_compartida_ERalpha_ERbeta"] = True
    row["proxy_perdida_ERE_sin_filtro_anotacion"] = (
        round(max(0.0, ere_delta_value), 4) if ere_drop else 0.0
    )
    row["proxy_ganancia_ERE_sin_filtro_anotacion"] = (
        round(max(0.0, -ere_delta_value), 4) if ere_gain else 0.0
    )
    row["proxy_cambio_ERE_sin_direccion"] = abs(ere_delta_value)
    prioridades[priority] += 1
    resultados.append(row)
  deleciones_anotadas, efectos_deleciones = [], []
  for d in deleciones:
    row = dict(d)
    p = int(d["punto_union_interbase_hg38_0based"])
    if enhancers:
      near = min(
          ((distancia_intervalo(p, x), x) for x in enhancers), key=lambda x: x[0]
      )
      row["distancia_enhancer_global_mas_cercano_bp"], row[
          "enhancer_global_mas_cercano"
      ] = near[0], near[1]["id"]
    for k, data in tracks.items():
      if data["intervalos"]:
        near = min(
            ((distancia_intervalo(p, x), x) for x in data["intervalos"]),
            key=lambda x: x[0],
        )
        row[f"distancia_{k}_mas_cercano_bp"], row[
            f"elemento_{k}_mas_cercano"
        ] = near[0], near[1]["id"]
      else:
        row[f"distancia_{k}_mas_cercano_bp"], row[
            f"elemento_{k}_mas_cercano"
        ] = "", ""
    resumen_del, efectos_del = evaluar_motivos_delecion(
        gene, d, genoma, inicio, motivos, enhancers, enh_idx, tracks, track_idx
    )
    row.update(resumen_del)
    efectos_deleciones.extend(efectos_del)
    row["nota"] = (
        "Un hueco MAF es candidato, no confirma una deleción evolutiva; los"
        " puntajes de motivo comparan una secuencia de simios inferida con la"
        " unión de flancos humanos."
    )
    deleciones_anotadas.append(row)
  efectos.extend(efectos_deleciones)
  interactions = escanear_interacciones(
      gene, variantes, genoma, inicio, motivos, posiciones_regulatorias
  )
  guardar_csv(CARPETA / f"variantes_{gene}_anotadas.csv", resultados)
  guardar_csv(CARPETA / f"efectos_motivo_{gene}.csv", efectos)
  guardar_csv(CARPETA / f"efectos_motivo_deleciones_{gene}.csv", efectos_deleciones)
  guardar_csv(CARPETA / f"deleciones_{gene}_anotadas.csv", deleciones_anotadas)
  guardar_csv(CARPETA / f"interacciones_motivos_{gene}.csv", interactions)
  guardar_csv(CARPETA / f"segmentos_e_MAF_{gene}.csv", segmentos_e)
  resumen = resumir_modelos(
      gene, resultados, efectos, interactions, fin - inicio, roles_modelo
  )
  print(
      f"{gene}: variantes anotadas {len(resultados):,}; efectos de motivo"
      f" incluidos SNP+deleciones {len(efectos):,}; interacciones no aditivas"
      f" candidatas {len(interactions):,}"
  )
  return resumen, resultados, efectos, deleciones_anotadas, interactions, estado, (inicio, fin)


def resumir_modelos(
    gene, variantes, efectos, interacciones, longitud_region, roles_modelo
):
  def numero(fila, campo):
    try:
      valor = fila.get(campo, "")
      return float(valor) if valor not in (None, "") else None
    except (TypeError, ValueError):
      return None

  def anotada(fila):
    return bool(
        fila.get("solapa_enhancer_global")
        or fila.get("overlap_cualquier_anotacion_celular")
        or fila.get("solapa_mamario_primario")
        or fila.get("tracks_celulares_solapados")
    )

  def perdida_ar_pwm(fila):
    anc = numero(fila, "AR_PWM_ancestral")
    delta = numero(fila, "AR_PWM_delta_ancestral_menos_humano")
    return (
        anc is not None
        and delta is not None
        and anc >= UMBRAL_MOTIVO
        and delta >= DELTA_MOTIVO
    )

  def perdida_are(fila):
    anc = numero(fila, "ARE_consenso_ancestral")
    delta = numero(fila, "ARE_delta_ancestral_menos_humano")
    return (
        anc is not None
        and delta is not None
        and anc >= 0.75
        and delta >= (1 / 12)
    )

  def ganancia_esr1(fila):
    hum = numero(fila, "ESR1_PWM_humano")
    delta = numero(fila, "ESR1_PWM_delta_ancestral_menos_humano")
    return (
        hum is not None
        and delta is not None
        and hum >= UMBRAL_MOTIVO
        and delta <= -DELTA_MOTIVO
    )

  def ganancia_esr2(fila):
    hum = numero(fila, "ESR2_PWM_humano")
    delta = numero(fila, "ESR2_PWM_delta_ancestral_menos_humano")
    return (
        hum is not None
        and delta is not None
        and hum >= UMBRAL_MOTIVO
        and delta <= -DELTA_MOTIVO
    )

  regulatory = [v for v in variantes if anotada(v)]
  ar_pwm_all = [v for v in variantes if perdida_ar_pwm(v)]
  are_all = [v for v in variantes if perdida_are(v)]
  ar_pwm_reg = [v for v in ar_pwm_all if anotada(v)]
  are_reg = [v for v in are_all if anotada(v)]
  ar_pwm_unannotated = [v for v in ar_pwm_all if not anotada(v)]
  are_unannotated = [v for v in are_all if not anotada(v)]

  def loss_score(rows, which):
    vals = []
    for v in rows:
      ar_delta = numero(v, "AR_PWM_delta_ancestral_menos_humano") or 0.0
      are_delta = numero(v, "ARE_delta_ancestral_menos_humano") or 0.0
      ar_ok = perdida_ar_pwm(v)
      are_ok = perdida_are(v)
      if which == "AR_PWM":
        value = ar_delta if ar_ok else 0.0
      elif which == "ARE":
        value = are_delta if are_ok else 0.0
      else:
        value = max(
            ar_delta if ar_ok else 0.0, are_delta if are_ok else 0.0
        )
      vals.append(max(0.0, value))
    return round(sum(vals), 4)

  esr1_gains_all = [v for v in variantes if ganancia_esr1(v)]
  esr1_gains_reg = [v for v in esr1_gains_all if anotada(v)]
  esr1_gains_unannotated = [v for v in esr1_gains_all if not anotada(v)]
  esr1_losses = [
      v
      for v in variantes
      if (numero(v, "ESR1_PWM_ancestral") or 0) >= UMBRAL_MOTIVO
      and (numero(v, "ESR1_PWM_delta_ancestral_menos_humano") or 0)
      >= DELTA_MOTIVO
  ]
  esr2_gains_all = [v for v in variantes if ganancia_esr2(v)]
  esr2_losses_all = [
      v
      for v in variantes
      if (numero(v, "ESR2_PWM_ancestral") or 0) >= UMBRAL_MOTIVO
      and (numero(v, "ESR2_PWM_delta_ancestral_menos_humano") or 0)
      >= DELTA_MOTIVO
  ]

  def score_esr1(rows):
    return round(
        sum(
            max(
                0.0,
                -(
                    numero(v, "ESR1_PWM_delta_ancestral_menos_humano") or 0.0
                ),
            )
            for v in rows
            if ganancia_esr1(v)
        ),
        4,
    )

  esr2_uncertain = round(
      sum(
          abs(
              numero(v, "ESR2_PWM_delta_ancestral_menos_humano") or 0.0
          )
          for v in esr2_gains_all + esr2_losses_all
      ),
      4,
  )

  def perdida_ere(fila):
    anc = numero(fila, "ERE_consenso_ancestral")
    delta = numero(fila, "ERE_delta_ancestral_menos_humano")
    return (
        anc is not None
        and delta is not None
        and anc >= UMBRAL_CONSENSO_ERE
        and delta >= DELTA_CONSENSO_ERE
    )

  def ganancia_ere(fila):
    hum = numero(fila, "ERE_consenso_humano")
    delta = numero(fila, "ERE_delta_ancestral_menos_humano")
    return (
        hum is not None
        and delta is not None
        and hum >= UMBRAL_CONSENSO_ERE
        and delta <= -DELTA_CONSENSO_ERE
    )

  ere_losses_all = [v for v in variantes if perdida_ere(v)]
  ere_gains_all = [v for v in variantes if ganancia_ere(v)]
  ere_losses_reg = [v for v in ere_losses_all if anotada(v)]
  ere_gains_reg = [v for v in ere_gains_all if anotada(v)]
  ere_losses_unannotated = [v for v in ere_losses_all if not anotada(v)]
  ere_gains_unannotated = [v for v in ere_gains_all if not anotada(v)]

  def ere_score(rows, direction):
    score = 0.0
    for v in rows:
      delta = numero(v, "ERE_delta_ancestral_menos_humano") or 0.0
      if direction == "gain" and ganancia_ere(v):
        score += max(0.0, -delta)
      elif direction == "loss" and perdida_ere(v):
        score += max(0.0, delta)
    return round(score, 4)

  ere_gain_reg_score = ere_score(ere_gains_reg, "gain")
  ere_gain_all_score = ere_score(ere_gains_all, "gain")
  ere_gain_unannotated_score = ere_score(ere_gains_unannotated, "gain")
  ere_loss_all_score = ere_score(ere_losses_all, "loss")
  ere_loss_reg_score = ere_score(ere_losses_reg, "loss")
  ere_loss_unannotated_score = ere_score(ere_losses_unannotated, "loss")
  ere_interactions = [
      x
      for x in interacciones
      if x.get("motif_id") == "ERE_CONSENSUS_SHARED_ER_ALPHA_BETA"
  ]
  ere_deletion_effects = [
      e
      for e in efectos
      if e.get("motif_id") == "ERE_CONSENSUS_SHARED_ER_ALPHA_BETA"
      and e.get("tipo_evento") == "delecion_candidata"
  ]
  ere_deletion_losses = [
      e
      for e in ere_deletion_effects
      if e.get("efecto_predicho")
      == "posible_debilitamiento_motivo_por_delecion"
  ]
  ere_deletion_gains = [
      e
      for e in ere_deletion_effects
      if e.get("efecto_predicho") == "posible_ganancia_motivo_en_junction"
  ]
  ere_deletion_loss_score = round(
      sum(
          max(0.0, float(e.get("delta_ancestral_menos_humano") or 0))
          for e in ere_deletion_losses
      ),
      4,
  )
  ere_deletion_gain_score = round(
      sum(
          max(0.0, -float(e.get("delta_ancestral_menos_humano") or 0))
          for e in ere_deletion_gains
      ),
      4,
  )
  cyp19_unsigned = sum(1 for x in regulatory) if gene == "CYP19A1" else 0
  ar_role = roles_modelo.get("AR", {"peso": 1.0, "activo": True, "rol": "freno"})
  esr1_role = roles_modelo.get(
      "ESR1_ERalpha", {"peso": 1.0, "activo": True, "rol": "acelerador"}
  )
  esr2_role = roles_modelo.get(
      "ESR2_ERbeta", {"peso": 0.0, "activo": False, "rol": "contexto_no_direccionado"}
  )
  ere_role = roles_modelo.get(
      "ERE_compartido_ERalpha_ERbeta",
      {"peso": 0.0, "activo": False, "rol": "contexto_no_direccionado"},
  )

  loss_reg_score = loss_score(regulatory, "COMBINADO")
  loss_all_score = loss_score(variantes, "COMBINADO")
  loss_unannotated_score = loss_score(
      [v for v in variantes if not anotada(v)], "COMBINADO"
  )
  loss_ar_pwm_all_score = loss_score(variantes, "AR_PWM")
  loss_are_all_score = loss_score(variantes, "ARE")
  gain_reg_score = score_esr1(esr1_gains_reg)
  gain_all_score = score_esr1(esr1_gains_all)
  gain_unannotated_score = score_esr1(esr1_gains_unannotated)

  loss_weighted = (
      loss_reg_score * ar_role["peso"]
      if ar_role["activo"] and ar_role["rol"] == "freno"
      else 0.0
  )
  gain_weighted = (
      gain_reg_score * esr1_role["peso"]
      if esr1_role["activo"] and esr1_role["rol"] == "acelerador"
      else 0.0
  )
  if esr2_role["activo"] and esr2_role["rol"] == "acelerador":
    gain_weighted += (
        score_esr1([v for v in esr2_gains_all if anotada(v)])
        * esr2_role["peso"]
    )
  if esr2_role["activo"] and esr2_role["rol"] == "freno":
    loss_weighted += (
        loss_score([v for v in esr2_losses_all if anotada(v)], "AR_PWM")
        * esr2_role["peso"]
    )
  gain_weighted_without_ere = gain_weighted
  ere_gain_weighted = (
      ere_gain_reg_score * ere_role["peso"]
      if ere_role["activo"] and ere_role["rol"] == "acelerador"
      else 0.0
  )
  gain_weighted += ere_gain_weighted
  gain_weighted_with_ere_sensitivity = round(
      gain_weighted_without_ere + ere_gain_reg_score, 4
  )

  ar_effect_losses = [
      e
      for e in efectos
      if e.get("es_motivo_AR")
      and e.get("motif_id") != "ARE_CONSENSUS_AR"
      and str(e.get("efecto_predicho", "")).startswith(
          "posible_debilitamiento_motivo"
      )
  ]
  ar_effect_losses_reg = [e for e in ar_effect_losses if anotada(e)]
  esr1_effect_gains = [
      e
      for e in efectos
      if e.get("es_motivo_ESR1")
      and e.get("efecto_predicho") == "posible_ganancia_motivo_humano"
  ]
  esr1_effect_gains_reg = [e for e in esr1_effect_gains if anotada(e)]

  return {
      "gen_locus": gene,
      "cromosoma": GENES[gene]["chrom"],
      "longitud_region_bp": longitud_region,
      "variantes_candidatas": len(variantes),
      "variantes_en_anotaciones_reguladoras": len(regulatory),
      "perdidas_motivo_AR_PWM_total_n": len(ar_pwm_all),
      "perdidas_motivo_AR_PWM_en_anotacion_n": len(ar_pwm_reg),
      "perdidas_motivo_AR_PWM_fuera_anotacion_n": len(ar_pwm_unannotated),
      "caidas_ARE_consenso_total_n": len(are_all),
      "caidas_ARE_consenso_en_anotacion_n": len(are_reg),
      "caidas_ARE_consenso_fuera_anotacion_n": len(are_unannotated),
      "proxy_perdida_AR_PWM_todos_candidatos": loss_ar_pwm_all_score,
      "proxy_perdida_ARE_todos_candidatos": loss_are_all_score,
      "perdidas_ERE_compartido_total_n": len(ere_losses_all),
      "perdidas_ERE_compartido_en_anotacion_n": len(ere_losses_reg),
      "perdidas_ERE_compartido_fuera_anotacion_n": len(ere_losses_unannotated),
      "ganancias_ERE_compartido_total_n": len(ere_gains_all),
      "ganancias_ERE_compartido_en_anotacion_n": len(ere_gains_reg),
      "ganancias_ERE_compartido_fuera_anotacion_n": len(ere_gains_unannotated),
      "proxy_perdida_ERE_sin_filtro_anotacion": ere_loss_all_score,
      "proxy_perdida_ERE_en_anotacion": ere_loss_reg_score,
      "proxy_perdida_ERE_fuera_anotacion": ere_loss_unannotated_score,
      "proxy_ganancia_ERE_sin_filtro_anotacion": ere_gain_all_score,
      "proxy_ganancia_ERE_en_anotacion": ere_gain_reg_score,
      "proxy_ganancia_ERE_fuera_anotacion": ere_gain_unannotated_score,
      "perdidas_ERE_deleciones_candidatas_n": len(ere_deletion_losses),
      "ganancias_ERE_deleciones_candidatas_n": len(ere_deletion_gains),
      "interacciones_ERE_consenso_no_aditivas_candidatas_n": len(
          ere_interactions
      ),
      "proxy_perdida_ERE_deleciones_candidatas": ere_deletion_loss_score,
      "proxy_ganancia_ERE_deleciones_candidatas": ere_deletion_gain_score,
      "ERE_peso_aplicado": ere_role["peso"]
      if ere_role["activo"] and ere_role["rol"] == "acelerador"
      else 0.0,
      "proxy_ganancia_ERE_incluida_en_modelo_ponderado": round(
          ere_gain_weighted, 4
      ),
      "proxy_perdida_frenos_AR_o_ARE_sin_filtro_anotacion": loss_all_score,
      "proxy_perdida_frenos_AR_o_ARE_en_anotacion": loss_reg_score,
      "proxy_perdida_frenos_AR_o_ARE_fuera_anotacion": loss_unannotated_score,
      "perdidas_motivo_AR_PWM_en_tabla_efectos_n": len(ar_effect_losses),
      "perdidas_motivo_AR_PWM_en_tabla_efectos_anotadas_n": len(
          ar_effect_losses_reg
      ),
      "ganancias_motivo_ESR1_total_n": len(esr1_gains_all),
      "ganancias_motivo_ESR1_en_anotacion_n": len(esr1_gains_reg),
      "ganancias_motivo_ESR1_fuera_anotacion_n": len(esr1_gains_unannotated),
      "proxy_ganancia_aceleradores_ESR1_sin_filtro_anotacion": gain_all_score,
      "proxy_ganancia_aceleradores_ESR1_en_anotacion": gain_reg_score,
      "proxy_ganancia_aceleradores_ESR1_fuera_anotacion": gain_unannotated_score,
      "perdidas_motivo_ESR1_en_anotacion_n": sum(
          1 for v in esr1_losses if anotada(v)
      ),
      "ganancias_motivo_ESR2_total_n": len(esr2_gains_all),
      "perdidas_motivo_ESR2_total_n": len(esr2_losses_all),
      "proxy_ESR2_direccion_no_asumida_suma_delta_abs": esr2_uncertain,
      "peso_AR_aplicado": ar_role["peso"]
      if ar_role["activo"] and ar_role["rol"] == "freno"
      else 0.0,
      "peso_ESR1_aplicado": esr1_role["peso"]
      if esr1_role["activo"] and esr1_role["rol"] == "acelerador"
      else 0.0,
      "peso_ESR2_aplicado": esr2_role["peso"]
      if esr2_role["activo"] and esr2_role["rol"] in {"acelerador", "freno"}
      else 0.0,
      "modelo_perdida_frenos_ponderado": round(loss_weighted, 4),
      "modelo_ganancia_aceleradores_ponderado": round(gain_weighted, 4),
      "modelo_ganancia_aceleradores_sin_ERE": round(
          gain_weighted_without_ere, 4
      ),
      "modelo_ganancia_aceleradores_con_ERE_sensibilidad": (
          gain_weighted_with_ere_sensitivity
      ),
      "proxy_modelo_combinado_ponderado": round(
          loss_weighted + gain_weighted, 4
      ),
      "proxy_modelo_combinado_con_ERE_sensibilidad": round(
          loss_weighted + gain_weighted_with_ere_sensitivity, 4
      ),
      "cambios_cis_region_CYP19A1_sin_direccion_n": cyp19_unsigned,
      "variantes_por_Mb": round(
          len(variantes) / (longitud_region / 1_000_000), 3
      )
      if longitud_region
      else "",
      "nota": (
          "Los candidatos de secuencia se informan aunque no haya una"
          " anotación conocida. La puntuación con anotación se presenta"
          " aparte. Ningún PWM/consenso demuestra unión in vivo, dirección"
          " reguladora, expresión ni fenotipo. Un solapamiento tampoco"
          " identifica por sí solo el gen diana."
      ),
  }


def calcular_ratio_e2_t(e2_pg_ml, testosterona_ng_dl):
  if (
      e2_pg_ml is None
      or testosterona_ng_dl is None
      or e2_pg_ml < 0
      or testosterona_ng_dl <= 0
  ):
    return None
  e2_molar = (e2_pg_ml * 1e-9) / 272.38
  t_molar = (testosterona_ng_dl * 1e-8) / 288.42
  return e2_molar / t_molar


def crear_plantillas():
  ruta_hormonas = CARPETA / "entradas_modelo_hormonal.csv"
  if not ruta_hormonas.exists():
    datos_escenarios = {
        "humano_adulto_baseline": {
            "e2": "50.0",
            "t": "400.0",
            "er": "1.0",
            "ar": "1.0",
            "nota": "Baseline adulto estándar.",
        },
        "chimpance_adulto_baseline": {
            "e2": "50.0",
            "t": "400.0",
            "er": "1.0",
            "ar": "1.0",
            "nota": "Baseline adulto simio de referencia.",
        },
        "Tanner_I": {
            "e2": "10.0",
            "t": "20.0",
            "er": "1.0",
            "ar": "1.0",
            "nota": "Prepubertal; niveles basales bajos.",
        },
        "Tanner_II": {
            "e2": "25.0",
            "t": "150.0",
            "er": "1.0",
            "ar": "1.0",
            "nota": "Inicio puberal.",
        },
        "ginecomastia": {
            "e2": "90.0",
            "t": "250.0",
            "er": "1.2",
            "ar": "0.8",
            "nota": "Desequilibrio estrógeno/andrógeno favorable a estrógenos.",
        },
        "humano_embarazo": {
            "e2": "12000.0",
            "t": "300.0",
            "er": "1.5",
            "ar": "1.0",
            "nota": "Estrógenos placentarios masivamente elevados.",
        },
        "humano_lactancia": {
            "e2": "20.0",
            "t": "100.0",
            "er": "1.0",
            "ar": "1.0",
            "nota": "Prolactina alta, estrógenos deprimidos.",
        },
        "humano_post_destete": {
            "e2": "40.0",
            "t": "350.0",
            "er": "1.0",
            "ar": "1.0",
            "nota": "Retorno al baseline tras lactancia.",
        },
        "chimpance_embarazo": {
            "e2": "10000.0",
            "t": "300.0",
            "er": "1.5",
            "ar": "1.0",
            "nota": "Embarazo en simios con alta producción placentaria.",
        },
        "chimpance_lactancia": {
            "e2": "20.0",
            "t": "100.0",
            "er": "1.0",
            "ar": "1.0",
            "nota": "Lactancia en simios.",
        },
        "menopausia_atrofia_mamaria": {
            "e2": "10.0",
            "t": "150.0",
            "er": "0.5",
            "ar": "1.0",
            "nota": "Caída ovárica de estrógenos y atrofia tisular.",
        },
        "CAIS_Morris": {
            "e2": "50.0",
            "t": "500.0",
            "er": "1.0",
            "ar": "0.0",
            "nota": (
                "Insensibilidad completa a andrógenos; receptor AR nulo"
                " (funcionalidad 0.0)."
            ),
        },
    }

    filas_h = []
    for escenario, vals in datos_escenarios.items():
      filas_h.append({
          "escenario": escenario,
          "estradiol_pg_ml": vals["e2"],
          "testosterona_ng_dl": vals["t"],
          "capacidad_funcional_ER_relativa": vals["er"],
          "capacidad_funcional_AR_relativa": vals["ar"],
          "epsilon_funcional_en_unidades_de_actividad": "0.000001",
          "notas": vals["nota"],
      })
    guardar_csv(ruta_hormonas, filas_h)

  ruta_roles = CARPETA / "supuestos_modelos_regulatorios.csv"
  roles_predeterminados = [
      {
          "factor": "AR",
          "rol_hipotetico": "freno",
          "peso": 1.0,
          "activado": "si",
          "supuesto": (
              "La pérdida de un motivo AR en un elemento regulador puede reducir"
              " una acción represora; no se asume para todos los sitios."
          ),
      },
      {
          "factor": "ESR1_ERalpha",
          "rol_hipotetico": "acelerador",
          "peso": 1.0,
          "activado": "si",
          "supuesto": (
              "Una ganancia de motivo no implica activación; depende del"
              " tejido, ligando, cofactores y cromatina."
          ),
      },
      {
          "factor": "ESR2_ERbeta",
          "rol_hipotetico": "contexto_no_direccionado",
          "peso": 0.0,
          "activado": "no",
          "supuesto": (
              "No asignar por defecto acción pro o anti-ERalpha."
          ),
      },
      {
          "factor": "ERE_compartido_ERalpha_ERbeta",
          "rol_hipotetico": "contexto_no_direccionado",
          "peso": 0.0,
          "activado": "no",
          "supuesto": (
              "ERE es reconocido por ERalpha y ERbeta; una ganancia predice un"
              " posible sitio de unión, no necesariamente activación. Puede"
              " activarse como acelerador candidato en un análisis de"
              " sensibilidad."
          ),
      },
      {
          "factor": "CYP19A1_aromatasa",
          "rol_hipotetico": "cambio_cis_sin_direccion",
          "peso": 0.0,
          "activado": "no",
          "supuesto": (
              "La aromatasa es una enzima, no tiene PWM propio; se examinan"
              " elementos cis y motivos de TF."
          ),
      },
  ]
  if not ruta_roles.exists():
    guardar_csv(ruta_roles, roles_predeterminados)
  else:
    with open(ruta_roles, newline="", encoding="utf-8-sig") as archivo_roles:
      existentes = list(csv.DictReader(archivo_roles))
    factores_existentes = {
        str(f.get("factor", "")).strip() for f in existentes
    }
    nuevos = [
        f for f in roles_predeterminados if f["factor"] not in factores_existentes
    ]
    if nuevos:
      guardar_csv(ruta_roles, existentes + nuevos)

  ruta_fenotipos = CARPETA / "plantilla_fenotipos_mamarios.csv"
  if not ruta_fenotipos.exists():
    fenotipos = [
        "Tanner_I",
        "Tanner_II",
        "ginecomastia",
        "senos_humanos_persistentes_no_embarazo_lactancia",
        "humano_embarazo",
        "humano_lactancia",
        "humano_post_destete",
        "chimpance_embarazo",
        "chimpance_lactancia",
        "menopausia_atrofia_mamaria",
    ]
    campos = [
        "condicion",
        "especie",
        "edad",
        "medida_morfologica",
        "unidad",
        "estradiol",
        "testosterona",
        "progesterona",
        "fase_ciclo_o_condicion",
        "metodo_medicion",
        "fuente_bibliografica",
        "notas",
    ]
    guardar_csv(
        ruta_fenotipos,
        [{k: (x if k == "condicion" else "") for k in campos} for x in fenotipos],
        campos,
    )


def cargar_roles_modelo():
  ruta = CARPETA / "supuestos_modelos_regulatorios.csv"
  configuracion = {}
  with open(ruta, newline="", encoding="utf-8-sig") as archivo:
    for fila in csv.DictReader(archivo):
      try:
        peso = float(fila.get("peso", 0) or 0)
      except ValueError:
        peso = 0.0
      activo = (
          str(fila.get("activado", "no")).strip().lower()
          in {"si", "sí", "yes", "true", "1"}
      )
      configuracion[fila.get("factor", "")] = {
          "peso": peso,
          "activo": activo,
          "rol": str(fila.get("rol_hipotetico", "")).strip().lower(),
      }
  return configuracion


def procesar_entradas_hormonales():
  ruta = CARPETA / "entradas_modelo_hormonal.csv"
  if not ruta.exists():
    return
  salida = []
  with open(ruta, newline="", encoding="utf-8-sig") as archivo:
    for fila in csv.DictReader(archivo):
      try:
        e2 = float(fila["estradiol_pg_ml"])
        t = float(fila["testosterona_ng_dl"])
        er = float(fila["capacidad_funcional_ER_relativa"])
        ar = float(fila["capacidad_funcional_AR_relativa"])
        epsilon = float(
            fila["epsilon_funcional_en_unidades_de_actividad"] or 1e-6
        )
        ratio = calcular_ratio_e2_t(e2, t)
        e2_molar = e2 * 1e-9 / 272.38
        t_molar = t * 1e-8 / 288.42
        e_act, a_act = e2_molar * er, t_molar * ar
        ratio_func = e_act / (a_act + epsilon)
        salida.append({
            "escenario": fila["escenario"],
            "E2_T_molar": ratio if ratio is not None else "",
            "E2_molar": e2_molar,
            "T_molar": t_molar,
            "actividad_estrogenica_proxy": e_act,
            "actividad_androgenica_proxy": a_act,
            "ratio_funcional_exploratorio": ratio_func,
            "epsilon": epsilon,
            "estado": "modelo_exploratorio_no_validado",
            "nota": (
                "Capacidades relativas no equivalen a afinidad ni actividad"
                " fisiológica medida; si AR=0 el ratio depende de epsilon y no"
                " es un valor clínico."
            ),
        })
      except (ValueError, TypeError, KeyError):
        continue
  guardar_csv(CARPETA / "resultados_modelo_hormonal.csv", salida)
  print(f"Escenarios hormonales calculados con entradas completas: {len(salida)}")


def main():
  crear_plantillas()
  print("ANÁLISIS MULTILOCUS AR / ESR1 / ESR2 / CYP19A1")
  print(
      "Las regiones se obtienen de las coordenadas hg38 presentes en cada MAF;"
      " no se fuerza una ventana común entre cromosomas."
  )
  motivos = cargar_motivos_jaspar()
  roles_modelo = cargar_roles_modelo()
  todos_resumen, todas_variantes, todos_efectos, todas_deleciones, todas_interacciones, todos_estados = (
      [],
      [],
      [],
      [],
      [],
      [],
  )
  disponibilidad_motivos = [
      {
          "matriz": "AR",
          "matrices_detectadas": sum(m["is_ar"] for m in motivos),
          "base_JASPAR": "MA0007",
          "nota": "Motivo de unión de AR; no mide expresión del gen AR.",
      },
      {
          "matriz": "ERalpha_ESR1",
          "matrices_detectadas": sum(m["is_esr1"] for m in motivos),
          "base_JASPAR": "MA0112",
          "nota": "PWM asociado a ERalpha; no representa un sitio exclusivo de ERalpha.",
      },
      {
          "matriz": "ERbeta_ESR2",
          "matrices_detectadas": sum(m["is_esr2"] for m in motivos),
          "base_JASPAR": "detección por nombre/ID en JASPAR",
          "nota": (
              "La falta de PWM ESR2 específico no implica que ERbeta no"
              " reconozca EREs compartidos."
          ),
      },
      {
          "matriz": "ERE_compartido_ERalpha_ERbeta",
          "matrices_detectadas": 1,
          "base_JASPAR": ERE_CONSENSO,
          "nota": (
              "Consenso de 15 bp AGGTCANNN TGACCT; se puntúan 12 posiciones"
              " conservadas y se evalúa como sitio compartido, no exclusivo de"
              " ESR1 o ESR2."
          ),
      },
      {
          "matriz": "CYP19A1_aromatasa",
          "matrices_detectadas": "no_aplica",
          "base_JASPAR": "enzima",
          "nota": (
              "Se analizan cambios cis y motivos de TF en el locus; no existe"
              " PWM propio de la aromatasa."
          ),
      },
  ]
  for gene in GENES:
    try:
      ruta_maf = resolver_maf(gene)
      variantes, deleciones, segmentos_e, region, fuentes = (
          extraer_variantes_y_deleciones(gene, ruta_maf)
      )
      guardar_csv(CARPETA / f"sustituciones_candidatas_{gene}.csv", variantes)
      guardar_csv(
          CARPETA / f"deleciones_candidatas_{gene}_por_bloque.csv", deleciones
      )
      guardar_csv(CARPETA / f"segmentos_e_MAF_{gene}.csv", segmentos_e)
      resumen, vout, eout, dout, iout, estado, real_region = anotar_un_locus(
          gene,
          variantes,
          deleciones,
          segmentos_e,
          region,
          motivos,
          roles_modelo,
      )
      resumen["archivo_MAF"] = ruta_maf.name
      resumen["region_hg38_0based_inicio"] = real_region[0]
      resumen["region_hg38_0based_fin_exclusivo"] = real_region[1]
      resumen["prefijos_MAF_detectados"] = ";".join(
          f"{k}:{v}" for k, v in sorted(fuentes.items())
      )
      resumen["estado_ejecucion"] = "completado"
      todos_resumen.append(resumen)
      todas_variantes.extend(vout)
      todos_efectos.extend(eout)
      todas_deleciones.extend(dout)
      todas_interacciones.extend(iout)
      todos_estados.extend(estado)
    except Exception as error_gene:
      print(f"[ERROR EN {gene}] {error_gene}")
      todos_resumen.append({
          "gen_locus": gene,
          "cromosoma": GENES[gene]["chrom"],
          "estado_ejecucion": "error",
          "detalle_error": str(error_gene),
          "nota": (
              "Este locus no se pudo procesar; revisa archivo MAF, ensamblaje"
              " o disponibilidad de datos."
          ),
      })
  if todos_resumen:
    numeric_fields = [
        "longitud_region_bp",
        "variantes_candidatas",
        "variantes_en_anotaciones_reguladoras",
        "perdidas_motivo_AR_PWM_total_n",
        "perdidas_motivo_AR_PWM_en_anotacion_n",
        "perdidas_motivo_AR_PWM_fuera_anotacion_n",
        "caidas_ARE_consenso_total_n",
        "caidas_ARE_consenso_en_anotacion_n",
        "caidas_ARE_consenso_fuera_anotacion_n",
        "perdidas_ERE_compartido_total_n",
        "perdidas_ERE_compartido_en_anotacion_n",
        "perdidas_ERE_compartido_fuera_anotacion_n",
        "ganancias_ERE_compartido_total_n",
        "ganancias_ERE_compartido_en_anotacion_n",
        "ganancias_ERE_compartido_fuera_anotacion_n",
        "perdidas_ERE_deleciones_candidatas_n",
        "ganancias_ERE_deleciones_candidatas_n",
        "interacciones_ERE_consenso_no_aditivas_candidatas_n",
        "ganancias_motivo_ESR1_total_n",
        "ganancias_motivo_ESR1_en_anotacion_n",
        "ganancias_motivo_ESR1_fuera_anotacion_n",
        "perdidas_motivo_ESR1_en_anotacion_n",
        "ganancias_motivo_ESR2_total_n",
        "perdidas_motivo_ESR2_total_n",
        "interacciones_PWM_no_aditivas_candidatas_n",
        "cambios_cis_region_CYP19A1_sin_direccion_n",
    ]
    total = {
        "gen_locus": "TOTAL_LOCI_PROCESADOS",
        "cromosoma": "multiple",
        "nota": (
            "Suma de loci; las tasas y los puntajes no son fenotipos ni"
            " actividad hormonal medida."
        ),
    }
    for field in numeric_fields:
      total[field] = sum(int(x.get(field, 0) or 0) for x in todos_resumen)
    for field in (
        "proxy_perdida_AR_PWM_todos_candidatos",
        "proxy_perdida_ARE_todos_candidatos",
        "proxy_perdida_ERE_sin_filtro_anotacion",
        "proxy_perdida_ERE_en_anotacion",
        "proxy_perdida_ERE_fuera_anotacion",
        "proxy_ganancia_ERE_sin_filtro_anotacion",
        "proxy_ganancia_ERE_en_anotacion",
        "proxy_ganancia_ERE_fuera_anotacion",
        "proxy_perdida_ERE_deleciones_candidatas",
        "proxy_ganancia_ERE_deleciones_candidatas",
        "proxy_ganancia_ERE_incluida_en_modelo_ponderado",
        "proxy_perdida_frenos_AR_o_ARE_sin_filtro_anotacion",
        "proxy_perdida_frenos_AR_o_ARE_en_anotacion",
        "proxy_perdida_frenos_AR_o_ARE_fuera_anotacion",
        "proxy_ganancia_aceleradores_ESR1_sin_filtro_anotacion",
        "proxy_ganancia_aceleradores_ESR1_en_anotacion",
        "proxy_ganancia_aceleradores_ESR1_fuera_anotacion",
        "proxy_ESR2_direccion_no_asumida_suma_delta_abs",
        "modelo_perdida_frenos_ponderado",
        "modelo_ganancia_aceleradores_ponderado",
        "modelo_ganancia_aceleradores_sin_ERE",
        "modelo_ganancia_aceleradores_con_ERE_sensibilidad",
        "proxy_modelo_combinado_con_ERE_sensibilidad",
    ):
      total[field] = round(
          sum(float(x.get(field, 0) or 0) for x in todos_resumen), 4
      )
    total["proxy_modelo_combinado_ponderado"] = round(
        total["modelo_perdida_frenos_ponderado"]
        + total["modelo_ganancia_aceleradores_ponderado"],
        4,
    )
    total["proxy_modelo_combinado_con_ERE_sensibilidad"] = round(
        total["modelo_perdida_frenos_ponderado"]
        + total["modelo_ganancia_aceleradores_con_ERE_sensibilidad"],
        4,
    )
    total["variantes_por_Mb"] = "no_comparable_suma_de_regiones_distintas"
    todos_resumen.append(total)
  guardar_csv(
      CARPETA / "comparacion_modelos_3_escenarios.csv", todos_resumen
  )
  guardar_csv(CARPETA / "variantes_todos_loci.csv", todas_variantes)
  guardar_csv(CARPETA / "efectos_motivo_todos_loci.csv", todos_efectos)
  guardar_csv(CARPETA / "deleciones_todos_loci.csv", todas_deleciones)
  guardar_csv(
      CARPETA / "interacciones_no_aditivas_todos_loci.csv", todas_interacciones
  )
  guardar_csv(CARPETA / "estado_anotaciones_todos_loci.csv", todos_estados)
  guardar_csv(
      CARPETA / "disponibilidad_matrices_jaspar.csv", disponibilidad_motivos
  )
  crear_plantillas()
  procesar_entradas_hormonales()
  print("\n=== COMPARACIÓN DE MODELOS EXPLORATORIOS ===")
  for row in todos_resumen:
    print(
        f"{row.get('gen_locus')}: "
        f"pérdidas AR PWM={row.get('perdidas_motivo_AR_PWM_total_n', 'NA')} "
        f"(anotadas={row.get('perdidas_motivo_AR_PWM_en_anotacion_n', 'NA')}, "
        f"sin"
        f" anotación={row.get('perdidas_motivo_AR_PWM_fuera_anotacion_n', 'NA')});"
        f" caídas ARE={row.get('caidas_ARE_consenso_total_n', 'NA')} "
        f"(anotadas={row.get('caidas_ARE_consenso_en_anotacion_n', 'NA')}, sin"
        f" anotación={row.get('caidas_ARE_consenso_fuera_anotacion_n', 'NA')});"
        f" proxy AR/ARE total="
        f"{row.get('proxy_perdida_frenos_AR_o_ARE_sin_filtro_anotacion', 'NA')};"
        f" proxy AR/ARE anotado="
        f"{row.get('proxy_perdida_frenos_AR_o_ARE_en_anotacion', 'NA')}; "
        f"ganancias ERE compartido="
        f"{row.get('ganancias_ERE_compartido_total_n', 'NA')} "
        f"(anotadas={row.get('ganancias_ERE_compartido_en_anotacion_n', 'NA')});"
        f" pérdidas ERE={row.get('perdidas_ERE_compartido_total_n', 'NA')};"
        f" proxy ERE gain anotado="
        f"{row.get('proxy_ganancia_ERE_en_anotacion', 'NA')}; modelo con ERE"
        f" sensibilidad="
        f"{row.get('proxy_modelo_combinado_con_ERE_sensibilidad', 'NA')};"
        f" ganancias ESR1={row.get('ganancias_motivo_ESR1_total_n', 'NA')} "
        f"(anotadas={row.get('ganancias_motivo_ESR1_en_anotacion_n', 'NA')});"
        f" modelo pérdida de frenos ponderado="
        f"{row.get('modelo_perdida_frenos_ponderado', 'NA')}; ganancia"
        f" aceleradores ponderada="
        f"{row.get('modelo_ganancia_aceleradores_ponderado', 'NA')};"
        f" combinado={row.get('proxy_modelo_combinado_ponderado', 'NA')}"
    )
  print("\nArchivos principales:")
  for nombre in (
      "comparacion_modelos_3_escenarios.csv",
      "variantes_todos_loci.csv",
      "efectos_motivo_todos_loci.csv",
      "deleciones_todos_loci.csv",
      "interacciones_no_aditivas_todos_loci.csv",
      "estado_anotaciones_todos_loci.csv",
      "disponibilidad_matrices_jaspar.csv",
      "entradas_modelo_hormonal.csv",
      "supuestos_modelos_regulatorios.csv",
      "plantilla_fenotipos_mamarios.csv",
  ):
    print(nombre)


if __name__ == "__main__":
  try:
    main()
  except Exception as error:
    raise SystemExit(
        f"\nERROR: {error}\nVerifica nombres MAF, conexión y que los MAF sean"
        " GRCh38 con hg38.chrX/chr6/chr14/chr15. Los tracks que no se descarguen"
        " se registran como no disponibles; eso no implica ausencia de"
        " actividad."
    ) from error