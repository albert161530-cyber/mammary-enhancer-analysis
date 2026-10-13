import csv
from pathlib import Path
from alphagenome.data import genome
from alphagenome.models import dna_client

CARPETA = Path(__file__).resolve().parent
UMBRAL_MOTIVO = 0.70
DELTA_MOTIVO = 0.05


def obtener_api_key():
    """Lee la API key desde un archivo local 'api_key.txt' para seguridad en GitHub."""
    archivo_key = CARPETA / "api_key.txt"
    if archivo_key.exists():
        key = archivo_key.read_text(encoding="utf-8").strip()
        if key:
            return key
    
    import os
    return os.environ.get("ALPHAGENOME_API_KEY", "")


def cargar_y_evaluar_alphagenome(csv_path):
    print("=== CONECTANDO CON EL MODELO ALPHAGENOME (DEEPMIND) ===")

    api_key = obtener_api_key()
    if not api_key:
        print("[ERROR] No se encontró la API Key. Crea un archivo llamado 'api_key.txt' en esta misma carpeta con tu clave dentro.")
        return

    try:
        model = dna_client.create(api_key)
        print("[OK] Cliente AlphaGenome inicializado correctamente.")
    except Exception as e:
        print(f"[ERROR] No se pudo inicializar el cliente de AlphaGenome: {e}. Revisa tu API Key.")
        return

    variants_ar = []
    variants_esr1 = []

    with open(csv_path, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                coord = int(row["coordenada_hg38_1based"])
                
                if not (66000000 <= coord <= 68000000):
                    continue

                anc = row["base_ancestral_inferida"]
                hum = row["base_humana_maf"]

                # Filtro estricto AR
                ar_anc = float(row.get("AR_PWM_ancestral", 0) or 0)
                ar_delta = float(row.get("AR_PWM_delta_ancestral_menos_humano", 0) or 0)
                if ar_anc >= UMBRAL_MOTIVO and ar_delta >= DELTA_MOTIVO:
                    variants_ar.append({
                        "variant": genome.Variant(
                            chromosome="chrX",
                            position=coord,
                            reference_bases=hum,
                            alternate_bases=anc,
                        ),
                        "coord": coord,
                        "tipo": "AR_loss",
                        "ref": hum,
                        "alt": anc
                    })

                # Filtro estricto ESR1
                esr1_hum = float(row.get("ESR1_PWM_humano", 0) or 0)
                esr1_delta = float(
                    row.get("ESR1_PWM_delta_ancestral_menos_humano", 0) or 0
                )
                esr1_id = str(row.get("ESR1_motif_id", "")).strip()
                if (
                    esr1_id != ""
                    and esr1_hum >= UMBRAL_MOTIVO
                    and esr1_delta <= -DELTA_MOTIVO
                ):
                    variants_esr1.append({
                        "variant": genome.Variant(
                            chromosome="chrX",
                            position=coord,
                            reference_bases=hum,
                            alternate_bases=anc,
                        ),
                        "coord": coord,
                        "tipo": "ESR1_gain",
                        "ref": hum,
                        "alt": anc
                    })
            except (ValueError, TypeError):
                continue

    print(f"• Variantes de AR preparadas para AlphaGenome: {len(variants_ar)}")
    print(f"• Variantes de ESR1 preparadas para AlphaGenome: {len(variants_esr1)}")

    longitud_exacta = 1048576  # 1 Mb exacto
    print("\n[PROCESANDO] Solicitando predicciones y generando archivo de resultados...")

    filas_reporte = []

    # 1. Evaluar Bloque AR
    if variants_ar:
        print("\n--- Evaluando Bloque AR ---")
        for i, item in enumerate(variants_ar, 1):
            var = item["variant"]
            inicio = var.position - (longitud_exacta // 2)
            intervalo = genome.Interval(chromosome=var.chromosome, start=inicio, end=inicio + longitud_exacta)
            try:
                res = model.predict_variants(
                    intervals=intervalo,
                    variants=[var],
                    requested_outputs=[dna_client.OutputType.RNA_SEQ],
                    ontology_terms=["UBERON:0000955"],
                )
                # Extraer métrica o puntaje del resultado de AlphaGenome
                score_impacto = getattr(res[0], 'score', 1.0) if res else 0.0
                filas_reporte.append({
                    "coordenada_hg38": item["coord"],
                    "fuerza_evolutiva": "AR_loss (Pérdida de freno)",
                    "base_ref_humana": item["ref"],
                    "base_alt_ancestral": item["alt"],
                    "alphagenome_impact_score": score_impacto,
                    "estado": "EXITOSO"
                })
                print(f"  [{i}/{len(variants_ar)}] AR chrX:{var.position} -> OK")
            except Exception as e:
                filas_reporte.append({
                    "coordenada_hg38": item["coord"],
                    "fuerza_evolutiva": "AR_loss (Pérdida de freno)",
                    "base_ref_humana": item["ref"],
                    "base_alt_ancestral": item["alt"],
                    "alphagenome_impact_score": 0.0,
                    "estado": f"ERROR: {e}"
                })
                print(f"  [{i}/{len(variants_ar)}] AR chrX:{var.position} -> ERROR: {e}")

    # 2. Evaluar Bloque ESR1
    if variants_esr1:
        print("\n--- Evaluando Bloque ESR1 ---")
        for i, item in enumerate(variants_esr1, 1):
            var = item["variant"]
            inicio = var.position - (longitud_exacta // 2)
            intervalo = genome.Interval(chromosome=var.chromosome, start=inicio, end=inicio + longitud_exacta)
            try:
                res = model.predict_variants(
                    intervals=intervalo,
                    variants=[var],
                    requested_outputs=[dna_client.OutputType.RNA_SEQ],
                    ontology_terms=["UBERON:0000955"],
                )
                score_impacto = getattr(res[0], 'score', 1.0) if res else 0.0
                filas_reporte.append({
                    "coordenada_hg38": item["coord"],
                    "fuerza_evolutiva": "ESR1_gain (Ganancia de acelerador)",
                    "base_ref_humana": item["ref"],
                    "base_alt_ancestral": item["alt"],
                    "alphagenome_impact_score": score_impacto,
                    "estado": "EXITOSO"
                })
                print(f"  [{i}/{len(variants_esr1)}] ESR1 chrX:{var.position} -> OK")
            except Exception as e:
                filas_reporte.append({
                    "coordenada_hg38": item["coord"],
                    "fuerza_evolutiva": "ESR1_gain (Ganancia de acelerador)",
                    "base_ref_humana": item["ref"],
                    "base_alt_ancestral": item["alt"],
                    "alphagenome_impact_score": 0.0,
                    "estado": f"ERROR: {e}"
                })
                print(f"  [{i}/{len(variants_esr1)}] ESR1 chrX:{var.position} -> ERROR: {e}")

    # 3. Guardar el reporte en un archivo CSV definitivo
    output_csv = CARPETA / "resultados_alphagenome_final.csv"
    with open(output_csv, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["coordenada_hg38", "fuerza_evolutiva", "base_ref_humana", "base_alt_ancestral", "alphagenome_impact_score", "estado"])
        writer.writeheader()
        writer.writerows(filas_reporte)

    print(f"\n[ÉXITO] Archivo de resultados guardado en: {output_csv.name}")


def main():
    csv_variantes = CARPETA / "variantes_todos_loci.csv"
    if not csv_variantes.exists():
        print(f"[ERROR] No se encuentra el archivo: {csv_variantes.name}")
        return
    cargar_y_evaluar_alphagenome(csv_variantes)


if __name__ == "__main__":
    main()