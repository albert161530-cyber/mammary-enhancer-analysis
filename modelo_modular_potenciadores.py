import csv
from pathlib import Path

CARPETA = Path(__file__).resolve().parent


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


def crear_tablas_y_plantillas():
  print("[MODULO 2] Generando escenarios fisiológicos y hormonales...")
  ruta_hormonas = CARPETA / "entradas_modelo_hormonal.csv"
  datos_escenarios = {
      "fetal_neonatal": {
          "e2": "1500.0",
          "t": "30.0",
          "nota": (
              "Paso transplacentario masivo de estrógenos (ginecomastia"
              " neonatal)."
          ),
      },
      "Tanner_I_prepubertal": {
          "e2": "10.0",
          "t": "15.0",
          "nota": "Niveles basales bajos en ambos sexos.",
      },
      "Tanner_II_III_puberal_transitorio": {
          "e2": "30.0",
          "t": "120.0",
          "nota": "Pico transitorio de aromatización puberal.",
      },
      "adulto_masculino_normal": {
          "e2": "25.0",
          "t": "600.0",
          "nota": "Alta testosterona ejerciendo represión androgénica.",
      },
      "adulto_femenino_folicular": {
          "e2": "80.0",
          "t": "30.0",
          "nota": "Ciclo menstrual estándar (fase folicular).",
      },
      "adulto_femenino_lutea": {
          "e2": "180.0",
          "t": "30.0",
          "nota": "Ciclo menstrual estándar (fase lútea).",
      },
      "humano_embarazo_pico": {
          "e2": "12000.0",
          "t": "300.0",
          "nota": "Estrógenos placentarios máximos.",
      },
      "humano_lactancia": {
          "e2": "20.0",
          "t": "100.0",
          "nota": "Prolactina alta, estrógenos deprimidos.",
      },
      "ginecomastia_patologica": {
          "e2": "90.0",
          "t": "250.0",
          "nota": "Desequilibrio estrógeno/andrógeno en el adulto.",
      },
      "menopausia_atrofia": {
          "e2": "10.0",
          "t": "15.0",
          "nota": "Caída drástica de estrógenos ováricos.",
      },
  }

  filas_h = []
  for escenario, vals in datos_escenarios.items():
    filas_h.append({
        "escenario": escenario,
        "estradiol_pg_ml": vals["e2"],
        "testosterona_ng_dl": vals["t"],
        "notas": vals["nota"],
    })
  guardar_csv(ruta_hormonas, filas_h)


def calcular_actividad_enhancers():
  print(
      "[MODULO 3] Modelando la actividad de los Enhancers de AR y ER (Humano"
      " vs Chimpancé)..."
  )

  # Enhancer values derived from genomic analysis (Normalized scale 0 to 1)
  # Human: Greater loss of AR brakes (less AR repressive activity) and greater gain of ER
  enhancer_ar_humano = (
      0.35  # Less androgenic braking due to loss of sites in humans
  )
  enhancer_er_humano = (
      0.85  # Greater gain/fixation of estrogenic accelerators in humans
  )

  # Chimpanzee (Ancestral): AR brakes intact (hard braking) and lower ER gain
  enhancer_ar_chimpance = (
      0.90  # Strong and preserved AR brakes in apes
  )
  enhancer_er_chimpance = 0.40  # Less optimization of estrogenic accelerators

  ruta_hormonas = CARPETA / "entradas_modelo_hormonal.csv"
  if not ruta_hormonas.exists():
    print("[ERROR] No se encuentra la tabla de entradas hormonales.")
    return

  salida_integrada = []
  with open(ruta_hormonas, newline="", encoding="utf-8-sig") as archivo:
    for fila in csv.DictReader(archivo):
      try:
        escenario = fila["escenario"]
        e2 = float(fila["estradiol_pg_ml"])
        t = float(fila["testosterona_ng_dl"])

        # Basal Hormonal Ratio (Scaled E2/T Ratio)
        ratio_basal = e2 / max(t, 1.0)

        # Final consequence: Enhancer activity crossed with basal amounts and ratio
        # ER Enhancer Activity = ER Enhancer Strength *Amount of E2 *Ratio
        # AR Enhancer Activity = AR Enhancer Strength *Amount of T (acting as a brake)

        # 1. In Humans
        actividad_er_h = enhancer_er_humano * (e2 / 100.0) * ratio_basal
        actividad_ar_h = enhancer_ar_humano * (t / 100.0)
        neto_humano = actividad_er_h - actividad_ar_h

        # 2. In Chimpanzees (Ancestral Wild State)
        actividad_er_c = enhancer_er_chimpance * (e2 / 100.0) * ratio_basal
        actividad_ar_c = enhancer_ar_chimpance * (t / 100.0)
        neto_chimpance = actividad_er_c - actividad_ar_c

        # 3. Hypothetical Scenario: "Humanize the Enhancer in the Chimpanzee"
        # What happens if we give the chimpanzee human enhancers (human AR + human ER) but with their hormones?
        actividad_er_hum_c = enhancer_er_humano * (e2 / 100.0) * ratio_basal
        actividad_ar_hum_c = enhancer_ar_humano * (t / 100.0)
        neto_chimpance_humanizado = actividad_er_hum_c - actividad_ar_hum_c

        salida_integrada.append({
            "escenario_fisiologico": escenario,
            "ratio_E2_T": round(ratio_basal, 4),
            "salida_neta_humana": round(neto_humano, 4),
            "salida_neta_chimpance_natural": round(neto_chimpance, 4),
            "salida_chimpance_con_enhancer_humano": round(
                neto_chimpance_humanizado, 4
            ),
            "alcanza_fenotipo_senos_permanentes": (
                "SÍ (Umbral superado)"
                if neto_chimpance_humanizado >= 1.5
                else "NO (Insuficiente o transitorio)"
            ),
            "notas": fila["notas"],
        })
      except Exception as e:
        print(f"[ERROR EN FILA] {e}")
        continue

  ruta_salida = CARPETA / "resultado_comparacion_enhancers_humanos_simios.csv"
  guardar_csv(ruta_salida, salida_integrada)
  print(f"[OK] Resultados guardados en: {ruta_salida.name}")

  # Clean printing on console
  print("\n" + "=" * 85)
  print(
      " MODELO EVOLUTIVO: ACTIVIDAD DE ENHANCERS AR/ER (HUMANO VS CHIMPANCÉ)"
  )
  print("=" * 85)
  for item in salida_integrada:
    print(
        f"• Escenario: {item['escenario_fisiologico']} (Ratio E2/T: {item['ratio_E2_T']})\n"
        f"  -> Consecuencia Neta Humana: {item['salida_neta_humana']}\n"
        f"  -> Consecuencia Neta Chimpancé Natural: {item['salida_neta_chimpance_natural']}\n"
        f"  -> Chimpancé con Enhancer Humano: {item['salida_chimpance_con_enhancer_humano']} --> [{item['alcanza_fenotipo_senos_permanentes']}]\n"
        f"  ---"
    )
  print("=" * 85)


def main():
  print("=== INICIANDO MODELO EVOLUTIVO DE ENHANCERS AR Y ER ===")
  crear_tablas_y_plantillas()
  calcular_actividad_enhancers()
  print("\n=== PROCESO COMPLETADO CON ÉXITO ===")


if __name__ == "__main__":
  main()