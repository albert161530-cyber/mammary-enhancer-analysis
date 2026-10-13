import csv
from pathlib import Path

CARPETA = Path(__file__).resolve().parent


def calcular_veredicto():
    archivo_csv = CARPETA / "resultados_alphagenome_final.csv"
    if not archivo_csv.exists():
        print(f"[ERROR] No se encuentra el archivo: {archivo_csv.name}")
        print(
            "Asegúrate de haber corrido 'automatizar_bloques.py' primero para"
            " generarlo."
        )
        return

    print("=== CALCULANDO VEREDICTO EVOLUTIVO (ALPHAGENOME) ===")

    total_ar = 0.0
    count_ar = 0
    total_esr1 = 0.0
    count_esr1 = 0
    errores = 0

    with open(archivo_csv, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if "EXITOSO" not in row.get("estado", ""):
                errores += 1
                continue

            try:
                score = float(row["alphagenome_impact_score"])
            except ValueError:
                continue

            fuerza = row["fuerza_evolutiva"]
            if "AR_loss" in fuerza:
                total_ar += score
                count_ar += 1
            elif "ESR1_gain" in fuerza:
                total_esr1 += score
                count_esr1 += 1

    print(f"\n--- REPORTE CONSOLIDADOS ---")
    print(
        f"• Bloque AR (Pérdida de frenos): {count_ar} variantes evaluadas |"
        f" Impacto Acumulado: {total_ar:.4f}"
    )
    print(
        f"• Bloque ESR1 (Ganancia de aceleradores): {count_esr1} variantes"
        f" evaluadas | Impacto Acumulado: {total_esr1:.4f}"
    )
    if errores > 0:
        print(f"• Variantes con errores omitidas: {errores}")

    print("\n" + "=" * 45)
    if total_ar > total_esr1:
        proporcion = total_ar / total_esr1 if total_esr1 > 0 else float("inf")
        print(
            " [VEREDICTO] La pérdida de enlaces AR domina el transcriptoma"
            f" mamario, siendo {proporcion:.2f} veces más impactante que la"
            " ganancia de ESR1."
        )
    else:
        proporcion = total_esr1 / total_ar if total_ar > 0 else float("inf")
        print(
            " [VEREDICTO] La ganancia de enlaces ESR1 domina el transcriptoma"
            f" mamario, siendo {proporcion:.2f} veces más impactante que la"
            " pérdida de AR."
        )
    print("=" * 45)


if __name__ == "__main__":
    calcular_veredicto()