
# prueba_gseapy_archs4.py
"""
Módulo 1 (Windows): Validación transcriptómica y enriquecimiento funcional
utilizando GSEApy (Enrichr API) y ARCHS4 REST API (/sigpy/).

Genera:
- gseapy_enrichr.csv
- archs4_correlaciones_breast.csv
- archs4_correlaciones_adipose.csv
- archs4_knn_signature.json
- resumen_transcriptomico.json

NOTA:
GEPIA2 no se consulta en esta versión. La importación opcional de
'gepia' no implica que se esté realizando una consulta a ese servicio.
"""

import json
from pathlib import Path
from typing import Dict, List, Any

import pandas as pd
import requests

from config_modelos import GENE_NETWORK_SIGNATURE, ARCHS4_API_URLS


try:
    import gseapy as gp
    HAS_GSEAPY = True
except ImportError:
    gp = None
    HAS_GSEAPY = False


try:
    import gepia
    HAS_GEPIA = True
except ImportError:
    gepia = None
    HAS_GEPIA = False


class TranscriptomicCloudValidator:
    def __init__(self, timeout: int = 45):
        self.timeout = timeout

        self.genes_full = (
            GENE_NETWORK_SIGNATURE["UP_IN_BROKEN_BRAKE"]
            + GENE_NETWORK_SIGNATURE["DOWN_IN_BROKEN_BRAKE"]
            + GENE_NETWORK_SIGNATURE["PIONEER_COFACTORS"]
        )

        self.genes_full = sorted(
            set(
                str(gene).strip().upper()
                for gene in self.genes_full
                if gene and str(gene).strip()
            )
        )

        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json",
            "User-Agent": "TranscriptomicCloudValidator/1.0"
        })

    def run_gseapy_enrichr(self) -> pd.DataFrame:
        """Ejecuta enriquecimiento funcional mediante GSEApy/Enrichr."""
        if not HAS_GSEAPY:
            print("  [!] GSEApy no está instalado en este Python.")
            print("      Instálalo con: py -m pip install gseapy")
            return pd.DataFrame()

        if not self.genes_full:
            print("  [!] La lista de genes está vacía.")
            return pd.DataFrame()

        libraries = [
            "MSigDB_Hallmark_2020",
            "KEGG_2021_Human",
            "GO_Biological_Process_2021",
        ]

        try:
            enrichment = gp.enrichr(
                gene_list=self.genes_full,
                gene_sets=libraries,
                organism="human",
                outdir=None,
                no_plot=True,
                verbose=False,
            )

            df = enrichment.res2d

            if df is None or df.empty:
                print("  [!] Enrichr no devolvió resultados.")
                return pd.DataFrame()

            columns = [
                "Gene_set",
                "Term",
                "Overlap",
                "P-value",
                "Adjusted P-value",
                "Combined Score",
                "Genes",
            ]

            available_columns = [
                column for column in columns if column in df.columns
            ]

            if "Adjusted P-value" in df.columns:
                df = df.sort_values(
                    "Adjusted P-value",
                    ascending=True,
                )

            return df[available_columns].head(15).copy()

        except Exception as exc:
            print(
                "  [!] Error consultando Enrichr mediante GSEApy: "
                f"{type(exc).__name__}: {exc}"
            )
            return pd.DataFrame()

    def _parse_archs4_correlation_json(
        self,
        data: Any,
        target_gene: str,
        meta_filter: str,
    ) -> List[Dict[str, Any]]:
        """
        Intenta interpretar formatos JSON conocidos.

        Importante: solo interpreta los campos reconocidos. Si ARCHS4
        devuelve otra estructura, no inventa correlaciones.
        """
        rows = []

        if not isinstance(data, dict):
            return rows

        positive = data.get(
            "positive_correlated_genes",
            data.get("positive", []),
        )
        negative = data.get(
            "negative_correlated_genes",
            data.get("negative", []),
        )

        def add_gene_row(item, direction):
            if not isinstance(item, dict):
                return

            gene = item.get("gene", item.get("symbol"))
            correlation = item.get("correlation", item.get("r"))

            if gene is None or correlation is None:
                return

            try:
                r_value = float(correlation)
            except (TypeError, ValueError):
                return

            if not -1.0 <= r_value <= 1.0:
                return

            expression = item.get("mean_log_expression")

            try:
                expression = (
                    float(expression) if expression is not None else None
                )
            except (TypeError, ValueError):
                expression = None

            rows.append({
                "Query_Gene": target_gene.upper(),
                "Tissue_Meta_Filter": meta_filter,
                "Direction": direction,
                "Correlated_Gene": str(gene),
                "Pearson_r": round(r_value, 4),
                "Mean_Log_Expression": (
                    round(expression, 3)
                    if expression is not None
                    else None
                ),
            })

        if isinstance(positive, list):
            for item in positive[:10]:
                add_gene_row(item, "Positive_Correlation")

        if isinstance(negative, list):
            for item in negative[:10]:
                add_gene_row(item, "Negative_Correlation")

        # Formato alternativo: diccionario plano con coeficientes.
        # Se acepta solo si el JSON no contenía listas reconocibles.
        if not rows:
            numeric_items = {}

            for key, value in data.items():
                if str(key).upper() == target_gene.upper():
                    continue

                if isinstance(value, (int, float)) and not isinstance(
                    value, bool
                ):
                    try:
                        r_value = float(value)
                    except (TypeError, ValueError):
                        continue

                    if -1.0 <= r_value <= 1.0:
                        numeric_items[str(key)] = r_value

            positive_items = sorted(
                (
                    (gene, value)
                    for gene, value in numeric_items.items()
                    if value > 0
                ),
                key=lambda item: item[1],
                reverse=True,
            )[:10]

            negative_items = sorted(
                (
                    (gene, value)
                    for gene, value in numeric_items.items()
                    if value < 0
                ),
                key=lambda item: item[1],
            )[:10]

            for gene, value in positive_items:
                rows.append({
                    "Query_Gene": target_gene.upper(),
                    "Tissue_Meta_Filter": meta_filter,
                    "Direction": "Positive_Correlation",
                    "Correlated_Gene": gene,
                    "Pearson_r": round(value, 4),
                    "Mean_Log_Expression": None,
                })

            for gene, value in negative_items:
                rows.append({
                    "Query_Gene": target_gene.upper(),
                    "Tissue_Meta_Filter": meta_filter,
                    "Direction": "Negative_Correlation",
                    "Correlated_Gene": gene,
                    "Pearson_r": round(value, 4),
                    "Mean_Log_Expression": None,
                })

        return rows

    def query_archs4_gene_correlation(
        self,
        target_gene: str = "CYP19A1",
        meta_filter: str = "breast",
        k: int = 25,
    ) -> pd.DataFrame:
        """
        Consulta el endpoint de correlación configurado para ARCHS4.

        Primero prueba el filtro tisular y luego una consulta sin filtro.
        La segunda consulta no se considera específica del tejido indicado.
        """
        url = ARCHS4_API_URLS["correlation"]

        payloads_to_try = [
            {
                "gene": target_gene.upper(),
                "species": "human",
                "meta": meta_filter,
                "k": k,
            },
            {
                "gene": target_gene.upper(),
                "species": "human",
                "k": k,
            },
        ]

        for index, payload in enumerate(payloads_to_try):
            try:
                response = self.session.post(
                    url,
                    json=payload,
                    timeout=self.timeout,
                )

                print(
                    f"  ARCHS4 correlation [{meta_filter}], "
                    f"intento {index + 1}: HTTP {response.status_code}"
                )

                if response.status_code != 200:
                    body = response.text[:300].replace("\n", " ")
                    print(f"  Respuesta: {body}")
                    continue

                try:
                    data = response.json()
                except ValueError:
                    print("  [!] La respuesta no contiene JSON válido.")
                    continue

                rows = self._parse_archs4_correlation_json(
                    data,
                    target_gene,
                    meta_filter,
                )

                if rows:
                    df = pd.DataFrame(rows)

                    # El resultado del segundo intento es global; no se
                    # presenta como si fuera específico del tejido.
                    if index == 1:
                        df["Tissue_Meta_Filter"] = (
                            "unfiltered_fallback_requested_" + meta_filter
                        )

                    return df

                print(
                    "  [!] HTTP 200, pero no se reconocieron "
                    "correlaciones en el formato recibido."
                )

            except requests.RequestException as exc:
                print(
                    "  [!] Error de conexión con ARCHS4: "
                    f"{type(exc).__name__}: {exc}"
                )
            except Exception as exc:
                print(
                    "  [!] Error procesando la respuesta de ARCHS4: "
                    f"{type(exc).__name__}: {exc}"
                )

        # No se convierten resultados de enriquecimiento en coeficientes
        # Pearson: un Combined Score no es una correlación génica.
        print(
            "  [!] No se obtuvieron correlaciones interpretables "
            "para el gen y filtro solicitados."
        )
        return pd.DataFrame()

    def query_archs4_knn_signature(self) -> Dict[str, Any]:
        """Busca firmas similares mediante el endpoint k-NN de ARCHS4."""
        url = ARCHS4_API_URLS["knn_signature"]

        payload = {
            "signatures": [
                {
                    "up_genes": GENE_NETWORK_SIGNATURE[
                        "UP_IN_BROKEN_BRAKE"
                    ],
                    "down_genes": GENE_NETWORK_SIGNATURE[
                        "DOWN_IN_BROKEN_BRAKE"
                    ],
                }
            ],
            "species": "human",
            "k": 50,
            "signame": "AR_Loss_Broken_Brake_Mammary_Signature",
        }

        try:
            response = self.session.post(
                url,
                json=payload,
                timeout=self.timeout,
            )

            print(f"  ARCHS4 k-NN: HTTP {response.status_code}")

            if response.status_code != 200:
                print(f"  Respuesta: {response.text[:300]}")
                return {}

            data = response.json()

            if not isinstance(data, dict):
                print("  [!] ARCHS4 k-NN devolvió un formato inesperado.")
                return {}

            samples = data.get("samples", [])
            distances = data.get("distances", [])

            return {
                "signame": data.get(
                    "signame",
                    "AR_Loss_Broken_Brake_Mammary_Signature",
                ),
                "series_count": data.get("series_count", 0),
                "top_matching_gsm_samples": (
                    samples[:15] if isinstance(samples, list) else []
                ),
                "top_distances": (
                    distances[:15] if isinstance(distances, list) else []
                ),
            }

        except requests.RequestException as exc:
            print(
                "  [!] Error de conexión con ARCHS4 k-NN: "
                f"{type(exc).__name__}: {exc}"
            )
        except (ValueError, TypeError) as exc:
            print(
                "  [!] Error interpretando la respuesta k-NN: "
                f"{type(exc).__name__}: {exc}"
            )
        except Exception as exc:
            print(
                "  [!] Error inesperado en ARCHS4 k-NN: "
                f"{type(exc).__name__}: {exc}"
            )

        return {}


def main():
    """Ejecuta las consultas y guarda los resultados en disco."""
    output_dir = (
        Path(__file__).resolve().parent
        / "resultados_transcriptomicos"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    validator = TranscriptomicCloudValidator()

    summary = {
        "genes_input_count": len(validator.genes_full),
        "gseapy_installed": HAS_GSEAPY,
        "gepia_python_module_installed": HAS_GEPIA,
        "tasks": {},
    }

    print("=" * 65)
    print("VALIDACIÓN TRANSCRIPTÓMICA: GSEApy + ARCHS4")
    print("=" * 65)
    print(f"Genes de entrada: {len(validator.genes_full)}")
    print(f"GSEApy instalado: {HAS_GSEAPY}")
    print(f"Módulo Python gepia instalado: {HAS_GEPIA}")
    print(f"Directorio de salida: {output_dir}")
    print()

    # 1. Enriquecimiento funcional.
    print("[1/4] Enriquecimiento funcional mediante Enrichr...")

    try:
        df_enrichr = validator.run_gseapy_enrichr()
        path = output_dir / "gseapy_enrichr.csv"
        df_enrichr.to_csv(
            path,
            index=False,
            encoding="utf-8-sig",
        )

        summary["tasks"]["gseapy_enrichr"] = {
            "rows": len(df_enrichr),
            "file": path.name,
            "status": "results" if not df_enrichr.empty else "empty",
        }

        print(f"  Filas: {len(df_enrichr)}")
        print(f"  Archivo: {path}")

    except Exception as exc:
        summary["tasks"]["gseapy_enrichr"] = {
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }
        print(f"  ERROR: {type(exc).__name__}: {exc}")

    # 2. Correlaciones solicitadas para mama.
    print("\n[2/4] Correlaciones ARCHS4 para mama...")

    try:
        df_breast = validator.query_archs4_gene_correlation(
            target_gene="CYP19A1",
            meta_filter="breast",
            k=25,
        )

        path = output_dir / "archs4_correlaciones_breast.csv"
        df_breast.to_csv(
            path,
            index=False,
            encoding="utf-8-sig",
        )

        summary["tasks"]["archs4_breast"] = {
            "rows": len(df_breast),
            "file": path.name,
            "status": "results" if not df_breast.empty else "empty",
        }

        print(f"  Filas: {len(df_breast)}")
        print(f"  Archivo: {path}")

    except Exception as exc:
        summary["tasks"]["archs4_breast"] = {
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }
        print(f"  ERROR: {type(exc).__name__}: {exc}")

    # 3. Correlaciones solicitadas para tejido adiposo.
    print("\n[3/4] Correlaciones ARCHS4 para tejido adiposo...")

    try:
        df_adipose = validator.query_archs4_gene_correlation(
            target_gene="CYP19A1",
            meta_filter="adipose",
            k=25,
        )

        path = output_dir / "archs4_correlaciones_adipose.csv"
        df_adipose.to_csv(
            path,
            index=False,
            encoding="utf-8-sig",
        )

        summary["tasks"]["archs4_adipose"] = {
            "rows": len(df_adipose),
            "file": path.name,
            "status": "results" if not df_adipose.empty else "empty",
        }

        print(f"  Filas: {len(df_adipose)}")
        print(f"  Archivo: {path}")

    except Exception as exc:
        summary["tasks"]["archs4_adipose"] = {
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }
        print(f"  ERROR: {type(exc).__name__}: {exc}")

    # 4. Búsqueda de firmas similares.
    print("\n[4/4] ARCHS4 k-NN Signature...")

    try:
        knn_result = validator.query_archs4_knn_signature()
        path = output_dir / "archs4_knn_signature.json"

        with path.open("w", encoding="utf-8") as file:
            json.dump(
                knn_result,
                file,
                indent=2,
                ensure_ascii=False,
            )

        sample_count = len(
            knn_result.get("top_matching_gsm_samples", [])
        )

        summary["tasks"]["archs4_knn_signature"] = {
            "matching_samples": sample_count,
            "file": path.name,
            "status": "results" if knn_result else "empty",
        }

        print(f"  Muestras devueltas: {sample_count}")
        print(f"  Archivo: {path}")

    except Exception as exc:
        summary["tasks"]["archs4_knn_signature"] = {
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }
        print(f"  ERROR: {type(exc).__name__}: {exc}")

    # Resumen general.
    summary_path = output_dir / "resumen_transcriptomico.json"

    try:
        with summary_path.open("w", encoding="utf-8") as file:
            json.dump(
                summary,
                file,
                indent=2,
                ensure_ascii=False,
            )

        print(f"\nResumen guardado: {summary_path}")

    except OSError as exc:
        print(f"\n[!] No se pudo guardar el resumen: {exc}")

    print("\n" + "=" * 65)
    print("EJECUCIÓN TERMINADA")
    print("=" * 65)
    print("Comprueba el estado de cada tarea y las filas obtenidas.")
    print("Un CSV vacío no significa que la consulta haya sido exitosa.")


if __name__ == "__main__":
    main()
