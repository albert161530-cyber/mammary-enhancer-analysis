
# prueba_encode_bigwig.py
"""
Windows coordina el análisis y Ubuntu/WSL ejecuta pyBigWig.

Genera:
    prueba_bigwig.csv
    prueba_bigwig_gene_summary.csv
    prueba_bigwig_locus_vs_gene.csv

Requisitos en Windows:
    pip install requests pandas

Requisitos en WSL:
    - Distribución Ubuntu.
    - /home/precision_5560/encode-env/bin/python
    - pyBigWig instalado en ese entorno.
    - ubuntu_bigwig_worker.py en la misma carpeta que este script.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests

from config_modelos import ENCODE_SEARCH_URL, ENCODE_BIOSAMPLES


BIGWIG_URL = (
    "https://www.encodeproject.org/files/"
    "ENCFF816WPO/@@download/ENCFF816WPO.bigWig"
)

WINDOW_BP = 150
WORKER_TIMEOUT = 600


class EncodeEpigenomicValidator:
    """Coordina las consultas ENCODE y la ejecución remota de BigWig."""

    def __init__(
        self,
        timeout: int = 25,
        container_mode: str = "wsl",
        container_name: str = "ubuntu_bio",
        conda_env_Or_python: str | None = None,
        worker_path: str | None = None,
        wsl_distribution: str = "Ubuntu",
    ) -> None:
        self.timeout = timeout
        self.container_mode = container_mode.strip().lower()
        self.container_name = container_name
        self.wsl_distribution = wsl_distribution

        self.conda_python = conda_env_Or_python or (
            "/home/precision_5560/encode-env/bin/python"
            if self.container_mode == "wsl"
            else "python3"
        )

        self.worker_path = (
            Path(worker_path).resolve()
            if worker_path
            else Path(__file__).resolve().parent
            / "ubuntu_bigwig_worker.py"
        )

        self.headers = {
            "accept": "application/json",
            "user-agent": "EncodeEpigenomicValidator/1.3",
        }

        self.session = requests.Session()
        self.session.headers.update(self.headers)

    # ---------------------------------------------------------
    # CONSULTAS ENCODE
    # ---------------------------------------------------------

    @staticmethod
    def _encode_search_endpoint() -> str:
        base_url = str(ENCODE_SEARCH_URL).strip().rstrip("/")

        if not base_url:
            raise ValueError(
                "ENCODE_SEARCH_URL está vacío en config_modelos.py."
            )

        if not base_url.endswith("/search"):
            base_url += "/search"

        return base_url + "/"

    def _get_encode_json(self, params: dict) -> dict:
        endpoint = self._encode_search_endpoint()

        response = self.session.get(
            endpoint,
            params=params,
            timeout=self.timeout,
            allow_redirects=True,
        )
        response.raise_for_status()

        payload = response.json()

        if not isinstance(payload, dict):
            raise ValueError("ENCODE no devolvió un objeto JSON.")

        return payload

    def search_encode_experiments(
        self,
        biosample_term: str,
        assay_title: str = "ATAC-seq",
        limit: int = 25,
    ) -> list[dict]:
        params = {
            "type": "Experiment",
            "biosample_ontology.term_name": biosample_term,
            "assay_title": assay_title,
            "status": "released",
            "frame": "object",
            "format": "json",
            "limit": str(limit),
        }

        try:
            payload = self._get_encode_json(params)
            graph = payload.get("@graph", [])
            return graph if isinstance(graph, list) else []

        except (requests.RequestException, ValueError) as exc:
            print(
                f"[ENCODE] Error consultando {biosample_term}, "
                f"{assay_title}: {exc}",
                flush=True,
            )
            return []

    def build_encode_inventory(self) -> pd.DataFrame:
        assays = [
            "ATAC-seq",
            "DNase-seq",
            "TF ChIP-seq",
            "Histone ChIP-seq",
        ]

        rows = []

        for category, biosamples in ENCODE_BIOSAMPLES.items():
            for biosample in biosamples:
                for assay in assays:
                    experiments = self.search_encode_experiments(
                        biosample_term=biosample,
                        assay_title=assay,
                        limit=25,
                    )

                    rows.append({
                        "Clinical_Category": category,
                        "Biosample_Term": biosample,
                        "Assay_Type": assay,
                        "Released_Experiments_Found": len(experiments),
                        "Sample_Accessions": ", ".join(
                            str(exp["accession"])
                            for exp in experiments[:3]
                            if exp.get("accession")
                        ),
                    })

        return pd.DataFrame(rows)

    def get_released_bigwig_url(
        self,
        biosample_term: str = "MCF-7",
        assay_title: str | None = None,
        output_type: str = "fold change over control",
        assembly: str = "GRCh38",
        limit: int = 10,
    ) -> dict | None:
        params = {
            "type": "File",
            "file_format": "bigWig",
            "output_type": output_type,
            "assembly": assembly,
            "biosample_ontology.term_name": biosample_term,
            "status": "released",
            "frame": "object",
            "format": "json",
            "limit": str(limit),
        }

        if assay_title:
            params["assay_title"] = assay_title

        try:
            payload = self._get_encode_json(params)

            for item in payload.get("@graph", []):
                if not isinstance(item, dict):
                    continue

                accession = item.get("accession")
                href = item.get("href")

                if not accession:
                    continue

                if href:
                    url = (
                        href if href.startswith("https://")
                        else "https://www.encodeproject.org" + href
                    )
                else:
                    url = (
                        "https://www.encodeproject.org/files/"
                        f"{accession}/@@download/{accession}.bigWig"
                    )

                parsed = urlparse(url)

                if (
                    parsed.scheme == "https"
                    and parsed.hostname
                    and (
                        parsed.hostname.endswith("encodeproject.org")
                        or parsed.hostname.endswith(".amazonaws.com")
                    )
                ):
                    return {
                        "biosample": biosample_term,
                        "accession": str(accession),
                        "url": url,
                        "assembly": str(
                            item.get("assembly") or assembly
                        ),
                        "output_type": str(
                            item.get("output_type") or output_type
                        ),
                        "assay_title": str(
                            item.get("assay_title") or assay_title or ""
                        ),
                    }

        except (requests.RequestException, ValueError) as exc:
            print(f"[ENCODE] Error buscando BigWig: {exc}", flush=True)

        return None

    # ---------------------------------------------------------
    # UTILIDADES DE EJECUCIÓN
    # ---------------------------------------------------------

    @staticmethod
    def _run_command(
        cmd: list[str],
        timeout: int,
    ) -> subprocess.CompletedProcess:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )

    def _windows_path_to_wsl(self, path: str) -> str:
        result = self._run_command(
            [
                "wsl",
                "-d", self.wsl_distribution,
                "--", "wslpath", "-a",
                os.path.abspath(path),
            ],
            timeout=30,
        )

        if result.returncode != 0:
            raise RuntimeError(
                result.stderr.strip()
                or f"No se pudo convertir la ruta: {path}"
            )

        converted = result.stdout.strip()

        if not converted:
            raise RuntimeError(
                f"WSL devolvió una ruta vacía para {path}"
            )

        return converted

    # ---------------------------------------------------------
    # EJECUCIÓN DEL WORKER EN WSL
    # ---------------------------------------------------------

    def _execute_wsl(
        self,
        bigwig_url: str,
        output_path: str,
        window_bp: int,
        timeout: int,
    ) -> bool:
        if not self.worker_path.is_file():
            raise FileNotFoundError(
                f"No existe el worker: {self.worker_path}"
            )

        worker_wsl = self._windows_path_to_wsl(
            str(self.worker_path)
        )
        output_wsl = self._windows_path_to_wsl(output_path)

        cmd = [
            "wsl",
            "-d", self.wsl_distribution,
            "--",
            self.conda_python,
            worker_wsl,
            bigwig_url,
            output_wsl,
            str(window_bp),
        ]

        print("[Windows] Ejecutando worker en WSL...", flush=True)
        print(f"[Windows] Distribución: {self.wsl_distribution}", flush=True)
        print(f"[Windows] Python WSL: {self.conda_python}", flush=True)

        result = self._run_command(cmd, timeout=timeout)

        if result.stdout.strip():
            print(result.stdout.strip(), flush=True)

        if result.stderr.strip():
            print(result.stderr.strip(), file=sys.stderr, flush=True)

        if result.returncode != 0:
            print(
                f"[Windows] El worker falló. Código: "
                f"{result.returncode}",
                flush=True,
            )
            return False

        return True

    # ---------------------------------------------------------
    # INTERFAZ PRINCIPAL
    # ---------------------------------------------------------

    def run_bigwig_via_ubuntu_container(
        self,
        bigwig_url: str,
        output_csv: str,
        window_bp: int = 150,
        force_refresh: bool = False,
        worker_timeout: int = 600,
    ) -> pd.DataFrame:
        if window_bp < 1:
            raise ValueError("window_bp debe ser positivo.")

        parsed_url = urlparse(bigwig_url)

        if parsed_url.scheme != "https" or not parsed_url.netloc:
            raise ValueError("La URL BigWig debe ser HTTPS.")

        output_path = os.path.abspath(output_csv)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        if not self.worker_path.is_file():
            raise FileNotFoundError(
                f"No se encontró el worker: {self.worker_path}"
            )

        if self.container_mode == "manual":
            print("[Windows] Modo manual: no se ejecuta el worker.")
            return pd.DataFrame()

        if self.container_mode != "wsl":
            raise ValueError(
                "Este script está configurado para WSL. "
                "Usa container_mode='wsl'."
            )

        if os.path.isfile(output_path) and not force_refresh:
            print(
                f"[Windows] Usando CSV existente: {output_path}",
                flush=True,
            )
            return pd.read_csv(output_path)

        succeeded = self._execute_wsl(
            bigwig_url=bigwig_url,
            output_path=output_path,
            window_bp=window_bp,
            timeout=worker_timeout,
        )

        if not succeeded:
            return pd.DataFrame()

        if not os.path.isfile(output_path):
            print(
                "[Windows] ERROR: el worker terminó, pero no existe "
                f"el CSV principal: {output_path}",
                flush=True,
            )
            return pd.DataFrame()

        df = pd.read_csv(output_path)

        print(
            f"[Windows] CSV principal recuperado: {len(df)} filas.",
            flush=True,
        )

        return df

    def close(self) -> None:
        self.session.close()

    def __enter__(self) -> "EncodeEpigenomicValidator":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


# -------------------------------------------------------------
# PUNTO DE ENTRADA: SIN ESTE BLOQUE EL SCRIPT NO HACE NADA
# -------------------------------------------------------------

def main() -> int:
    print("=" * 60, flush=True)
    print(" ENCODE: COMPARACIÓN DE LOCI CONTRA GENES", flush=True)
    print("=" * 60, flush=True)

    base_dir = Path(__file__).resolve().parent
    output_csv = base_dir / "prueba_bigwig.csv"

    expected_files = [
        output_csv,
        output_csv.with_name(
            f"{output_csv.stem}_gene_summary.csv"
        ),
        output_csv.with_name(
            f"{output_csv.stem}_locus_vs_gene.csv"
        ),
    ]

    try:
        with EncodeEpigenomicValidator(
            container_mode="wsl",
            wsl_distribution="Ubuntu",
            conda_env_Or_python=(
                "/home/precision_5560/encode-env/bin/python"
            ),
        ) as validator:
            df = validator.run_bigwig_via_ubuntu_container(
                bigwig_url=BIGWIG_URL,
                output_csv=str(output_csv),
                window_bp=WINDOW_BP,
                force_refresh=True,
                worker_timeout=WORKER_TIMEOUT,
            )

        if df.empty:
            print(
                "[ERROR] No se obtuvieron resultados. "
                "Revisa los mensajes anteriores.",
                flush=True,
            )
            return 1

        print("\n=== SEÑALES DE LOS LOCI ===", flush=True)
        print(df.to_string(index=False), flush=True)

        print("\n=== ARCHIVOS DE SALIDA ===", flush=True)

        missing = []

        for path in expected_files:
            if path.is_file():
                print(f"[OK] {path}", flush=True)
            else:
                print(f"[FALTA] {path}", flush=True)
                missing.append(path.name)

        if missing:
            print(
                "\n[AVISO] El worker no produjo todos los archivos "
                "esperados. Comprueba que ubuntu_bigwig_worker.py "
                "contenga la versión que genera los tres CSV.",
                flush=True,
            )
            return 1

        print("\nAnálisis terminado.", flush=True)
        return 0

    except subprocess.TimeoutExpired:
        print(
            f"[ERROR] WSL superó el límite de {WORKER_TIMEOUT} segundos.",
            flush=True,
        )
        return 1

    except Exception as exc:
        print(
            f"[ERROR] {type(exc).__name__}: {exc}",
            file=sys.stderr,
            flush=True,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
