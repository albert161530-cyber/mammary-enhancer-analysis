# prueba_encode_bigwig.py
"""
Módulo 2 (Windows): Consulta The ENCODE REST API desde Windows y delega
únicamente la lectura binaria de archivos .bigWig al contenedor Ubuntu.
"""

import os
import subprocess
import requests
import pandas as pd
from typing import List, Dict, Optional
from config_modelos import ENCODE_SEARCH_URL, ENCODE_BIOSAMPLES


class EncodeEpigenomicValidator:
    def __init__(
        self,
        timeout: int = 25,
        container_mode: str = "manual",
        container_name: str = "ubuntu_bio",
        conda_env_Or_python: str = "python3"
    ):
        self.headers = {"accept": "application/json"}
        self.timeout = timeout
        self.container_mode = container_mode
        self.container_name = container_name
        self.conda_python = conda_env_Or_python

    def search_encode_experiments(self, biosample_term: str, assay_title: str = "ATAC-seq", limit: int = 25) -> List[Dict]:
        params = {
            "type": "Experiment",
            "biosample_ontology.term_name": biosample_term,
            "assay_title": assay_title,
            "status": "released",
            "frame": "object",
            "format": "json",
            "limit": str(limit)
        }
        try:
            resp = requests.get(ENCODE_SEARCH_URL, headers=self.headers, params=params, timeout=self.timeout)
            if resp.status_code == 200:
                return resp.json().get("@graph", [])
        except Exception as e:
            print(f"  [!] Error consultando ENCODE en Windows para {biosample_term} ({assay_title}): {e}")
        return []

    def build_encode_inventory(self) -> pd.DataFrame:
        assays = ["ATAC-seq", "DNase-seq", "TF ChIP-seq", "Histone ChIP-seq"]
        rows = []
        for category, biosamples in ENCODE_BIOSAMPLES.items():
            for bio in biosamples:
                for assay in assays:
                    exps = self.search_encode_experiments(bio, assay_title=assay, limit=25)
                    rows.append({
                        "Clinical_Category": category,
                        "Biosample_Term": bio,
                        "Assay_Type": assay,
                        "Released_Experiments_Found": len(exps),
                        "Sample_Accessions": ", ".join([e.get("accession", "") for e in exps[:3]])
                    })
        return pd.DataFrame(rows)

    def get_released_bigwig_url(self, biosample_term: str = "MCF-7") -> Optional[Dict[str, str]]:
        """Busca desde Windows el enlace directo a un archivo bigWig GRCh38 liberado en ENCODE."""
        params_bw = {
            "type": "File",
            "file_format": "bigWig",
            "output_type": "fold change over control",
            "assembly": "GRCh38",
            "biosample_ontology.term_name": biosample_term,
            "status": "released",
            "frame": "object",
            "format": "json",
            "limit": "1"
        }
        try:
            r = requests.get(ENCODE_SEARCH_URL, headers=self.headers, params=params_bw, timeout=self.timeout)
            if r.status_code == 200 and r.json().get("@graph"):
                file_obj = r.json()["@graph"][0]
                return {
                    "biosample": biosample_term,
                    "accession": file_obj.get("accession", "ENCFF_UNKNOWN"),
                    "url": "https://www.encodeproject.org" + file_obj["href"]
                }
        except Exception as e:
            print(f"  [!] Error buscando URL de bigWig en ENCODE ({biosample_term}): {e}")
        return None

    def run_bigwig_via_ubuntu_container(self, bigwig_url: str, output_csv: str, window_bp: int = 150) -> pd.DataFrame:
        if os.path.exists(output_csv):
            print(f"  ✔ Archivo '{output_csv}' detectado en Windows. Cargando señales de pyBigWig...")
            return pd.read_csv(output_csv)

        cmd = None
        if self.container_mode == "wsl":
            cmd = ["wsl", self.conda_python, "ubuntu_bigwig_worker.py", bigwig_url, output_csv, str(window_bp)]
        elif self.container_mode == "docker":
            cmd = ["docker", "exec", self.container_name, self.conda_python, "/workspace/ubuntu_bigwig_worker.py", bigwig_url, f"/workspace/{output_csv}", str(window_bp)]

        if cmd:
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                if res.returncode == 0 and os.path.exists(output_csv):
                    print(res.stdout.strip())
                    return pd.read_csv(output_csv)
            except Exception:
                pass

        print(f"  [Comando para tu contenedor Ubuntu ({output_csv})]:")
        print(f'  python ubuntu_bigwig_worker.py "{bigwig_url}" "{output_csv}" {window_bp}\n')
        return pd.DataFrame()