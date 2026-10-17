# prueba_scenic_celloracle.py

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from config_modelos import DRIVER_VARIANTS


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

LOGGER = logging.getLogger("single_cell_grn")


class SingleCellGRNAdapter:
    def __init__(self, enhancer_half_width: int = 250):
        if enhancer_half_width < 1:
            raise ValueError("enhancer_half_width debe ser mayor que cero.")

        self.half_width = enhancer_half_width

    def build_candidate_regions(self) -> pd.DataFrame:
        rows = []

        for variant in DRIVER_VARIANTS:
            if variant.pos < 1:
                raise ValueError(
                    f"Posición genómica inválida: {variant.variant_id}"
                )

            if not variant.chrom.startswith("chr"):
                raise ValueError(
                    f"Cromosoma inválido: {variant.chrom}"
                )

            start = variant.pos - self.half_width - 1
            end = variant.pos + self.half_width

            if start < 0 or end <= start:
                raise ValueError(
                    f"Intervalo genómico inválido: {variant.variant_id}"
                )

            rows.append({
                "Variant_ID": variant.variant_id,
                "Chromosome": variant.chrom,
                "Position_1based": variant.pos,
                "Start_0based": start,
                "End_0based_exclusive": end,
                "Window_Length_bp": end - start,
                "Configured_Target_Gene": variant.gene,
                "AlphaGenome_L1_Weight_Imported": variant.l1_score,
                "Reference_Allele_From_Config": variant.ref_human,
                "Ancestral_Allele_From_Config": variant.alt_ancestral,
                "Analysis_Status": "candidate_not_experimentally_validated"
            })

        if not rows:
            raise ValueError("DRIVER_VARIANTS no contiene variantes.")

        return pd.DataFrame(rows)

    def export_scenicplus_eregulon_metadata(self) -> pd.DataFrame:
        """
        Exporta metadatos de regiones candidatas con columnas compatibles
        con el orquestador. No representa eRegulons inferidos por SCENIC+.

        Las relaciones TF-gen y región-gen se dejan vacías porque no deben
        inventarse antes de ejecutar y revisar el análisis real de SCENIC+.
        Las coordenadas son BED (0-based, end-exclusive).
        """
        annotations = self.build_candidate_regions()

        result = pd.DataFrame({
            "Region_Name_SCENIC": (
                annotations["Chromosome"].astype(str)
                + ":"
                + annotations["Start_0based"].astype(str)
                + "-"
                + annotations["End_0based_exclusive"].astype(str)
            ),
            "Chromosome": annotations["Chromosome"],
            "Start": annotations["Start_0based"],
            "End": annotations["End_0based_exclusive"],
            "Variant_ID": annotations["Variant_ID"],
            "TF_Regulator": pd.NA,
            "Target_Gene": annotations["Configured_Target_Gene"],
            "TF_to_Gene_Relationship": pd.NA,
            "Region_to_Gene_Relationship": pd.NA,
            "Evolutionary_State_hg38": pd.NA,
            "AlphaGenome_L1_Weight": (
                annotations["AlphaGenome_L1_Weight_Imported"]
            ),
            "Analysis_Status": (
                "candidate_region_not_SCENICplus_eregulon"
            ),
        })

        return result

    def export_celloracle_base_grn_bed(self) -> pd.DataFrame:
        """
        Exporta regiones candidatas con las columnas de identificación
        que espera el orquestador. Esto NO es una base_GRN validada de
        CellOracle: no incluye TFs ni motivos hipotéticos.
        """
        annotations = self.build_candidate_regions()

        result = pd.DataFrame({
            "peak_id": (
                annotations["Chromosome"].astype(str)
                + "_"
                + annotations["Start_0based"].astype(str)
                + "_"
                + annotations["End_0based_exclusive"].astype(str)
            ),
            "gene_short_name": (
                annotations["Configured_Target_Gene"].astype(str)
            ),
            "Variant_ID": annotations["Variant_ID"],
            "Analysis_Status": (
                "candidate_peak_not_validated_celloracle_base_grn"
            ),
        })

        return result

    @staticmethod
    def _bed_name(variant_id: str) -> str:
        return (
            variant_id
            .replace(":", "_")
            .replace(">", "_to_")
            .replace("/", "_")
        )

    def export_candidate_files(
        self,
        output_dir: Path
    ) -> tuple[pd.DataFrame, Path, Path]:

        output_dir.mkdir(parents=True, exist_ok=True)

        annotations = self.build_candidate_regions()

        annotations_path = output_dir / "driver_region_annotations.csv"
        annotations.to_csv(annotations_path, index=False)

        bed_rows = []

        for row in annotations.itertuples(index=False):
            bed_rows.append({
                "chrom": row.Chromosome,
                "start": row.Start_0based,
                "end": row.End_0based_exclusive,
                "name": self._bed_name(row.Variant_ID)
            })

        bed = pd.DataFrame(bed_rows)

        bed_path = output_dir / "driver_variants.bed"
        bed.to_csv(
            bed_path,
            sep="\t",
            index=False,
            header=False
        )

        scenicplus_regions = (
            output_dir
            / "scenicplus_candidate_region_sets"
            / "AR_loss_candidates"
        )
        scenicplus_regions.mkdir(parents=True, exist_ok=True)

        scenicplus_bed = (
            scenicplus_regions / "AR_loss_candidates.bed"
        )

        bed.to_csv(
            scenicplus_bed,
            sep="\t",
            index=False,
            header=False
        )

        LOGGER.info("Variantes exportadas: %s", len(annotations))
        LOGGER.info("Anotaciones: %s", annotations_path)
        LOGGER.info("BED de variantes: %s", bed_path)

        return annotations, bed_path, scenicplus_regions.parent

    @staticmethod
    def _read_matrix(path: Path) -> pd.DataFrame:
        suffix = path.suffix.lower()

        if suffix == ".parquet":
            frame = pd.read_parquet(path)
        elif suffix in {".tsv", ".bed"}:
            frame = pd.read_csv(path, sep="\t")
        else:
            frame = pd.read_csv(path)

        frame = frame.loc[
            :,
            ~frame.columns.astype(str).str.match(r"^Unnamed")
        ]

        return frame

    @staticmethod
    def _parse_peak_id(peak_id: str):
        parts = str(peak_id).rsplit("_", 2)

        if len(parts) != 3:
            return None

        chrom, start, end = parts

        try:
            start = int(start)
            end = int(end)
        except ValueError:
            return None

        if start < 0 or end <= start:
            return None

        if not chrom.startswith("chr"):
            return None

        return chrom, start, end

    def export_candidate_peak_overlaps(
        self,
        annotations: pd.DataFrame,
        base_grn: pd.DataFrame,
        output_dir: Path
    ) -> pd.DataFrame:
        required = {"peak_id", "gene_short_name"}

        missing = required - set(base_grn.columns)

        if missing:
            raise ValueError(
                "La matriz CellOracle no contiene las columnas "
                f"necesarias: {sorted(missing)}"
            )

        peak_coordinates = []

        for index, peak_id in enumerate(
            base_grn["peak_id"].astype(str)
        ):
            parsed = self._parse_peak_id(peak_id)

            if parsed is not None:
                chrom, start, end = parsed
                peak_coordinates.append({
                    "row_index": index,
                    "peak_id": peak_id,
                    "peak_chrom": chrom,
                    "peak_start": start,
                    "peak_end": end
                })

        peaks = pd.DataFrame(peak_coordinates)

        results = []

        for variant in annotations.itertuples(index=False):
            if peaks.empty:
                matches = peaks
            else:
                matches = peaks[
                    (peaks["peak_chrom"] == variant.Chromosome)
                    & (
                        peaks["peak_start"]
                        < variant.End_0based_exclusive
                    )
                    & (
                        peaks["peak_end"]
                        > variant.Start_0based
                    )
                ]

            if matches.empty:
                results.append({
                    "Variant_ID": variant.Variant_ID,
                    "Configured_Target_Gene":
                        variant.Configured_Target_Gene,
                    "Peak_ID": None,
                    "Peak_Target_Gene": None,
                    "Overlap_bp": 0,
                    "TF_Motifs_Present_In_Base_GRN": "",
                    "Peak_Overlap_Status":
                        "no_matching_peak_in_supplied_base_GRN"
                })
                continue

            for peak in matches.itertuples(index=False):
                peak_row = base_grn.iloc[peak.row_index]

                overlap_bp = max(
                    0,
                    min(
                        peak.peak_end,
                        variant.End_0based_exclusive
                    )
                    - max(
                        peak.peak_start,
                        variant.Start_0based
                    )
                )

                motif_columns = [
                    column
                    for column in base_grn.columns
                    if column not in {
                        "peak_id",
                        "gene_short_name"
                    }
                ]

                present_tfs = []

                for column in motif_columns:
                    value = pd.to_numeric(
                        pd.Series([peak_row[column]]),
                        errors="coerce"
                    ).iloc[0]

                    if pd.notna(value) and value > 0:
                        present_tfs.append(str(column))

                peak_gene = str(peak_row["gene_short_name"])

                results.append({
                    "Variant_ID": variant.Variant_ID,
                    "Configured_Target_Gene":
                        variant.Configured_Target_Gene,
                    "Peak_ID": peak.peak_id,
                    "Peak_Target_Gene": peak_gene,
                    "Overlap_bp": overlap_bp,
                    "TF_Motifs_Present_In_Base_GRN":
                        ";".join(present_tfs),
                    "Peak_Overlap_Status":
                        "overlapping_peak_found"
                })

        result = pd.DataFrame(results)

        output_path = (
            output_dir / "celloracle_candidate_peak_overlap.csv"
        )

        result.to_csv(output_path, index=False)

        LOGGER.info(
            "Correspondencias de regiones exportadas: %s",
            output_path
        )

        return result

    def run_celloracle(
        self,
        h5ad_path: Path,
        base_grn_path: Path,
        annotations: pd.DataFrame,
        output_dir: Path,
        cluster_column: str,
        embedding_key: str,
        counts_layer: str | None,
        n_pcs: int,
        knn: int,
        alpha: float,
        pvalue: float,
        max_edges: int,
        n_jobs: int,
        simulate_gene: str | None,
        simulate_value: float,
        grn_unit: str,
        propagation_steps: int
    ) -> dict:

        try:
            import scanpy as sc
            import celloracle as co
        except ImportError as exc:
            raise RuntimeError(
                "Falta CellOracle o Scanpy en este entorno de Python."
            ) from exc

        if not h5ad_path.is_file():
            raise FileNotFoundError(
                f"No existe el archivo de expresión: {h5ad_path}"
            )

        if not base_grn_path.is_file():
            raise FileNotFoundError(
                f"No existe la matriz reguladora: {base_grn_path}"
            )

        LOGGER.info("Cargando datos celulares...")
        adata = sc.read_h5ad(str(h5ad_path))

        if cluster_column not in adata.obs.columns:
            raise ValueError(
                f"No existe la columna '{cluster_column}' en adata.obs. "
                f"Columnas disponibles: {list(adata.obs.columns)}"
            )

        if embedding_key not in adata.obsm:
            raise ValueError(
                f"No existe '{embedding_key}' en adata.obsm. "
                f"Embeddings disponibles: {list(adata.obsm.keys())}"
            )

        if counts_layer:
            if counts_layer not in adata.layers:
                raise ValueError(
                    f"No existe la capa de conteos '{counts_layer}'. "
                    f"Capas disponibles: {list(adata.layers.keys())}"
                )

            adata.X = adata.layers[counts_layer].copy()
        else:
            LOGGER.warning(
                "Se utilizará adata.X como matriz de conteos brutos. "
                "Comprueba que contenga conteos sin normalizar."
            )

        LOGGER.info("Cargando matriz reguladora...")
        base_grn = self._read_matrix(base_grn_path)

        required_columns = {"peak_id", "gene_short_name"}
        missing = required_columns - set(base_grn.columns)

        if missing:
            raise ValueError(
                f"Faltan columnas en la matriz reguladora: {missing}"
            )

        synthetic_columns = {
            "AR_motif_ancestral",
            "AR_motif_human_hg38",
            "AR_motif_CRISPR_rescue",
            "ESR1_motif",
            "FOXA1_pioneer_motif"
        }

        if synthetic_columns.intersection(base_grn.columns):
            raise ValueError(
                "La matriz contiene columnas de motivos hipotéticos "
                "del antiguo script. No constituyen una matriz base "
                "validada de CellOracle. Proporciona una matriz "
                "generada a partir de datos reales de motivos/regiones."
            )

        tf_columns = [
            column
            for column in base_grn.columns
            if column not in required_columns
        ]

        if not tf_columns:
            raise ValueError(
                "La matriz base no contiene columnas de factores "
                "de transcripción."
            )

        if base_grn.empty:
            raise ValueError("La matriz reguladora está vacía.")

        self.export_candidate_peak_overlaps(
            annotations=annotations,
            base_grn=base_grn,
            output_dir=output_dir
        )

        if adata.n_obs < 3 or adata.n_vars < 3:
            raise ValueError(
                "Los datos celulares tienen muy pocas células o genes."
            )

        effective_pcs = min(
            n_pcs,
            50,
            adata.n_obs - 1,
            adata.n_vars - 1
        )

        effective_knn = min(knn, adata.n_obs - 1)

        LOGGER.info(
            "CellOracle: %s células, %s genes y %s regiones reguladoras.",
            adata.n_obs,
            adata.n_vars,
            len(base_grn)
        )

        oracle = co.Oracle()

        oracle.import_anndata_as_raw_count(
            adata=adata,
            cluster_column_name=cluster_column,
            embedding_name=embedding_key
        )

        oracle.import_TF_data(
            TF_info_matrix=base_grn
        )

        LOGGER.info("Calculando PCA...")
        oracle.perform_PCA()

        LOGGER.info("Calculando imputación KNN...")
        oracle.knn_imputation(
            n_pca_dims=effective_pcs,
            k=effective_knn,
            balanced=True
        )

        # CellOracle espera "whole" o el nombre concreto de un grupo celular,
        # no el nombre de la columna que contiene los grupos.
        if grn_unit != "whole":
            available_clusters = set(adata.obs[cluster_column].dropna().astype(str))
            if grn_unit not in available_clusters:
                raise ValueError(
                    f"--grn-unit '{grn_unit}' no coincide con ningún valor de "
                    f"adata.obs['{cluster_column}']. Valores disponibles: "
                    f"{sorted(available_clusters)}. Usa 'whole' para una red global."
                )

        LOGGER.info("Infiriendo la red reguladora para la unidad %s...", grn_unit)
        links = oracle.get_links(
            cluster_name_for_GRN_unit=grn_unit,
            alpha=alpha,
            n_jobs=n_jobs
        )

        links.filter_links(
            p=pvalue,
            weight="coef_abs",
            threshold_number=max_edges
        )

        links_path = (
            output_dir / "celloracle_filtered_links.celloracle.links"
        )

        links.to_hdf5(file_path=str(links_path))

        filtered_links = getattr(links, "filtered_links", None)

        if isinstance(filtered_links, dict):
            for cluster_name, frame in filtered_links.items():
                if not isinstance(frame, pd.DataFrame):
                    continue

                safe_name = re.sub(
                    r"[^A-Za-z0-9_.-]",
                    "_",
                    str(cluster_name)
                )

                frame.to_csv(
                    output_dir
                    / f"celloracle_links_{safe_name}.csv",
                    index=True
                )

        oracle.get_cluster_specific_TFdict_from_Links(
            links_object=links
        )

        oracle.fit_GRN_for_simulation(
            GRN_unit=grn_unit,
            alpha=alpha,
            use_cluster_specific_TFdict=True,
            verbose_level=1
        )

        simulation_performed = False

        if simulate_gene is not None:
            if simulate_value < 0:
                raise ValueError(
                    "El valor de perturbación debe ser no negativo."
                )

            available_genes = set(map(str, oracle.adata.var_names))

            if simulate_gene not in available_genes:
                raise ValueError(
                    f"El gen '{simulate_gene}' no existe en los datos "
                    "de expresión de CellOracle. Comprueba sus "
                    "identificadores en adata.var_names."
                )

            LOGGER.info(
                "Simulando perturbación de %s con valor %s...",
                simulate_gene,
                simulate_value
            )

            oracle.simulate_shift(
                perturb_condition={
                    simulate_gene: simulate_value
                },
                GRN_unit=grn_unit,
                n_propagation=propagation_steps,
                clip_delta_X=True
            )

            if "simulated_count" not in oracle.adata.layers:
                raise RuntimeError(
                    "CellOracle no generó la capa simulated_count."
                )

            baseline = oracle.adata.X
            simulated = oracle.adata.layers["simulated_count"]

            if hasattr(baseline, "toarray"):
                baseline = baseline.toarray()

            if hasattr(simulated, "toarray"):
                simulated = simulated.toarray()

            baseline = np.asarray(baseline)
            simulated = np.asarray(simulated)

            gene_summary = pd.DataFrame({
                "Gene": list(map(str, oracle.adata.var_names)),
                "Mean_Baseline": baseline.mean(axis=0),
                "Mean_Simulated": simulated.mean(axis=0),
                "Mean_Change": (
                    simulated.mean(axis=0)
                    - baseline.mean(axis=0)
                )
            }).sort_values(
                "Mean_Change",
                key=lambda values: values.abs(),
                ascending=False
            )

            gene_summary.to_csv(
                output_dir / "celloracle_simulation_gene_changes.csv",
                index=False
            )

            simulation_performed = True

        oracle_path = (
            output_dir / "celloracle_model.celloracle.oracle"
        )

        oracle.to_hdf5(file_path=str(oracle_path))

        summary = {
            "cells": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "base_grn_regions": int(len(base_grn)),
            "base_grn_tf_columns": int(len(tf_columns)),
            "pca_components": int(effective_pcs),
            "knn_neighbors": int(effective_knn),
            "alpha": alpha,
            "filtered_links_file": str(links_path),
            "oracle_model_file": str(oracle_path),
            "simulation_performed": simulation_performed,
            "simulated_gene": simulate_gene
        }

        with (output_dir / "celloracle_summary.json").open(
            "w",
            encoding="utf-8"
        ) as file:
            json.dump(summary, file, indent=2, ensure_ascii=False)

        LOGGER.info("CellOracle terminó su análisis.")

        return summary

    @staticmethod
    def _is_inside(path: Path, parent: Path) -> bool:
        try:
            path.resolve().relative_to(parent.resolve())
            return True
        except ValueError:
            return False

    def _build_scenicplus_region_overlay(
        self,
        existing_regions: Path,
        candidate_regions: Path,
        overlay_dir: Path,
        output_dir: Path
    ) -> Path:
        """Copia los BED originales y añade candidatos sin sobrescribirlos."""
        if not existing_regions.is_dir():
            raise NotADirectoryError(
                f"No existe el directorio de regiones original: {existing_regions}"
            )
        if not candidate_regions.is_dir():
            raise NotADirectoryError(
                f"No existe el directorio de regiones candidatas: {candidate_regions}"
            )

        source_beds = [
            path for path in existing_regions.rglob("*.bed")
            if path.is_file() and not self._is_inside(path, output_dir)
        ]
        candidate_beds = [
            path for path in candidate_regions.rglob("*.bed") if path.is_file()
        ]
        if not source_beds:
            raise ValueError(
                f"No se encontraron archivos BED originales en: {existing_regions}"
            )
        if not candidate_beds:
            raise ValueError(
                f"No se encontraron archivos BED candidatos en: {candidate_regions}"
            )

        # Se trabaja sobre una copia independiente: el directorio original no se modifica.
        if overlay_dir.exists():
            shutil.rmtree(overlay_dir)
        overlay_dir.mkdir(parents=True, exist_ok=True)
        used_paths = set()

        for source in source_beds:
            relative = source.relative_to(existing_regions)
            destination = overlay_dir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            used_paths.add(relative.as_posix().casefold())

        # Si hay colisión de nombre/ruta, se renombra el candidato, nunca el original.
        for source in candidate_beds:
            relative = source.relative_to(candidate_regions)
            destination = overlay_dir / relative
            if relative.as_posix().casefold() in used_paths or destination.exists():
                destination = destination.with_name(
                    f"{source.stem}__candidate{source.suffix}"
                )
                counter = 2
                while destination.exists():
                    destination = destination.with_name(
                        f"{source.stem}__candidate_{counter}{source.suffix}"
                    )
                    counter += 1
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            used_paths.add(destination.relative_to(overlay_dir).as_posix().casefold())

        LOGGER.info(
            "BED originales copiados y BED candidatos añadidos sin sobrescritura: %s",
            overlay_dir
        )
        LOGGER.warning(
            "Comprueba que el workflow SCENIC+ reconozca los BED candidatos añadidos "
            "como conjuntos de regiones válidos; copiarlos no garantiza que el workflow los use."
        )
        return overlay_dir

    def run_scenicplus(
        self,
        config_path: Path,
        candidate_regions: Path,
        output_dir: Path,
        cores: int
    ) -> dict:

        try:
            import yaml
        except ImportError as exc:
            raise RuntimeError(
                "Falta PyYAML para leer la configuración de SCENIC+."
            ) from exc

        if not config_path.is_file():
            raise FileNotFoundError(
                f"No existe la configuración de SCENIC+: {config_path}"
            )

        pipeline_root = config_path.parent.parent

        snakefile = pipeline_root / "workflow" / "Snakefile"

        if not snakefile.is_file():
            raise FileNotFoundError(
                "No se encontró workflow/Snakefile. "
                "La ruta de configuración debe apuntar a "
                "Snakemake/config/config.yaml de un pipeline inicializado."
            )

        snakemake_executable = shutil.which("snakemake")

        if snakemake_executable is None:
            raise RuntimeError(
                "No se encuentra snakemake en este entorno de Python. "
                "Activa el entorno donde está instalado SCENIC+."
            )

        with config_path.open("r", encoding="utf-8") as file:
            config = yaml.safe_load(file)

        if not isinstance(config, dict):
            raise ValueError("La configuración YAML no es válida.")

        input_data = config.get("input_data", {})

        required_inputs = [
            "cisTopic_obj_fname",
            "GEX_anndata_fname",
            "region_set_folder",
            "ctx_db_fname",
            "dem_db_fname",
            "path_to_motif_annotations"
        ]

        resolved_inputs = {}

        for key in required_inputs:
            value = input_data.get(key)

            if not value:
                raise ValueError(
                    f"Falta el campo input_data.{key} en SCENIC+."
                )

            path = Path(value).expanduser()

            if not path.is_absolute():
                candidates = [
                    pipeline_root / path,
                    config_path.parent / path
                ]

                existing = next(
                    (candidate for candidate in candidates
                     if candidate.exists()),
                    None
                )

                if existing is None:
                    raise FileNotFoundError(
                        f"No existe la entrada de SCENIC+ '{key}': "
                        f"{value}"
                    )

                path = existing

            if not path.exists():
                raise FileNotFoundError(
                    f"No existe la entrada de SCENIC+ '{key}': {path}"
                )

            resolved_inputs[key] = path.resolve()

        original_region_folder = resolved_inputs["region_set_folder"]

        overlay_dir = (
            output_dir / "scenicplus_region_sets_with_candidates"
        )

        overlay_dir = self._build_scenicplus_region_overlay(
            existing_regions=original_region_folder,
            candidate_regions=candidate_regions,
            overlay_dir=overlay_dir,
            output_dir=output_dir
        )

        for key, path in resolved_inputs.items():
            input_data[key] = str(path)

        input_data["region_set_folder"] = str(overlay_dir.resolve())
        config["input_data"] = input_data

        output_data = config.get("output_data", {})

        for key, value in output_data.items():
            if not isinstance(value, str) or not value.strip():
                continue

            path = Path(value).expanduser()

            if not path.is_absolute():
                path = pipeline_root / path

            path = path.resolve()
            path.parent.mkdir(parents=True, exist_ok=True)

            output_data[key] = str(path)

        config["output_data"] = output_data

        general_params = config.get("params_general", {})
        temp_dir = general_params.get("temp_dir")

        if isinstance(temp_dir, str) and temp_dir.strip():
            temp_path = Path(temp_dir).expanduser()

            if not temp_path.is_absolute():
                temp_path = pipeline_root / temp_path

            temp_path = temp_path.resolve()
            temp_path.mkdir(parents=True, exist_ok=True)
            general_params["temp_dir"] = str(temp_path)

        config["params_general"] = general_params

        effective_config = (
            output_dir / "scenicplus_config_with_candidates.yaml"
        )

        with effective_config.open("w", encoding="utf-8") as file:
            yaml.safe_dump(
                config,
                file,
                sort_keys=False,
                allow_unicode=True
            )

        command = [
            snakemake_executable,
            "--snakefile",
            str(snakefile),
            "--configfile",
            str(effective_config),
            "--cores",
            str(cores)
        ]

        LOGGER.info("Iniciando el pipeline real de SCENIC+...")
        LOGGER.info("Configuración: %s", effective_config)

        completed = subprocess.run(
            command,
            cwd=str(pipeline_root),
            check=False
        )

        if completed.returncode != 0:
            raise RuntimeError(
                "SCENIC+ terminó con código de error "
                f"{completed.returncode}. Revisa los registros de "
                "Snakemake y los archivos de entrada."
            )

        summary = {
            "status": "completed",
            "configuration": str(effective_config),
            "region_sets": str(overlay_dir),
            "snakefile": str(snakefile),
            "cores": cores
        }

        with (output_dir / "scenicplus_summary.json").open(
            "w",
            encoding="utf-8"
        ) as file:
            json.dump(summary, file, indent=2, ensure_ascii=False)

        LOGGER.info("SCENIC+ terminó su pipeline.")

        return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Exporta variantes candidatas y ejecuta SCENIC+ "
            "y CellOracle con sus datos reales."
        )
    )

    parser.add_argument(
        "--output-dir",
        default="resultados_scenic_celloracle",
        help="Carpeta de resultados."
    )

    parser.add_argument(
        "--half-width",
        type=int,
        default=250,
        help="Semiancho de las regiones candidatas en pares de bases."
    )

    parser.add_argument(
        "--celloracle-h5ad",
        help="Archivo AnnData H5AD con expresión celular."
    )

    parser.add_argument(
        "--celloracle-base-grn",
        help="Matriz base reguladora de CellOracle en CSV o Parquet."
    )

    parser.add_argument(
        "--cluster-column",
        default="cell_type",
        help="Columna de agrupación en adata.obs."
    )

    parser.add_argument(
        "--embedding-key",
        default="X_umap",
        help="Embedding existente en adata.obsm."
    )

    parser.add_argument(
        "--counts-layer",
        default=None,
        help="Capa de conteos brutos en adata.layers."
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Número máximo de componentes principales."
    )

    parser.add_argument(
        "--knn",
        type=int,
        default=30,
        help="Número de vecinos para imputación KNN."
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=10.0,
        help="Regularización para la inferencia de redes."
    )

    parser.add_argument(
        "--pvalue",
        type=float,
        default=0.001,
        help="Umbral p para filtrar las conexiones."
    )

    parser.add_argument(
        "--max-edges",
        type=int,
        default=2000,
        help="Número máximo de conexiones retenidas por filtro."
    )

    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Procesos usados por CellOracle."
    )

    parser.add_argument(
        "--scenicplus-config",
        help=(
            "Ruta a Snakemake/config/config.yaml del pipeline SCENIC+."
        )
    )

    parser.add_argument(
        "--scenicplus-cores",
        type=int,
        default=4,
        help="Núcleos de CPU para Snakemake."
    )

    parser.add_argument(
        "--simulate-gene",
        default=None,
        help=(
            "Gen cuya perturbación se simulará en CellOracle. "
            "Es una perturbación génica, no una edición de variante."
        )
    )

    parser.add_argument(
        "--simulate-value",
        type=float,
        default=0.0,
        help="Valor de expresión impuesto en la simulación."
    )

    parser.add_argument(
        "--grn-unit",
        default="whole",
        help=(
            "Unidad de red de CellOracle: 'whole' para la red global o el valor "
            "exacto de un grupo celular presente en adata.obs[cluster-column]."
        )
    )

    parser.add_argument(
        "--propagation-steps",
        type=int,
        default=3,
        help="Número de pasos de propagación en la simulación."
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if bool(args.celloracle_h5ad) != bool(args.celloracle_base_grn):
        parser.error(
            "--celloracle-h5ad y --celloracle-base-grn "
            "deben proporcionarse juntos."
        )

    if args.n_pcs < 1 or args.knn < 1:
        parser.error("--n-pcs y --knn deben ser positivos.")

    if args.alpha <= 0:
        parser.error("--alpha debe ser mayor que cero.")

    if not 0 < args.pvalue <= 1:
        parser.error("--pvalue debe estar entre 0 y 1.")

    if args.max_edges < 1 or args.scenicplus_cores < 1:
        parser.error("--max-edges y --scenicplus-cores deben ser positivos.")

    if args.propagation_steps < 1:
        parser.error("--propagation-steps debe ser positivo.")

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    adapter = SingleCellGRNAdapter(
        enhancer_half_width=args.half_width
    )

    try:
        annotations, bed_path, candidate_regions = (
            adapter.export_candidate_files(output_dir)
        )

        status = {
            "candidate_regions_exported": True,
            "celloracle": "not_requested",
            "scenicplus": "not_requested"
        }

        if args.celloracle_h5ad:
            adapter.run_celloracle(
                h5ad_path=Path(args.celloracle_h5ad).expanduser().resolve(),
                base_grn_path=(
                    Path(args.celloracle_base_grn).expanduser().resolve()
                ),
                annotations=annotations,
                output_dir=output_dir,
                cluster_column=args.cluster_column,
                embedding_key=args.embedding_key,
                counts_layer=args.counts_layer,
                n_pcs=args.n_pcs,
                knn=args.knn,
                alpha=args.alpha,
                pvalue=args.pvalue,
                max_edges=args.max_edges,
                n_jobs=args.n_jobs,
                simulate_gene=args.simulate_gene,
                simulate_value=args.simulate_value,
                grn_unit=args.grn_unit,
                propagation_steps=args.propagation_steps
            )

            status["celloracle"] = "completed"

        elif args.simulate_gene:
            parser.error(
                "--simulate-gene requiere los dos archivos de CellOracle."
            )

        if args.scenicplus_config:
            adapter.run_scenicplus(
                config_path=(
                    Path(args.scenicplus_config).expanduser().resolve()
                ),
                candidate_regions=candidate_regions,
                output_dir=output_dir,
                cores=args.scenicplus_cores
            )

            status["scenicplus"] = "completed"

        status_path = output_dir / "run_status.json"

        with status_path.open("w", encoding="utf-8") as file:
            json.dump(status, file, indent=2, ensure_ascii=False)

        if (
            status["celloracle"] == "not_requested"
            and status["scenicplus"] == "not_requested"
        ):
            LOGGER.warning(
                "Se exportaron las coordenadas, pero no se ejecutó "
                "ningún modelo. Proporciona las entradas de CellOracle "
                "o la configuración completa de SCENIC+."
            )
        else:
            LOGGER.info("Ejecución finalizada: %s", status)

        LOGGER.info("Resultados: %s", output_dir)
        return 0

    except Exception as exc:
        LOGGER.exception("Error durante la ejecución: %s", exc)
        status_path = output_dir / "run_status.json"
        try:
            if status_path.exists():
                with status_path.open("r", encoding="utf-8") as file:
                    saved_status = json.load(file)
            else:
                saved_status = {
                    "candidate_regions_exported": False,
                    "celloracle": "not_requested",
                    "scenicplus": "not_requested"
                }
            if args.celloracle_h5ad and saved_status.get("celloracle") != "completed":
                saved_status["celloracle"] = "failed"
            if args.scenicplus_config and saved_status.get("scenicplus") != "completed":
                saved_status["scenicplus"] = "failed"
            saved_status["error"] = f"{type(exc).__name__}: {exc}"
            with status_path.open("w", encoding="utf-8") as file:
                json.dump(saved_status, file, indent=2, ensure_ascii=False)
        except Exception:
            LOGGER.exception("No se pudo actualizar run_status.json tras el error.")
        return 1


if __name__ == "__main__":
    sys.exit(main())