"""Worker BigWig para Ubuntu/WSL.

Calcula señal en los 10 loci candidatos, cuerpo de cada gen, flancos
izquierdo/derecho de 1 Mb y controles locales equivalentes por locus.
Las coordenadas de variantes se conservan. Intervalos BigWig: 0-based semiabiertos.
"""
import sys
import csv
import urllib.request
import urllib.parse
from pathlib import Path
import pyBigWig

WINDOW_BP = 150
DISTAL_FLANK_BP = 1_000_000
LOCAL_CONTROL_OFFSET_BP = 10_000

DRIVER_LOCI = [
    ("chr6:150866313:T>C", "chr6", 150866313, "ESR1"),
    ("chr6:151452006:A>G", "chr6", 151452006, "ESR1"),
    ("chr6:151808589:A>G", "chr6", 151808589, "ESR1"),
    ("chr15:51569074:A>G", "chr15", 51569074, "CYP19A1"),
    ("chr6:150866751:T>C", "chr6", 150866751, "ESR1"),
    ("chr14:64280081:A>G", "chr14", 64280081, "ESR2"),
    ("chr6:151299349:A>G", "chr6", 151299349, "ESR1"),
    ("chr6:151014628:A>G", "chr6", 151014628, "ESR1"),
    ("chr6:151322799:T>C", "chr6", 151322799, "ESR1"),
    ("chr6:151763572:G>C", "chr6", 151763572, "ESR1"),
]

# GRCh38; coordenadas 1-based inclusivas.
GENE_COORDS = {
    "ESR1": ("chr6", 151656672, 152129619),
    "CYP19A1": ("chr15", 51208057, 51338596),
    "ESR2": ("chr14", 64226707, 64338613),
}

LOCUS_FIELDS = [
    "Variant_ID", "Gene", "Chrom", "Locus_Window_0based",
    "Locus_Mean_Signal", "Locus_Max_Peak", "Gene_Mean_Signal", "Gene_Max_Peak",
    "Left_Distal_Mean_Signal", "Left_Distal_Max_Peak",
    "Right_Distal_Mean_Signal", "Right_Distal_Max_Peak",
    "Local_Left_Control_Window_0based", "Local_Left_Control_Mean", "Local_Left_Control_Max",
    "Local_Right_Control_Window_0based", "Local_Right_Control_Mean", "Local_Right_Control_Max",
    "Local_Control_Mean_Average", "Locus_vs_Gene_Mean_Ratio",
    "Locus_vs_Left_Distal_Mean_Ratio", "Locus_vs_Right_Distal_Mean_Ratio",
    "Locus_vs_Local_Left_Control_Ratio", "Locus_vs_Local_Right_Control_Ratio",
    "Locus_vs_Local_Control_Average_Ratio", "Locus_Region_Relative_To_Gene",
    "Locus_Within_Gene_Bounds",
]
GENE_FIELDS = [
    "Gene", "Chrom", "Gene_Window_0based", "Gene_Length_bp", "Gene_Mean_Signal", "Gene_Max_Peak",
    "Left_Distal_Window_0based", "Left_Distal_Length_bp", "Left_Distal_Mean_Signal", "Left_Distal_Max_Peak",
    "Right_Distal_Window_0based", "Right_Distal_Length_bp", "Right_Distal_Mean_Signal", "Right_Distal_Max_Peak",
]
COMPARISON_FIELDS = [
    "Variant_ID", "Gene", "Locus_Mean_Signal", "Gene_Mean_Signal", "Left_Distal_Mean_Signal",
    "Right_Distal_Mean_Signal", "Local_Left_Control_Mean", "Local_Right_Control_Mean",
    "Local_Control_Mean_Average", "Locus_vs_Gene_Mean_Ratio", "Locus_vs_Left_Distal_Mean_Ratio",
    "Locus_vs_Right_Distal_Mean_Ratio", "Locus_vs_Local_Left_Control_Ratio",
    "Locus_vs_Local_Right_Control_Ratio", "Locus_vs_Local_Control_Average_Ratio", "Locus_Max_Peak",
    "Gene_Max_Peak", "Left_Distal_Max_Peak", "Right_Distal_Max_Peak", "Local_Left_Control_Max",
    "Local_Right_Control_Max", "Locus_Region_Relative_To_Gene", "Locus_Within_Gene_Bounds",
]


def resolve_bigwig_url(url):
    """Resuelve la redirección de ENCODE sin descargar el BigWig entero."""
    request = urllib.request.Request(url, headers={
        "User-Agent": "BioBigWigWorker/1.0", "Range": "bytes=0-63"
    })
    print("[Worker] Resolviendo URL de ENCODE...", flush=True)
    with urllib.request.urlopen(request, timeout=60) as response:
        status = getattr(response, "status", 200)
        final_url = response.geturl()
        sample = response.read(64)
    if status not in (200, 206):
        raise RuntimeError(f"Respuesta HTTP inesperada: {status}")
    if len(sample) < 4 or sample[:4] != bytes.fromhex("26fc8f88"):
        raise RuntimeError(f"La respuesta no tiene firma BigWig. Bytes: {sample[:4].hex()}")
    print("[Worker] URL final:", urllib.parse.urlparse(final_url).netloc, flush=True)
    return final_url


def get_stats(bw, chrom, start0, end0):
    """Media y máximo exactos para un intervalo 0-based semiabierto."""
    if start0 is None or end0 is None or start0 >= end0:
        return None, None
    means = bw.stats(chrom, start0, end0, type="mean", exact=True)
    maxima = bw.stats(chrom, start0, end0, type="max", exact=True)
    return means[0] if means else None, maxima[0] if maxima else None


def rounded(value):
    return round(float(value), 4) if value is not None else None


def safe_ratio(numerator, denominator):
    if numerator is None or denominator is None or float(denominator) == 0:
        return None
    return rounded(float(numerator) / float(denominator))


def average_available(a, b):
    values = [float(v) for v in (a, b) if v is not None]
    return sum(values) / len(values) if values else None


def write_csv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[Worker] CSV guardado: {path}", flush=True)


def extract_signals(bigwig_url, output_csv, window_bp=WINDOW_BP,
                    local_offset_bp=LOCAL_CONTROL_OFFSET_BP):
    if window_bp <= 0:
        raise ValueError("window_bp debe ser positivo.")
    if local_offset_bp <= window_bp:
        raise ValueError("local_offset_bp debe ser mayor que window_bp para evitar solapamiento.")

    print("[Worker] Abriendo BigWig remoto...", flush=True)
    bw = pyBigWig.open(bigwig_url)
    if bw is None:
        raise RuntimeError("pyBigWig no pudo abrir el archivo.")

    locus_rows, gene_rows, comparison_rows = [], [], []
    gene_stats = {}
    try:
        available = bw.chroms()
        if not available:
            raise RuntimeError("No se encontraron cromosomas.")
        print(f"[Worker] Cromosomas disponibles: {len(available)}", flush=True)
        print(f"[Worker] Flancos distales: {DISTAL_FLANK_BP:,} bp por lado; controles locales a ±{local_offset_bp:,} bp", flush=True)

        # Cuerpo del gen y flancos genómicos de 1 Mb.
        for gene, (chrom, start1, end1) in GENE_COORDS.items():
            if chrom not in available:
                raise RuntimeError(f"Falta {chrom} en el BigWig para {gene}.")
            chrom_length = available[chrom]
            gene_start, gene_end = max(0, start1 - 1), min(end1, chrom_length)
            left_start, left_end = max(0, gene_start - DISTAL_FLANK_BP), gene_start
            right_start, right_end = gene_end, min(chrom_length, gene_end + DISTAL_FLANK_BP)
            gene_mean, gene_max = get_stats(bw, chrom, gene_start, gene_end)
            left_mean, left_max = get_stats(bw, chrom, left_start, left_end)
            right_mean, right_max = get_stats(bw, chrom, right_start, right_end)
            gene_stats[gene] = {
                "chrom": chrom, "start": gene_start, "end": gene_end,
                "mean": gene_mean, "max": gene_max,
                "left_mean": left_mean, "left_max": left_max,
                "right_mean": right_mean, "right_max": right_max,
            }
            gene_rows.append({
                "Gene": gene, "Chrom": chrom, "Gene_Window_0based": f"{gene_start}-{gene_end}",
                "Gene_Length_bp": gene_end - gene_start, "Gene_Mean_Signal": rounded(gene_mean),
                "Gene_Max_Peak": rounded(gene_max), "Left_Distal_Window_0based": f"{left_start}-{left_end}",
                "Left_Distal_Length_bp": left_end - left_start, "Left_Distal_Mean_Signal": rounded(left_mean),
                "Left_Distal_Max_Peak": rounded(left_max), "Right_Distal_Window_0based": f"{right_start}-{right_end}",
                "Right_Distal_Length_bp": right_end - right_start, "Right_Distal_Mean_Signal": rounded(right_mean),
                "Right_Distal_Max_Peak": rounded(right_max),
            })
            print(f"[Worker] {gene}: gen={rounded(gene_mean)}, izq={rounded(left_mean)}, der={rounded(right_mean)}", flush=True)

        # Cada locus se compara con controles simétricos de igual tamaño.
        for index, (var_id, chrom, pos1, gene) in enumerate(DRIVER_LOCI, start=1):
            if chrom not in available:
                raise RuntimeError(f"Falta {chrom} en el BigWig para {var_id}.")
            chrom_length = available[chrom]
            locus_start = max(0, pos1 - window_bp - 1)
            locus_end = min(pos1 + window_bp, chrom_length)
            locus_mean, locus_max = get_stats(bw, chrom, locus_start, locus_end)

            center0 = pos1 - 1
            left_center0, right_center0 = center0 - local_offset_bp, center0 + local_offset_bp
            # Los controles son exactamente 2*window_bp+1 bases. Si no caben completos, quedan vacíos.
            if left_center0 - window_bp < 0 or left_center0 + window_bp >= chrom_length:
                lc_start = lc_end = None
                lc_mean = lc_max = None
            else:
                lc_start, lc_end = left_center0 - window_bp, left_center0 + window_bp + 1
                lc_mean, lc_max = get_stats(bw, chrom, lc_start, lc_end)
            if right_center0 - window_bp < 0 or right_center0 + window_bp >= chrom_length:
                rc_start = rc_end = None
                rc_mean = rc_max = None
            else:
                rc_start, rc_end = right_center0 - window_bp, right_center0 + window_bp + 1
                rc_mean, rc_max = get_stats(bw, chrom, rc_start, rc_end)

            info = gene_stats[gene]
            inside = chrom == info["chrom"] and locus_start >= info["start"] and locus_end <= info["end"]
            if inside:
                region = "gene_body"
            elif locus_end <= info["start"]:
                region = "left_genomic_flank"
            elif locus_start >= info["end"]:
                region = "right_genomic_flank"
            else:
                region = "overlaps_gene_boundary"

            local_avg = average_available(lc_mean, rc_mean)
            gene_ratio = safe_ratio(locus_mean, info["mean"])
            left_ratio = safe_ratio(locus_mean, info["left_mean"])
            right_ratio = safe_ratio(locus_mean, info["right_mean"])
            lc_ratio, rc_ratio = safe_ratio(locus_mean, lc_mean), safe_ratio(locus_mean, rc_mean)
            local_ratio = safe_ratio(locus_mean, local_avg)
            row = {
                "Variant_ID": var_id, "Gene": gene, "Chrom": chrom,
                "Locus_Window_0based": f"{locus_start}-{locus_end}",
                "Locus_Mean_Signal": rounded(locus_mean), "Locus_Max_Peak": rounded(locus_max),
                "Gene_Mean_Signal": rounded(info["mean"]), "Gene_Max_Peak": rounded(info["max"]),
                "Left_Distal_Mean_Signal": rounded(info["left_mean"]), "Left_Distal_Max_Peak": rounded(info["left_max"]),
                "Right_Distal_Mean_Signal": rounded(info["right_mean"]), "Right_Distal_Max_Peak": rounded(info["right_max"]),
                "Local_Left_Control_Window_0based": f"{lc_start}-{lc_end}" if lc_start is not None else "",
                "Local_Left_Control_Mean": rounded(lc_mean), "Local_Left_Control_Max": rounded(lc_max),
                "Local_Right_Control_Window_0based": f"{rc_start}-{rc_end}" if rc_start is not None else "",
                "Local_Right_Control_Mean": rounded(rc_mean), "Local_Right_Control_Max": rounded(rc_max),
                "Local_Control_Mean_Average": rounded(local_avg), "Locus_vs_Gene_Mean_Ratio": gene_ratio,
                "Locus_vs_Left_Distal_Mean_Ratio": left_ratio, "Locus_vs_Right_Distal_Mean_Ratio": right_ratio,
                "Locus_vs_Local_Left_Control_Ratio": lc_ratio, "Locus_vs_Local_Right_Control_Ratio": rc_ratio,
                "Locus_vs_Local_Control_Average_Ratio": local_ratio,
                "Locus_Region_Relative_To_Gene": region, "Locus_Within_Gene_Bounds": inside,
            }
            locus_rows.append(row)
            comparison_rows.append({
                "Variant_ID": var_id, "Gene": gene, "Locus_Mean_Signal": rounded(locus_mean),
                "Gene_Mean_Signal": rounded(info["mean"]), "Left_Distal_Mean_Signal": rounded(info["left_mean"]),
                "Right_Distal_Mean_Signal": rounded(info["right_mean"]),
                "Local_Left_Control_Mean": rounded(lc_mean), "Local_Right_Control_Mean": rounded(rc_mean),
                "Local_Control_Mean_Average": rounded(local_avg), "Locus_vs_Gene_Mean_Ratio": gene_ratio,
                "Locus_vs_Left_Distal_Mean_Ratio": left_ratio, "Locus_vs_Right_Distal_Mean_Ratio": right_ratio,
                "Locus_vs_Local_Left_Control_Ratio": lc_ratio, "Locus_vs_Local_Right_Control_Ratio": rc_ratio,
                "Locus_vs_Local_Control_Average_Ratio": local_ratio, "Locus_Max_Peak": rounded(locus_max),
                "Gene_Max_Peak": rounded(info["max"]), "Left_Distal_Max_Peak": rounded(info["left_max"]),
                "Right_Distal_Max_Peak": rounded(info["right_max"]), "Local_Left_Control_Max": rounded(lc_max),
                "Local_Right_Control_Max": rounded(rc_max), "Locus_Region_Relative_To_Gene": region,
                "Locus_Within_Gene_Bounds": inside,
            })
            print(f"[Worker] Locus {index}/{len(DRIVER_LOCI)}: {var_id}; región={region}; señal={rounded(locus_mean)}; control_local={rounded(local_avg)}; ratio_local={local_ratio}", flush=True)
    finally:
        bw.close()

    output_path = Path(output_csv)
    write_csv(output_path, LOCUS_FIELDS, locus_rows)
    write_csv(output_path.with_name(f"{output_path.stem}_gene_summary.csv"), GENE_FIELDS, gene_rows)
    write_csv(output_path.with_name(f"{output_path.stem}_locus_vs_gene.csv"), COMPARISON_FIELDS, comparison_rows)
    print(f"[Worker] Loci procesados: {len(locus_rows)}/{len(DRIVER_LOCI)}", flush=True)


def stream_signals(bigwig_url, output_csv, window_bp=WINDOW_BP, local_offset_bp=LOCAL_CONTROL_OFFSET_BP):
    final_url = resolve_bigwig_url(bigwig_url)
    extract_signals(final_url, output_csv, window_bp, local_offset_bp)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python ubuntu_bigwig_worker.py <BIGWIG_URL> <OUTPUT_CSV> [WINDOW_BP] [LOCAL_CONTROL_OFFSET_BP]")
        sys.exit(1)
    try:
        url_arg, output_arg = sys.argv[1], sys.argv[2]
        window_arg = int(sys.argv[3]) if len(sys.argv) > 3 else WINDOW_BP
        offset_arg = int(sys.argv[4]) if len(sys.argv) > 4 else LOCAL_CONTROL_OFFSET_BP
        stream_signals(url_arg, output_arg, window_arg, offset_arg)
    except Exception as error:
        print(f"[Worker] ERROR: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
