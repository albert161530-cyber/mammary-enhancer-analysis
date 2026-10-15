# ubuntu_bigwig_worker.py
"""
WORKER EXCLUSIVO PARA CONTENEDOR UBUNTU (Conda pyBigWig).
No requiere pandas ni matplotlib; solo pyBigWig y librería estándar de Python.
Uso en Ubuntu:
    python ubuntu_bigwig_worker.py "<BIGWIG_URL>" "<OUTPUT_CSV>" [WINDOW_BP]
"""

import sys
import csv
import pyBigWig

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

def stream_signals(bigwig_url: str, output_csv: str, window_bp: int = 150):
    print(f"  [Ubuntu Worker] Abriendo stream HTTP Range con pyBigWig: {bigwig_url}")
    bw = pyBigWig.open(bigwig_url)
    chroms_available = bw.chroms()

    rows = []
    for var_id, chrom, pos, gene in DRIVER_LOCI:
        start = max(0, pos - window_bp)
        end = pos + window_bp
        if chrom in chroms_available:
            mean_val = bw.stats(chrom, start, end, type="mean")[0]
            max_val = bw.stats(chrom, start, end, type="max")[0]
            rows.append({
                "Variant_ID": var_id,
                "Locus": gene,
                "Chrom": chrom,
                "Window": f"{start}-{end}",
                "BigWig_Mean_Signal": round(float(mean_val or 0.0), 4),
                "BigWig_Max_Peak": round(float(max_val or 0.0), 4)
            })
    bw.close()

    with open(output_csv, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Variant_ID", "Locus", "Chrom", "Window", "BigWig_Mean_Signal", "BigWig_Max_Peak"
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"  [Ubuntu Worker] Señales guardadas exitosamente en: {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python ubuntu_bigwig_worker.py <BIGWIG_URL> <OUTPUT_CSV> [WINDOW_BP]")
        sys.exit(1)
    url_arg = sys.argv[1]
    out_arg = sys.argv[2]
    win_arg = int(sys.argv[3]) if len(sys.argv) > 3 else 150
    stream_signals(url_arg, out_arg, win_arg)