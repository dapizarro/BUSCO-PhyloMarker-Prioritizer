#!/usr/bin/env python3
import os
import re
import argparse

def parse_busco_summary(summary_path):
    """
    Parse a BUSCO short_summary file and return a dict with key metrics.
    """
    data = {
        "busco_version": None,
        "lineage": None,
        "input_file": None,
        "C": None,
        "S": None,
        "D": None,
        "F": None,
        "M": None,
        "n": None,
    }

    with open(summary_path, "r") as fh:
        for line in fh:
            line = line.strip()

            # BUSCO version
            if line.startswith("# BUSCO version is:"):
                data["busco_version"] = line.split(":", 1)[1].strip()

            # Lineage
            elif line.startswith("# The lineage dataset is:"):
                lineage = line.split(":", 1)[1].strip()
                lineage = lineage.split("(", 1)[0].strip()
                data["lineage"] = lineage

            # Input file
            elif line.startswith("# Summarized benchmarking in BUSCO notation for file:"):
                data["input_file"] = line.split(":", 1)[1].strip()

            # Línea con los porcentajes C,S,D,F,M,n
            elif line.startswith("C:") and "n:" in line:
                clean = line.replace(",", " ").replace("[", " ").replace("]", " ")
                m = re.search(
                    r"C:(?P<C>[\d\.]+)%\s*S:(?P<S>[\d\.]+)%\s*D:(?P<D>[\d\.]+)%\s*F:(?P<F>[\d\.]+)%\s*M:(?P<M>[\d\.]+)%\s*n:(?P<n>\d+)",
                    clean
                )
                if m:
                    data["C"] = float(m.group("C"))
                    data["S"] = float(m.group("S"))
                    data["D"] = float(m.group("D"))
                    data["F"] = float(m.group("F"))
                    data["M"] = float(m.group("M"))
                    data["n"] = int(m.group("n"))

    return data


def main():
    parser = argparse.ArgumentParser(
        description="Parse all BUSCO short_summary files in a directory tree and produce a summary TSV."
    )
    parser.add_argument(
        "-i", "--busco_dir", required=True,
        help="Directorio raíz con los resultados de BUSCO (contiene las carpetas busco_XXX)."
    )
    parser.add_argument(
        "-o", "--out_tsv", default="busco_summary_all.tsv",
        help="Nombre del archivo de salida TSV (default: busco_summary_all.tsv)."
    )

    args = parser.parse_args()
    busco_dir = os.path.abspath(args.busco_dir)

    # Buscar todos los short_summary dentro de busco_dir
    summary_files = []
    for root, dirs, files in os.walk(busco_dir):
        for f in files:
            if f.startswith("short_summary") and f.endswith(".txt"):
                full_path = os.path.join(root, f)
                # ⬇️ FILTRO NUEVO: solo aceptar los que están directamente bajo una carpeta "busco_XXX"
                parent_dir = os.path.basename(os.path.dirname(full_path))
                if not parent_dir.startswith("busco_"):
                    continue
                summary_files.append(full_path)

    if not summary_files:
        raise SystemExit(f"No se encontraron archivos 'short_summary*.txt' válidos en {busco_dir}")

    # Procesar cada short_summary
    rows = []
    for s in sorted(summary_files):
        data = parse_busco_summary(s)

        # run_dir es la carpeta busco_XXX
        run_dir = os.path.basename(os.path.dirname(s))  # p.ej. "busco_Arthonio_Alyxoria_varia"
        if run_dir.startswith("busco_"):
            sample_id = run_dir[len("busco_"):]
        else:
            sample_id = run_dir

        # Cálculo aproximado de conteos absolutos
        n = data["n"] if data["n"] is not None else 0
        if n and data["C"] is not None:
            complete_count = round(n * data["C"] / 100.0)
        else:
            complete_count = None

        if n and data["M"] is not None:
            missing_count = round(n * data["M"] / 100.0)
        else:
            missing_count = None

        row = {
            "run_dir": run_dir,
            "sample_id": sample_id,
            "summary_file": s,
            "busco_version": data["busco_version"],
            "lineage": data["lineage"],
            "input_file": data["input_file"],
            "C_pct": data["C"],
            "S_pct": data["S"],
            "D_pct": data["D"],
            "F_pct": data["F"],
            "M_pct": data["M"],
            "n_busco": n,
            "n_complete_approx": complete_count,
            "n_missing_approx": missing_count,
        }
        rows.append(row)

    # Escribir la tabla TSV
    header = [
        "run_dir",
        "sample_id",
        "summary_file",
        "busco_version",
        "lineage",
        "input_file",
        "C_pct",
        "S_pct",
        "D_pct",
        "F_pct",
        "M_pct",
        "n_busco",
        "n_complete_approx",
        "n_missing_approx",
    ]

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(
                ["" if r[h] is None else str(r[h]) for h in header]
            ) + "\n")

    print(f"Escrito resumen de {len(rows)} runs en: {args.out_tsv}")


if __name__ == "__main__":
    main()

