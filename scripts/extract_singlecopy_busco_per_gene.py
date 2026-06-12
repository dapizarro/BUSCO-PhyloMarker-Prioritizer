#!/usr/bin/env python3
import argparse
from pathlib import Path

def read_fasta(path: Path):
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)

def write_wrapped(out, s, width=60):
    for i in range(0, len(s), width):
        out.write(s[i:i+width] + "\n")

def find_seq_dir(species_dir: Path, mode: str) -> Path | None:
    # Busca species_dir/run_*/busco_sequences/<mode>_busco_sequences
    for run_dir in sorted(species_dir.glob("run_*")):
        candidate = run_dir / "busco_sequences" / f"{mode}_busco_sequences"
        if candidate.is_dir():
            return candidate
    return None

def main():
    ap = argparse.ArgumentParser(
        description="Extrae BUSCO sequences de busco_runs/busco_*/run_*/busco_sequences/... y crea FASTAs por gen sin abrir miles de archivos."
    )
    ap.add_argument("-i", "--busco_runs_dir", required=True,
                    help="Directorio con carpetas busco_* (p.ej. .../data/busco_runs).")
    ap.add_argument("-o", "--out_dir", required=True,
                    help="Salida: un FASTA por BUSCO ID.")
    ap.add_argument("--mode", default="single_copy", choices=["single_copy", "multi_copy", "fragmented"],
                    help="Tipo de secuencias BUSCO (default: single_copy).")
    ap.add_argument("--ext", default=".faa",
                    help="Extensión a extraer (.faa proteínas, .fna nucleótidos). Default: .faa")
    ap.add_argument("--header", default="species", choices=["species", "species|orig", "orig"],
                    help="Header salida: species / species|orig / orig. Default: species")
    ap.add_argument("--overwrite", action="store_true",
                    help="Si se activa, limpia out_dir antes (solo archivos con esa extensión).")
    ap.add_argument("--report", default="busco_extract_report.tsv",
                    help="Nombre del TSV de reporte dentro de out_dir.")
    args = ap.parse_args()

    busco_runs_dir = Path(args.busco_runs_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ext = args.ext if args.ext.startswith(".") else f".{args.ext}"

    species_dirs = [p for p in sorted(busco_runs_dir.iterdir()) if p.is_dir() and p.name.startswith("busco_")]
    if not species_dirs:
        raise SystemExit(f"[ERROR] No se encontraron carpetas busco_* en: {busco_runs_dir}")

    # Si overwrite: borramos FASTAs per-gene previos (solo los de la extensión elegida)
    if args.overwrite:
        for f in out_dir.glob(f"*{ext}"):
            try:
                f.unlink()
            except Exception:
                pass

    n_species_total = len(species_dirs)
    n_species_ok = 0
    n_species_missing = 0
    n_records = 0
    genes_seen = set()

    report_lines = ["species\tstatus\tseq_dir\tn_files\tn_records\n"]

    for sp_dir in species_dirs:
        species = sp_dir.name[len("busco_"):]
        seq_dir = find_seq_dir(sp_dir, args.mode)

        if seq_dir is None:
            n_species_missing += 1
            report_lines.append(f"{species}\tMISSING\t\t0\t0\n")
            continue

        files = [f for f in sorted(seq_dir.iterdir()) if f.is_file() and f.suffix.lower() == ext.lower()]
        if not files:
            n_species_missing += 1
            report_lines.append(f"{species}\tNOFILES\t{seq_dir}\t0\t0\n")
            continue

        sp_records = 0
        for f in files:
            gene_id = f.stem
            genes_seen.add(gene_id)
            out_path = out_dir / f"{gene_id}{ext}"

            # Abrir SOLO este output, escribir, cerrar → evita "too many open files"
            with open(out_path, "a") as out:
                for orig_header, seq in read_fasta(f):
                    if args.header == "species":
                        hdr = species
                    elif args.header == "species|orig":
                        hdr = f"{species}|{orig_header}"
                    else:
                        hdr = orig_header

                    out.write(f">{hdr}\n")
                    write_wrapped(out, seq)
                    n_records += 1
                    sp_records += 1

        n_species_ok += 1
        report_lines.append(f"{species}\tOK\t{seq_dir}\t{len(files)}\t{sp_records}\n")

    report_path = out_dir / args.report
    with open(report_path, "w") as rep:
        rep.writelines(report_lines)

    print("Resumen")
    print("-------")
    print(f"busco_* encontrados: {n_species_total}")
    print(f"Especies OK:        {n_species_ok}")
    print(f"Especies faltan:    {n_species_missing}")
    print(f"Genes detectados:   {len(genes_seen)}")
    print(f"Registros escritos: {n_records}")
    print(f"Salida:            {out_dir}")
    print(f"Reporte:           {report_path}")

if __name__ == "__main__":
    main()

