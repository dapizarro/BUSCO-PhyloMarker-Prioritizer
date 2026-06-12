#cd /home/dapizarr/pezizo_pipeline
#mkdir -p scripts results/phylo_informativeness

#cat > scripts/compute_alignment_metrics.py << 'PY'
#!/usr/bin/env python3
import math, glob, argparse, os
from collections import Counter

def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs

def shannon_entropy(states):
    n = len(states)
    if n == 0:
        return 0.0
    freqs = Counter(states)
    h = 0.0
    for c in freqs.values():
        p = c / n
        h -= p * math.log(p)
    return h

def is_gap(x):
    return x in ("-", ".")

def is_unknown(x):
    return x in ("X", "B", "Z", "J", "*", "?")

def compute_metrics(seqs, treat_unknown_as_missing=True):
    if not seqs:
        return None
    taxa = list(seqs.keys())
    L = len(next(iter(seqs.values())))
    for t, s in seqs.items():
        if len(s) != L:
            raise ValueError(f"Different alignment lengths in {t} (expected {L}, got {len(s)}).")

    n_taxa = len(taxa)
    var_sites = 0
    pis_sites = 0
    ent_sum = 0.0
    ent_cols = 0
    gap_cols = 0
    gap_chars = 0
    non_missing_cells = 0

    for i in range(L):
        col = [seqs[t][i] for t in taxa]

        g = sum(1 for x in col if is_gap(x))
        gap_chars += g
        if g > 0:
            gap_cols += 1

        residues = []
        for x in col:
            if is_gap(x):
                continue
            if treat_unknown_as_missing and is_unknown(x):
                continue
            residues.append(x)

        non_missing_cells += len(residues)

        if len(residues) >= 2:
            states = Counter(residues)
            if len(states) >= 2:
                var_sites += 1
                if sum(1 for c in states.values() if c >= 2) >= 2:
                    pis_sites += 1

            ent_sum += shannon_entropy(residues)
            ent_cols += 1

    gap_fraction = gap_chars / (n_taxa * L) if (n_taxa * L) else 0.0
    occupancy = non_missing_cells / (n_taxa * L) if (n_taxa * L) else 0.0
    mean_entropy = ent_sum / ent_cols if ent_cols else 0.0

    return {
        "n_taxa": n_taxa,
        "aln_len": L,
        "occupancy": occupancy,
        "gap_fraction": gap_fraction,
        "gap_cols_fraction": (gap_cols / L) if L else 0.0,
        "variable_sites": var_sites,
        "PIS": pis_sites,
        "VS_per_len": (var_sites / L) if L else 0.0,
        "PIS_per_len": (pis_sites / L) if L else 0.0,
        "PIS_over_VS": (pis_sites / var_sites) if var_sites else 0.0,
        "mean_entropy": mean_entropy,
    }

def main():
    ap = argparse.ArgumentParser(description="Compute per-alignment phylo-informativeness metrics from protein FASTA alignments.")
    ap.add_argument("-i", "--input_glob", required=True, help="e.g. 'results/alignment_per_gene_trimmed/*.trim.faa'")
    ap.add_argument("-o", "--out_tsv", required=True, help="Output TSV path")
    ap.add_argument("--keep_unknowns", action="store_true",
                    help="Count X/B/Z/J/?/* as residues (NOT recommended for scoring).")
    args = ap.parse_args()

    paths = sorted(glob.glob(args.input_glob))
    if not paths:
        raise SystemExit(f"No files matched: {args.input_glob}")

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)

    header = ["gene","n_taxa","aln_len","occupancy","gap_fraction","gap_cols_fraction",
              "variable_sites","PIS","VS_per_len","PIS_per_len","PIS_over_VS","mean_entropy"]

    with open(args.out_tsv, "w") as out:
        out.write("\t".join(header) + "\n")
        for p in paths:
            seqs = read_fasta(p)
            m = compute_metrics(seqs, treat_unknown_as_missing=(not args.keep_unknowns))
            gene = os.path.basename(p)
            row = [gene] + [str(m[h]) for h in header[1:]]
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()




