#!/usr/bin/env bash
set -euo pipefail

BUSCO_DIR=${1:-data/busco_runs}
RESULTS=${2:-results}
THREADS=${THREADS:-1}

mkdir -p "$RESULTS" "$RESULTS/plots" "$RESULTS/panels"

echo "[1/7] Summarising BUSCO runs"
python scripts/make_busco_summary_table.py \
  -i "$BUSCO_DIR" \
  -o "$RESULTS/busco_summary_all.tsv"

echo "[2/7] Extracting single-copy BUSCO sequences per gene"
python scripts/extract_singlecopy_busco_per_gene.py \
  -i "$BUSCO_DIR" \
  -o "$RESULTS/busco_per_gene" \
  --mode single_copy \
  --ext .faa \
  --header species \
  --overwrite

echo "[3/7] Aligning and trimming per-gene FASTA files"
python scripts/align_and_trim.py \
  -i "$RESULTS/busco_per_gene/*.faa" \
  -o "$RESULTS/alignment_per_gene_trimmed" \
  --threads "$THREADS" \
  --mafft_mode auto \
  --gt 0.8 \
  --cons 60

echo "[4/7] Computing alignment-level informativeness metrics"
python scripts/compute_alignment_metrics.py \
  -i "$RESULTS/alignment_per_gene_trimmed/*.trim.faa" \
  -o "$RESULTS/per_gene_metrics.tsv"

echo "[5/7] Ranking genes by PCA"
python scripts/rank_genes_pca.py \
  -i "$RESULTS/per_gene_metrics.tsv" \
  -o "$RESULTS/per_gene_metrics_PCA.tsv"

echo "[6/7] Finding data-driven optimal panel size"
python scripts/find_optimal_panel_size.py \
  -i "$RESULTS/per_gene_metrics_PCA.tsv" \
  -o "$RESULTS/optimal_panel" \
  --maxN 300 \
  --minN 10 \
  --edge_frac 0.05

echo "[7/7] Plotting metrics and exporting marker panels"
python scripts/plot_metrics.py \
  -i "$RESULTS/per_gene_metrics_PCA.tsv" \
  -o "$RESULTS/plots" \
  --optimal_panel_tsv "$RESULTS/optimal_panel_optimal_panel_size.tsv" \
  --make_corr

python scripts/export_panels.py \
  -i "$RESULTS/per_gene_metrics_PCA.tsv" \
  -o "$RESULTS/panels" \
  --sizes 10,25,50,100,164,200,300,500

echo "[OK] Pipeline finished. Main output: $RESULTS/per_gene_metrics_PCA.tsv"
