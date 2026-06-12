#!/usr/bin/env python3
"""
find_optimal_panel_size.py
--------------------------

Goal
----
Derive a data-driven panel size (N loci) from cumulative improvements along a ranked gene list,
using a consensus (median) breakpoint across multiple alignment/phylogenetic-informativeness metrics.

Inputs
------
A TSV with per-gene metrics + ranking columns (typically produced by rank_genes_pca.py), e.g.:
  results/phylo_informativeness/per_gene_metrics_PCA.tsv

Required columns (minimum)
-------------------------
- gene
- rank  (1..n; smaller is better)
- occupancy
- gap_fraction
- PIS_per_len
- PIS_over_VS
- mean_entropy

Method (high level)
-------------------
1) Sort genes by rank (best -> worst).
2) For each N=1..maxN:
   - compute expanding (cumulative) MEDIAN for each metric over topN genes.
3) Normalize each cumulative curve to 0..1 across N (min-max on the curve).
   - gap_fraction is inverted (1 - normalized gap) because lower gaps are better.
4) For each normalized curve y(N), estimate a 2-segment breakpoint tau that minimizes SSE:
   y ~ a1 + b1*N  for N <= tau
   y ~ a2 + b2*N  for N >  tau
5) Return tau per metric, then consensus tau as the median across metrics.

Outputs
-------
- out_prefix + "_optimal_panel_size.tsv" : tau per metric + consensus
- out_prefix + "_cumulative_curves.tsv"  : cumulative medians + normalized curves + score
- out_prefix + "_optimal_panel_curve.png": plot of normalized curves + consensus tau

Usage
-----
python scripts/find_optimal_panel_size.py \
  -i results/phylo_informativeness/per_gene_metrics_PCA.tsv \
  -o results/phylo_informativeness/optimal_panel

Common options
--------------
--maxN 300        # evaluate N from 1..300 (or up to total genes)
--minN 10         # avoid extremely small N for breakpoint search
--edge_frac 0.05  # exclude breakpoint candidates too close to edges
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


METRICS = ["occupancy", "gap_fraction", "PIS_per_len", "PIS_over_VS", "mean_entropy"]


def _linreg_sse(x: np.ndarray, y: np.ndarray) -> float:
    """Return SSE of linear regression y ~ a + b*x for arrays x,y."""
    if len(x) < 2:
        return np.inf
    x = x.astype(float)
    y = y.astype(float)
    x_mean = x.mean()
    y_mean = y.mean()
    denom = np.sum((x - x_mean) ** 2)
    if denom == 0:
        return np.inf
    b = np.sum((x - x_mean) * (y - y_mean)) / denom
    a = y_mean - b * x_mean
    y_hat = a + b * x
    sse = np.sum((y - y_hat) ** 2)
    return float(sse)


def best_breakpoint(x: np.ndarray, y: np.ndarray, minN: int, edge_frac: float) -> int:
    """
    Find breakpoint tau (integer N) that minimizes SSE of two linear fits.
    Search tau in [minN .. maxN], excluding edges by edge_frac.
    """
    n = len(x)
    if n < max(20, minN * 2):
        # dataset too small to be meaningful
        return int(x[min(len(x) - 1, minN - 1)])

    # candidate taus (indices in x, using N values)
    lo = max(minN, int(np.ceil(edge_frac * n)))
    hi = min(n - minN, int(np.floor((1.0 - edge_frac) * n)))
    if lo >= hi:
        lo = minN
        hi = n - minN
    lo = max(2, lo)
    hi = max(lo + 1, hi)

    best_tau = int(x[lo - 1])
    best_sse = np.inf

    for tau_idx in range(lo, hi):
        # split: [0:tau_idx] and [tau_idx:n]
        sse1 = _linreg_sse(x[:tau_idx], y[:tau_idx])
        sse2 = _linreg_sse(x[tau_idx:], y[tau_idx:])
        sse = sse1 + sse2
        if sse < best_sse:
            best_sse = sse
            best_tau = int(x[tau_idx - 1])

    return best_tau


def minmax_norm(arr: np.ndarray) -> np.ndarray:
    """Min-max normalize to [0,1]. If constant, return zeros."""
    arr = np.asarray(arr, dtype=float)
    mn = np.nanmin(arr)
    mx = np.nanmax(arr)
    if not np.isfinite(mn) or not np.isfinite(mx) or mx == mn:
        return np.zeros_like(arr)
    return (arr - mn) / (mx - mn)


def main():
    ap = argparse.ArgumentParser(description="Find data-driven optimal panel size using consensus breakpoints across metrics.")
    ap.add_argument("-i", "--in_tsv", required=True, help="Input TSV with ranking + per-gene metrics (e.g., per_gene_metrics_PCA.tsv)")
    ap.add_argument("-o", "--out_prefix", required=True, help="Output prefix (dir/prefix), e.g. results/phylo_informativeness/optimal_panel")
    ap.add_argument("--maxN", type=int, default=300, help="Max panel size N to evaluate (default: 300)")
    ap.add_argument("--minN", type=int, default=10, help="Min N for breakpoint search (default: 10)")
    ap.add_argument("--edge_frac", type=float, default=0.05, help="Exclude breakpoint candidates within this fraction of edges (default: 0.05)")
    args = ap.parse_args()

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(args.in_tsv, sep="\t")
    missing = [c for c in (["gene", "rank"] + METRICS) if c not in df.columns]
    if missing:
        raise SystemExit(f"[ERROR] Missing required columns: {missing}")

    # Sort by rank (best first)
    df = df.sort_values("rank", ascending=True).reset_index(drop=True)

    total_genes = df.shape[0]
    maxN = min(args.maxN, total_genes)
    if maxN < args.minN:
        raise SystemExit(f"[ERROR] maxN={maxN} < minN={args.minN}. Increase maxN or lower minN.")

    # Use only top maxN for curve building (still derived from global ordering)
    df_top = df.iloc[:maxN].copy()

    # N axis
    N = np.arange(1, maxN + 1)

    # Expanding/cumulative medians for each metric
    curves = {"N": N}
    for m in METRICS:
        series = df_top[m].replace([np.inf, -np.inf], np.nan).astype(float)
        # expanding median is robust; min_periods=1 to start at N=1
        curves[m + "_cum_med"] = series.expanding(min_periods=1).median().values

    curves_df = pd.DataFrame(curves)

    # Normalize cumulative curves to 0..1
    norm = {}
    for m in METRICS:
        y = curves_df[m + "_cum_med"].values
        y_norm = minmax_norm(y)
        if m == "gap_fraction":
            # lower gaps are better -> invert after normalization
            y_norm = 1.0 - y_norm
        norm[m + "_norm"] = y_norm

    for k, v in norm.items():
        curves_df[k] = v

    # Composite score (equal weights across normalized metrics)
    norm_cols = [m + "_norm" for m in METRICS]
    curves_df["score"] = curves_df[norm_cols].mean(axis=1)

    # Breakpoints per metric and for score
    tau_rows = []
    for m in METRICS + ["score"]:
        y = curves_df[m + "_norm"].values if m != "score" else curves_df["score"].values
        tau = best_breakpoint(N, y, minN=args.minN, edge_frac=args.edge_frac)
        tau_rows.append({"metric": m, "tau": int(tau)})

    tau_df = pd.DataFrame(tau_rows)

    # Consensus tau: median across metrics (excluding score to avoid circularity)
    metric_taus = tau_df[tau_df["metric"].isin(METRICS)]["tau"].values
    tau_consensus = int(np.median(metric_taus))
    tau_df = pd.concat(
        [tau_df, pd.DataFrame([{"metric": "consensus_median", "tau": tau_consensus}])],
        ignore_index=True
    )

    # Write outputs
    curves_out = args.out_prefix + "_cumulative_curves.tsv"
    tau_out = args.out_prefix + "_optimal_panel_size.tsv"
    plot_out = args.out_prefix + "_optimal_panel_curve.png"

    curves_df.to_csv(curves_out, sep="\t", index=False)
    tau_df.to_csv(tau_out, sep="\t", index=False)

    # Plot: normalized curves + consensus
    plt.figure(figsize=(10, 6))
    for m in METRICS:
        plt.plot(curves_df["N"], curves_df[m + "_norm"], label=m)
    plt.plot(curves_df["N"], curves_df["score"], label="score", linewidth=2.5)

    plt.axvline(tau_consensus, linestyle="--", linewidth=2)
    plt.xlabel("Panel size (Top N loci along ranking)")
    plt.ylabel("Normalized cumulative metric (0–1)")
    plt.title(f"Data-driven panel size from saturation points (consensus tau = {tau_consensus})")
    plt.legend(loc="lower right", ncol=2, fontsize=9)
    plt.tight_layout()
    plt.savefig(plot_out, dpi=220)
    plt.close()

    print("[OK] Wrote:")
    print(" -", curves_out)
    print(" -", tau_out)
    print(" -", plot_out)
    print(f"[OK] Consensus panel size (median across metrics): N = {tau_consensus}")


if __name__ == "__main__":
    main()
