#!/usr/bin/env python3
import argparse
import os

import numpy as np
import pandas as pd

# HPC-safe matplotlib backend (no X server needed)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Metrics to consider if present in the per-gene table
METRICS = [
    "occupancy",
    "gap_fraction",
    "PIS_per_len",
    "PIS_over_VS",
    "mean_entropy",
    "VS_per_len",
    "score",
    "consensus_median",
]

TOPN_FIXED = [25, 50, 100, 200]


def _clean_series(s: pd.Series) -> pd.Series:
    return s.replace([np.inf, -np.inf], np.nan)


def save_hist(df: pd.DataFrame, col: str, out_png: str, bins: int = 50) -> None:
    x = _clean_series(df[col]).dropna().values
    if x.size == 0:
        return
    plt.figure()
    plt.hist(x, bins=bins)
    plt.xlabel(col)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()


def save_scatter(df: pd.DataFrame, xcol: str, ycol: str, out_png: str,
                 alpha: float = 0.25, s: int = 8) -> None:
    x = _clean_series(df[xcol])
    y = _clean_series(df[ycol])
    sub = pd.concat([x, y], axis=1).dropna()
    if sub.shape[0] == 0:
        return
    plt.figure()
    plt.scatter(sub.iloc[:, 0].values, sub.iloc[:, 1].values, alpha=alpha, s=s)
    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()


def save_boxplot(groups: list[tuple[str, pd.DataFrame]], metric: str, out_png: str) -> None:
    labels, data = [], []
    for name, gdf in groups:
        if metric not in gdf.columns:
            continue
        vals = _clean_series(gdf[metric]).dropna().values
        if vals.size == 0:
            continue
        labels.append(name)
        data.append(vals)

    if len(data) < 2:
        return

    plt.figure()
    plt.boxplot(data, labels=labels, showfliers=False)
    plt.ylabel(metric)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()


def corr_heatmap(df: pd.DataFrame, cols: list[str], out_png: str) -> None:
    sub = df[cols].replace([np.inf, -np.inf], np.nan).dropna()
    if sub.shape[0] < 3:
        return
    corr = sub.corr(numeric_only=True).values

    plt.figure()
    plt.imshow(corr, interpolation="nearest")
    plt.xticks(range(len(cols)), cols, rotation=45, ha="right")
    plt.yticks(range(len(cols)), cols)
    plt.colorbar(label="Pearson r")
    plt.tight_layout()
    plt.savefig(out_png, dpi=240)
    plt.close()


def read_optimal_N_from_metric_tau(optimal_panel_tsv: str) -> int:
    """
    Read optimal N from a TSV with columns: metric, tau
    The optimal N is stored as tau for metric == 'consensus_median'.

    Example:
    metric  tau
    ...
    consensus_median  164
    """
    df = pd.read_csv(optimal_panel_tsv, sep=None, engine="python")  # auto-detect separator
    df.columns = [c.strip() for c in df.columns]

    required = {"metric", "tau"}
    if not required.issubset(set(df.columns)):
        raise ValueError(
            f"[ERROR] '{optimal_panel_tsv}' must contain columns {sorted(required)}. "
            f"Found columns: {list(df.columns)}"
        )

    df["metric"] = df["metric"].astype(str).str.strip()
    hit = df.loc[df["metric"] == "consensus_median", "tau"]

    if hit.empty:
        raise ValueError(
            f"[ERROR] '{optimal_panel_tsv}' does not contain metric 'consensus_median'. "
            f"Available metrics: {sorted(df['metric'].unique().tolist())}"
        )

    optimalN = int(float(hit.iloc[0]))
    return optimalN


def save_tau_barplot(optimal_panel_tsv: str, out_png: str) -> None:
    """
    Plot tau by metric from a TSV: metric, tau
    """
    df = pd.read_csv(optimal_panel_tsv, sep=None, engine="python")
    df.columns = [c.strip() for c in df.columns]
    if not {"metric", "tau"}.issubset(df.columns):
        return

    df["metric"] = df["metric"].astype(str).str.strip()
    df["tau"] = pd.to_numeric(df["tau"], errors="coerce")
    df = df.dropna(subset=["tau"]).sort_values("tau", ascending=False)

    if df.shape[0] == 0:
        return

    plt.figure(figsize=(9, max(3.5, 0.35 * len(df))))
    plt.barh(df["metric"].values, df["tau"].values)
    plt.xlabel("tau")
    plt.ylabel("metric")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(out_png, dpi=240)
    plt.close()


def build_groups(df: pd.DataFrame, topNs: list[int], optimalN: int | None):
    """
    Create groups: All, Top25/50/100/200, Optimal (N=optimalN)
    Uses 'rank' column if present. Otherwise uses head(n) after sorting by PC1 desc if present.
    """
    groups = [("All", df)]

    if "rank" in df.columns:
        df_sorted = df.sort_values("rank", ascending=True)
        for n in topNs:
            groups.append((f"Top{n}", df_sorted[df_sorted["rank"] <= n]))
        if optimalN is not None:
            groups.append((f"Optimal (N={optimalN})", df_sorted[df_sorted["rank"] <= optimalN]))
        return groups

    # Fallback if no rank: sort by PC1 descending if present, otherwise keep as-is.
    if "PC1" in df.columns:
        df_sorted = df.sort_values("PC1", ascending=False)
    else:
        df_sorted = df.copy()

    for n in topNs:
        groups.append((f"Top{n}", df_sorted.head(n)))
    if optimalN is not None:
        groups.append((f"Optimal (N={optimalN})", df_sorted.head(optimalN)))
    return groups


def main():
    ap = argparse.ArgumentParser(
        description="Generate stats + plots for per-gene metrics/PCA ranking. "
                    "Includes fixed TopN panels and an Optimal panel inferred from a metric/tau TSV."
    )
    ap.add_argument("-i", "--in_tsv", required=True, help="Input TSV (e.g., per_gene_metrics_PCA.tsv)")
    ap.add_argument("-o", "--out_dir", required=True, help="Output directory for plots and tables")
    ap.add_argument(
        "--optimal_panel_tsv",
        default=None,
        help="TSV with columns: metric, tau. Optimal N read as tau where metric==consensus_median."
    )
    ap.add_argument("--make_corr", action="store_true", help="Also generate a correlation heatmap")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.in_tsv, sep="\t")

    # Infer optimal N (if provided)
    optimalN = None
    if args.optimal_panel_tsv:
        optimalN = read_optimal_N_from_metric_tau(args.optimal_panel_tsv)
        print(f"[INFO] Optimal N inferred from metric/tau file (consensus_median): N={optimalN}")

        # Also save a tau-by-metric plot (nice for the report)
        save_tau_barplot(args.optimal_panel_tsv, os.path.join(args.out_dir, "tau_by_metric.png"))

    # Descriptive stats (for columns present)
    base_cols = ["PC1", "PC2", "rank", "aln_len", "n_taxa"]
    stats_cols = [c for c in (base_cols + METRICS) if c in df.columns]
    if stats_cols:
        desc = df[stats_cols].replace([np.inf, -np.inf], np.nan).describe(
            percentiles=[.01, .05, .1, .25, .5, .75, .9, .95, .99]
        ).T
        desc.to_csv(os.path.join(args.out_dir, "metrics_descriptive_stats.tsv"), sep="\t")

    # Histograms
    for col in ["occupancy", "gap_fraction", "PIS_per_len", "PIS_over_VS", "mean_entropy",
                "score", "consensus_median", "PC1"]:
        if col in df.columns:
            save_hist(df, col, os.path.join(args.out_dir, f"hist_{col}.png"))

    # Scatter: PC1 vs metrics
    if "PC1" in df.columns:
        for m in ["occupancy", "gap_fraction", "PIS_per_len", "PIS_over_VS",
                  "mean_entropy", "score", "consensus_median"]:
            if m in df.columns:
                save_scatter(df, "PC1", m, os.path.join(args.out_dir, f"scatter_PC1_vs_{m}.png"))

    # Groups: All + TopN + Optimal
    groups = build_groups(df, topNs=TOPN_FIXED, optimalN=optimalN)

    # Boxplots comparing All vs TopN vs Optimal
    for m in ["occupancy", "gap_fraction", "PIS_per_len", "PIS_over_VS",
              "mean_entropy", "score", "consensus_median"]:
        if m in df.columns:
            save_boxplot(groups, m, os.path.join(args.out_dir, f"box_{m}_All_vs_TopN_and_Optimal.png"))

    # Optional correlation heatmap
    if args.make_corr:
        cols = [c for c in (["PC1"] + METRICS) if c in df.columns]
        if len(cols) >= 3:
            corr_heatmap(df, cols, os.path.join(args.out_dir, "corr_heatmap.png"))

    # Save a tiny provenance file
    with open(os.path.join(args.out_dir, "panel_sizes_used.txt"), "w") as fh:
        fh.write("TopN_fixed:\t" + ",".join(map(str, TOPN_FIXED)) + "\n")
        fh.write("OptimalN:\t" + (str(optimalN) if optimalN is not None else "NA") + "\n")

    print(f"[OK] Wrote plots + tables to: {args.out_dir}")


if __name__ == "__main__":
    main()

