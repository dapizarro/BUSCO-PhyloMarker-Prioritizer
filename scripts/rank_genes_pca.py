#!/usr/bin/env python3
import argparse
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def main():
    ap = argparse.ArgumentParser(description="Rank genes using PCA on phylogenetic informativeness metrics.")
    ap.add_argument("-i", "--in_tsv", required=True, help="Input TSV with per-gene metrics")
    ap.add_argument("-o", "--out_tsv", required=True, help="Output TSV with PCA scores and ranking")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t")

    metrics = [
        "occupancy",
        "VS_per_len",
        "PIS_per_len",
        "PIS_over_VS",
        "mean_entropy"
    ]

    X = df[metrics].fillna(0.0)

    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X_scaled)

    df["PC1"] = pcs[:, 0]
    df["PC2"] = pcs[:, 1]

    df["rank"] = df["PC1"].rank(ascending=False, method="dense").astype(int)
    df = df.sort_values("PC1", ascending=False)

    df.to_csv(args.out_tsv, sep="\t", index=False)

    print("Explained variance:")
    print(f"PC1: {pca.explained_variance_ratio_[0]:.3f}")
    print(f"PC2: {pca.explained_variance_ratio_[1]:.3f}")

if __name__ == "__main__":
    main()

