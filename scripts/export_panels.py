#!/usr/bin/env python3
"""Export ranked marker panels from PCA/ranking table."""
import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description='Export top-N gene lists from a ranked TSV.')
    ap.add_argument('-i','--ranked_tsv', required=True)
    ap.add_argument('-o','--out_dir', required=True)
    ap.add_argument('--sizes', default='10,25,50,100,164,200,300,500', help='Comma-separated panel sizes')
    ap.add_argument('--gene_col', default='gene')
    ap.add_argument('--rank_col', default='rank')
    args = ap.parse_args()
    df = pd.read_csv(args.ranked_tsv, sep='\t').sort_values(args.rank_col)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    sizes = [int(x) for x in args.sizes.split(',') if x.strip()]
    for n in sizes:
        sub = df.head(n)
        sub[[args.gene_col,args.rank_col]].to_csv(out / f'top{n}_genes.tsv', sep='\t', index=False)
        with open(out / f'top{n}_genes.txt','w') as fh:
            fh.write('\n'.join(sub[args.gene_col].astype(str)) + '\n')
    print(f'[OK] Exported panels to {out}')

if __name__ == '__main__':
    main()
