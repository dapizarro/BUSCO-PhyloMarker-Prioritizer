#!/usr/bin/env python3
"""Align per-gene BUSCO FASTA files with MAFFT and trim with trimAl."""
import argparse, glob, os, subprocess
from pathlib import Path


def run(cmd):
    print('[RUN]', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description='MAFFT + trimAl for per-gene protein FASTA files.')
    ap.add_argument('-i','--input_glob', required=True, help="Input FASTA glob, e.g. results/busco_per_gene/*.faa")
    ap.add_argument('-o','--out_dir', required=True, help='Output directory for trimmed alignments')
    ap.add_argument('--threads', type=int, default=1, help='MAFFT threads per gene')
    ap.add_argument('--mafft_mode', default='auto', choices=['auto','linsi','localpair'], help='MAFFT strategy')
    ap.add_argument('--gt', default='0.8', help='trimAl -gt threshold')
    ap.add_argument('--cons', default='60', help='trimAl -cons threshold')
    ap.add_argument('--keep_untrimmed', action='store_true')
    args = ap.parse_args()

    paths = sorted(glob.glob(args.input_glob))
    if not paths:
        raise SystemExit(f'No FASTA files matched: {args.input_glob}')
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    aln_dir = out / 'untrimmed'; aln_dir.mkdir(exist_ok=True)

    for p in paths:
        p = Path(p)
        stem = p.stem
        aln = aln_dir / f'{stem}.aln.faa'
        trim = out / f'{stem}.trim.faa'
        if args.mafft_mode == 'linsi':
            mafft_cmd = ['mafft','--localpair','--maxiterate','1000','--thread',str(args.threads),str(p)]
        elif args.mafft_mode == 'localpair':
            mafft_cmd = ['mafft','--localpair','--thread',str(args.threads),str(p)]
        else:
            mafft_cmd = ['mafft','--auto','--thread',str(args.threads),str(p)]
        print(f'[INFO] Aligning {p.name}')
        with open(aln, 'w') as fh:
            subprocess.run(mafft_cmd, check=True, stdout=fh)
        run(['trimal','-in',str(aln),'-out',str(trim),'-gt',str(args.gt),'-cons',str(args.cons)])
        if not args.keep_untrimmed:
            aln.unlink(missing_ok=True)
    print(f'[OK] Trimmed alignments written to {out}')

if __name__ == '__main__':
    main()
