#!/usr/bin/env python3

"""
For single sequence use `python L7.py --input sequence.fasta --out_dir results`
For multiple sequences in the same file use `python L7.py --input sequence.fasta --input_all --out_dir results_multi`
For a directory containing multiple fasta files with one single sequence inside,
use `python L7.py --dir <directory containing multiple fasta files> --input_all --out_dir results_multi`
"""

import argparse
import os
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO


def read_fasta_single(path):
    recs = list(SeqIO.parse(path, "fasta"))
    if len(recs) == 0:
        raise ValueError(f"No sequences found in {path}")
    seq = str(recs[0].seq).upper()
    return seq.replace("\n", "").replace("\r", "")


def read_all_fasta_files(path):
    recs = list(SeqIO.parse(path, "fasta"))
    if len(recs) == 0:
        raise ValueError(f"No sequences found in {path}")
    seqs = [str(r.seq).upper().replace("\n", "").replace("\r", "") for r in recs]
    return seqs


def count_substrings(seq, min_l=2, max_l=6, min_repeat_units=2):
    seq = seq.upper()
    n = len(seq)
    counts = defaultdict(int)

    # For each unit length, scan sequence using regex to find maximal runs of repeated unit
    for L in range(min_l, max_l + 1):
        # slide start position to allow overlapping tandem repeats with different frame
        for i in range(0, n - L * min_repeat_units + 1):
            unit = seq[i:i + L]
            if set(unit) - {"A", "C", "G", "T"}:
                continue
            k = 1
            j = i + L
            while j + L <= n and seq[j:j + L] == unit:
                k += 1
                j += L
            if k >= min_repeat_units:
                counts[unit] += k
                seq = seq[:i] + ("N" * (k * L)) + seq[j:]
                n = len(seq)
    return dict(counts)


def filter_repeats(counts, min_count=2):
    # Only keep repeats that occur at least twice (i.e., are repeated)
    return {k: v for k, v in counts.items() if v >= min_count}


def top_n_counts(counts_dict, n=50):
    return dict(Counter(counts_dict).most_common(n))


def plot_top_bar(counts, title, outpath, top_n=30):
    top = top_n_counts(counts, top_n)
    if len(top) == 0:
        print(f"No repeats >=2 found for {title}")
        return
    keys = list(top.keys())
    vals = list(top.values())
    plt.figure(figsize=(max(6, len(keys) * 0.25 * 2), 4))
    plt.bar(range(len(keys)), vals, color="tab:blue")
    plt.xticks(range(len(keys)), keys, rotation=90, fontsize=8)
    plt.ylabel("Count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def analyze_file(path, out_dir, label=None, top_n_plot=30):
    if label is None:
        label = os.path.splitext(os.path.basename(path))[0]
    seq = read_fasta_single(path)
    if not (100 <= len(seq) <= 3000):
        print(f"Warning: sequence length {len(seq)} for {label} is outside 100-3000")
    counts = count_substrings(seq)
    repeats = filter_repeats(counts, min_count=2)
    os.makedirs(out_dir, exist_ok=True)
    plotp = os.path.join(out_dir, f"{label}_top_repeats.png")
    plot_top_bar(repeats, f"Top repeats in {label}", plotp, top_n=top_n_plot)
    return repeats


def compare_genomes(repeats_per_genome, out_dir, top_k=50):
    # Construct dataframe of top_k repeats across all genomes
    all_counts = Counter()
    for g, d in repeats_per_genome.items():
        all_counts.update(d)
    top_repeats = [r for r, _ in all_counts.most_common(top_k)]
    df = pd.DataFrame(index=top_repeats)
    for g, d in repeats_per_genome.items():
        df[g] = [d.get(r, 0) for r in top_repeats]
    df = df.fillna(0).astype(int)
    # Heatmap-like plot
    plt.figure(figsize=(max(6, len(top_repeats) * 0.15), max(4, len(repeats_per_genome) * 0.4)))
    plt.imshow(df.T, aspect='auto', cmap='viridis')
    plt.colorbar(label='count')
    plt.yticks(range(len(df.columns)), df.columns)
    plt.xticks(range(len(df.index)), df.index, rotation=90, fontsize=8)
    plt.title("Repeat counts (rows=genomes, cols=repeats)")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "comparison_heatmap.png"), dpi=150)
    plt.close()


def main():
    p = argparse.ArgumentParser(description="Detect and plot repeats (2-6 nt) in FASTA sequences")
    p.add_argument("--input", help="Single FASTA file path")
    p.add_argument("--input_all", action="store_true",
                   help="When used with --input, process all records in that file (default is only the first record)")
    p.add_argument("--dir", help="Directory containing FASTA files (process all records in each file)")
    p.add_argument("--out_dir", default="results", help="Output directory")
    p.add_argument("--top_plot", type=int, default=30, help="Top-N repeats to plot per genome")
    p.add_argument("--compare_top", type=int, default=50, help="Top-K repeats for comparison plot")
    args = p.parse_args()

    files = []
    mode = None  # "input", "input_all", "dir"
    if args.input:
        files = [args.input]
        mode = "single_all" if args.input_all else "single_first"
    elif args.dir:
        for fname in sorted(os.listdir(args.dir)):
            if fname.lower().endswith((".fa", ".fasta", ".fna")):
                files.append(os.path.join(args.dir, fname))
        mode = "dir"
    else:
        p.error("Either --input or --dir must be provided")

    if len(files) == 0:
        p.error("No FASTA files found")

    os.makedirs(args.out_dir, exist_ok=True)
    repeats_per_genome = {}

    for f in files:
        basename = os.path.splitext(os.path.basename(f))[0]

        if mode == "single_first":
            # read only the first record from the provided FASTA file
            try:
                seq = read_fasta_single(f)
            except Exception as e:
                print(f"Error reading {f}: {e}")
                continue
            label = f"{basename}_rec1"
            print("Processing (single-first)", label)
            if not (800 <= len(seq) <= 3000):
                print(f"Warning: sequence length {len(seq)} for {label} is outside 800-3000")
            counts = count_substrings(seq)
            repeats = filter_repeats(counts, min_count=2)
            plot_top_bar(repeats, f"Top repeats in {label}", os.path.join(args.out_dir, f"{label}_top_repeats.png"),
                         top_n=args.top_plot)
            repeats_per_genome[label] = repeats

        else:
            # For "single_all" and "dir" we process all records in each file
            try:
                recs = list(SeqIO.parse(f, "fasta"))
            except Exception as e:
                print(f"Error reading {f}: {e}")
                continue
            if len(recs) == 0:
                print(f"No records found in {f}, skipping.")
                continue
            for idx, rec in enumerate(recs, start=1):
                seq = str(rec.seq).upper().replace("\n", "").replace("\r", "")
                rec_id = rec.id if hasattr(rec, "id") and rec.id else f"rec{idx}"
                label = f"{basename}_{rec_id}"
                print("Processing", label)
                if not (800 <= len(seq) <= 3000):
                    print(f"Warning: sequence length {len(seq)} for {label} is outside 800-3000")
                counts = count_substrings(seq)
                repeats = filter_repeats(counts, min_count=2)
                plot_top_bar(repeats, f"Top repeats in {label}", os.path.join(args.out_dir, f"{label}_top_repeats.png"),
                             top_n=args.top_plot)
                repeats_per_genome[label] = repeats

    # If multiple genomes/records processed, produce comparison plot
    if len(repeats_per_genome) > 1:
        compare_genomes(repeats_per_genome, args.out_dir, top_k=args.compare_top)


if __name__ == "__main__":
    main()
