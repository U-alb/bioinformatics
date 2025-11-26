#!/usr/bin/env python3
"""
Lab9_1.py

Usage:
  python Lab9_1.py input.fasta

Outputs to stdout:
 - cleavage counts and positions for each enzyme
 - fragment lengths per enzyme
 - simple ASCII gel simulation for single and combined digests

Enzymes supported: EcoRI, BamHI, HindIII, TaqI, HaeIII
"""

import sys

ENZYMES = {
    "EcoRI":  {"site":"GAATTC", "cut": 1},  # cuts between G|AATTC -> coordinate = start + 1
    "BamHI":  {"site":"GGATCC", "cut": 1},  # G|GATCC
    "HindIII":{"site":"AAGCTT", "cut": 1},  # A|AGCTT
    "TaqI":   {"site":"TCGA",   "cut": 1},  # T|CGA
    "HaeIII": {"site":"GGCC",   "cut": 2},  # GG|CC (cuts in middle) -> coordinate = start + 2
}

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"): continue
            seq.append(line.strip())
    return "".join(seq).upper()

def find_cuts(seq, site, cut_offset):
    starts = []
    i = seq.find(site, 0)
    while i != -1:
        cut_pos = i + cut_offset  # cut coordinate between bases cut_pos-1 and cut_pos
        starts.append(cut_pos)
        i = seq.find(site, i+1)
    return starts

def fragments_from_cuts(seq_len, cut_positions, circular=False):
    """Return list of fragment lengths in 5'->3' order for linear sequence.
       If circular=True, sequence is circular; cuts wrap around."""
    cuts = sorted(set(cut_positions))
    if not cuts:
        return [seq_len]
    if circular:
        # linearize starting at first cut, then compute gaps around circle
        ordered = cuts[:]
        frags = []
        for a,b in zip(ordered, ordered[1:]+ordered[:1]):
            if b> a:
                frags.append(b-a)
            else:
                frags.append(seq_len - a + b)
        return frags
    else:
        # linear: fragments are from start->first cut, between cuts, last cut->end
        frags = []
        prev = 0
        for c in cuts:
            frags.append(c - prev)
            prev = c
        frags.append(seq_len - prev)
        return frags

def report_single(seq, name, info):
    cuts = find_cuts(seq, info["site"], info["cut"])
    frags = fragments_from_cuts(len(seq), cuts, circular=False)
    print(f"== {name} ==")
    print(f"Recognition: {info['site']}, cut at offset {info['cut']}")
    print(f"Number of cleavages: {len(cuts)}")
    if cuts:
        print("Cleavage positions (1-based between-base coordinate):")
        print(", ".join(str(p) for p in cuts))
    else:
        print("Cleavage positions: none")
    print("Fragment lengths (bp):", ", ".join(str(x) for x in frags))
    print()

def combined_digest(seq, enames):
    # collect all cuts from listed enzymes and compute fragments
    cuts = []
    for e in enames:
        info = ENZYMES[e]
        cuts += find_cuts(seq, info["site"], info["cut"])
    frags = fragments_from_cuts(len(seq), cuts, circular=False)
    return sorted(frags, reverse=True)  # for gel display, show largest->smallest

def ascii_gel(band_sizes, max_width=40, gel_height=20, ladder=None):
    """
    Simple vertical gel simulation:
    - gel_height: number of rows (top is large fragments)
    - ladder: list of bp sizes (optional) to annotate left
    We'll map sizes log-scale to rows so large range looks reasonable.
    """
    if not band_sizes:
        band_sizes = []
    sizes = sorted(set(band_sizes), reverse=True)
    if not sizes:
        rows = ["|"+ " "*(max_width) + "|"] * gel_height
        return "\n".join(rows)
    mx = max(sizes)
    mn = min(sizes)
    import math
    def size_to_row(s):
        # map log(size) between 0..gel_height-1 (0 top)
        if s<=0: return gel_height-1
        if mx==mn:
            return gel_height//2
        return int((math.log(mx)-math.log(s))/(math.log(mx)-math.log(mn))*(gel_height-1))
    # build empty gel
    gel = [[" "]*max_width for _ in range(gel_height)]
    for s in sizes:
        r = size_to_row(s)
        # put a horizontal band: a run of '#' centered
        width = max(1, int(max_width * (0.5 * (s/mx) + 0.1)))  # heuristic width
        start = max(0, (max_width - width)//2)
        for c in range(start, min(max_width, start+width)):
            gel[r][c] = "#"
    lines = []
    for i,row in enumerate(gel):
        line = "|" + "".join(row) + "|"
        lines.append(line)
    return "\n".join(lines)

def main(path):
    seq = read_fasta(path)
    L = len(seq)
    print(f"Sequence length: {L} bp\n")
    # single enzyme reports
    for name, info in ENZYMES.items():
        report_single(seq, name, info)
    # Show combined digests: single enzymes, pairwise, all five
    combos = []
    # singles
    for e in ENZYMES:
        combos.append((e,))
    # pairs (first 4 pairs to keep output compact)
    names = list(ENZYMES.keys())
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            combos.append((names[i], names[j]))
    combos.append(tuple(names))  # full digest
    print("== Gel simulations (ASCII). Each gel: top = large fragments, bottom = small fragments ==\n")
    for combo in combos:
        label = "+".join(combo)
        frags = combined_digest(seq, combo)
        print(f"-- {label} ({len(frags)} fragments) --")
        print("Fragment sizes (bp):", ", ".join(str(x) for x in frags))
        print(ascii_gel(frags))
        print()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python digest_sim.py input.fasta")
        sys.exit(1)
    main(sys.argv[1])
