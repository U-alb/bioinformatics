import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
import numpy as np
from collections import defaultdict
import math
import os


def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip().upper()
    return sequence

def chunk_sequence(seq, chunk_size=1000, overlap=200):
    chunks = []
    positions = []
    for i in range(0, len(seq), chunk_size - overlap):
        chunk = seq[i:i + chunk_size]
        if len(chunk) >= 500:  # Minimum chunk size
            chunks.append(chunk)
            positions.append((i, i + len(chunk)))
    return chunks, positions

def find_kmer_matches(seq1, seq2, k=12):
    kmer_dict = defaultdict(list)

    for i in range(len(seq2) - k + 1):
        kmer = seq2[i:i + k]
        kmer_dict[kmer].append(i)

    matches = []
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i:i + k]
        if kmer in kmer_dict:
            for j in kmer_dict[kmer]:
                matches.append((i, j))

    return matches

def banded_smith_waterman(seq1, seq2, match=2, mismatch=-1, gap_open=-2, gap_extend=-1, bandwidth=50):
    m, n = len(seq1), len(seq2)

    H = np.zeros((m + 1, n + 1), dtype=np.int32)
    E = np.zeros((m + 1, n + 1), dtype=np.int32)
    F = np.zeros((m + 1, n + 1), dtype=np.int32)

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        band_start = max(1, i - bandwidth)
        band_end = min(n, i + bandwidth)

        for j in range(band_start, band_end + 1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch

            diag = H[i-1][j-1] + score
            h_up = H[i-1][j] + gap_open if E[i-1][j] == 0 else H[i-1][j] + gap_extend
            h_left = H[i][j-1] + gap_open if F[i][j-1] == 0 else H[i][j-1] + gap_extend

            H[i][j] = max(0, diag, h_up, h_left)

            E[i][j] = max(0, h_up)
            F[i][j] = max(0, h_left)

            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j)

    if max_score == 0:
        return "", "", 0

    i, j = max_pos
    align1, align2 = [], []

    while i > 0 and j > 0 and H[i][j] > 0:
        if H[i][j] == H[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif H[i][j] == (H[i-1][j] + gap_open if E[i-1][j] == 0 else H[i-1][j] + gap_extend):
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1

    align1 = ''.join(reversed(align1))
    align2 = ''.join(reversed(align2))

    return align1, align2, max_score

def hierarchical_alignment(seq1, seq2, chunk_size=2000):
    print(f"Aligning sequences: {len(seq1)}bp vs {len(seq2)}bp")

    print("Step 1: Finding seed matches...")
    seeds = find_kmer_matches(seq1[:10000], seq2[:10000], k=15)

    if not seeds:
        print("No seed matches found, using banded alignment...")
        align1, align2, score = banded_smith_waterman(
            seq1[:10000], seq2[:10000],
            bandwidth=100
        )
        return align1, align2, score

    print(f"Step 2: Grouping {len(seeds)} seeds...")
    clusters = []
    if seeds:
        seeds.sort(key=lambda x: x[0])
        current_cluster = [seeds[0]]

        for seed in seeds[1:]:
            last_seed = current_cluster[-1]
            if (seed[0] - last_seed[0] < 500 and
                    abs(seed[1] - last_seed[1]) < 500):
                current_cluster.append(seed)
            else:
                if len(current_cluster) >= 3:
                    clusters.append(current_cluster)
                current_cluster = [seed]

        if len(current_cluster) >= 3:
            clusters.append(current_cluster)

    print(f"Step 3: Aligning {len(clusters)} regions...")
    all_alignments = []

    for cluster in clusters[:10]:
        pos1_min = min(s[0] for s in cluster)
        pos1_max = max(s[0] for s in cluster)
        pos2_min = min(s[1] for s in cluster)
        pos2_max = max(s[1] for s in cluster)

        start1 = max(0, pos1_min - 500)
        end1 = min(len(seq1), pos1_max + 1500)
        start2 = max(0, pos2_min - 500)
        end2 = min(len(seq2), pos2_max + 1500)

        region1 = seq1[start1:end1]
        region2 = seq2[start2:end2]

        align1, align2, score = banded_smith_waterman(
            region1, region2,
            bandwidth=100
        )

        if score > 50:
            all_alignments.append({
                'align1': align1,
                'align2': align2,
                'score': score,
                'pos1': start1,
                'pos2': start2
            })

    if all_alignments:
        all_alignments.sort(key=lambda x: x['score'], reverse=True)

        top_alignments = all_alignments[:3]

        combined_align1 = "\n---\n".join([a['align1'] for a in top_alignments])
        combined_align2 = "\n---\n".join([a['align2'] for a in top_alignments])
        total_score = sum(a['score'] for a in top_alignments)

        return combined_align1, combined_align2, total_score

    return "", "", 0

def calculate_similarity_metrics(align1, align2):
    if not align1 or not align2:
        return {}

    seq1_no_gaps = align1.replace('-', '')
    seq2_no_gaps = align2.replace('-', '')

    matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    mismatches = sum(1 for a, b in zip(align1, align2) if a != b and a != '-' and b != '-')
    gaps = align1.count('-') + align2.count('-')

    total_positions = len(align1)

    identity = (matches / total_positions * 100) if total_positions > 0 else 0
    similarity = ((matches + 0.5 * mismatches) / total_positions * 100) if total_positions > 0 else 0
    gap_percentage = (gaps / (2 * total_positions) * 100) if total_positions > 0 else 0

    gc1 = (seq1_no_gaps.count('G') + seq1_no_gaps.count('C')) / len(seq1_no_gaps) * 100 if seq1_no_gaps else 0
    gc2 = (seq2_no_gaps.count('G') + seq2_no_gaps.count('C')) / len(seq2_no_gaps) * 100 if seq2_no_gaps else 0

    return {
        'identity': identity,
        'similarity': similarity,
        'gap_percentage': gap_percentage,
        'matches': matches,
        'mismatches': mismatches,
        'gaps': gaps,
        'total_positions': total_positions,
        'gc_content_seq1': gc1,
        'gc_content_seq2': gc2,
        'gc_difference': abs(gc1 - gc2)
    }

def create_dot_plot_matrix(seq1, seq2, window_size=20, threshold=0.7):
    max_size = 500
    if len(seq1) > max_size:
        seq1 = seq1[:max_size]
    if len(seq2) > max_size:
        seq2 = seq2[:max_size]

    m, n = len(seq1), len(seq2)
    matrix = np.zeros((m, n), dtype=np.float32)

    for i in range(m - window_size + 1):
        window1 = seq1[i:i + window_size]
        for j in range(n - window_size + 1):
            window2 = seq2[j:j + window_size]
            matches = sum(1 for a, b in zip(window1, window2) if a == b)
            similarity = matches / window_size
            if similarity >= threshold:
                for di in range(window_size):
                    for dj in range(window_size):
                        if i + di < m and j + dj < n:
                            matrix[i + di][j + dj] = max(matrix[i + di][j + dj], similarity)

    return matrix, seq1, seq2


class GenomeAlignmentGUI:
    def __init__(self, root):
        self.root = root
        root.title("Influenza vs COVID-19 Genome Alignment")
        root.geometry("1200x800")

        self.seq1 = ""
        self.seq2 = ""
        self.seq1_name = "Influenza"
        self.seq2_name = "COVID-19"

        self._setup_ui()

    def _setup_ui(self):
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.grid(row=0, column=0, sticky='nsew')

        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(3, weight=1)

        title = ttk.Label(main_frame, text="Genome Alignment: Influenza vs COVID-19",
                          font=("Arial", 16, "bold"))
        title.grid(row=0, column=0, columnspan=3, pady=(0, 20))

        file_frame = ttk.LabelFrame(main_frame, text="Genome Input", padding=10)
        file_frame.grid(row=1, column=0, columnspan=3, sticky='ew', pady=(0, 10))
        file_frame.columnconfigure(1, weight=1)
        file_frame.columnconfigure(3, weight=1)

        ttk.Label(file_frame, text="Influenza Genome (FASTA):").grid(row=0, column=0, sticky='w')
        self.influenza_path = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.influenza_path, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(file_frame, text="Browse",
                   command=lambda: self._browse_file(self.influenza_path, "influenza")).grid(row=0, column=2)

        ttk.Label(file_frame, text="COVID-19 Genome (FASTA):").grid(row=1, column=0, sticky='w', pady=(10, 0))
        self.covid_path = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.covid_path, width=50).grid(row=1, column=1, padx=5, pady=(10, 0))
        ttk.Button(file_frame, text="Browse",
                   command=lambda: self._browse_file(self.covid_path, "covid")).grid(row=1, column=2, pady=(10, 0))

        ttk.Button(file_frame, text="Load Demo Sequences",
                   command=self._load_demo_sequences).grid(row=2, column=1, pady=(10, 0))

        control_frame = ttk.Frame(main_frame)
        control_frame.grid(row=2, column=0, columnspan=3, pady=(0, 10), sticky='ew')

        ttk.Label(control_frame, text="Alignment Parameters:").pack(side='left', padx=(0, 10))

        ttk.Label(control_frame, text="Match:").pack(side='left')
        self.match_var = tk.IntVar(value=2)
        ttk.Spinbox(control_frame, from_=1, to=10, textvariable=self.match_var, width=5).pack(side='left', padx=(0, 10))

        ttk.Label(control_frame, text="Mismatch:").pack(side='left')
        self.mismatch_var = tk.IntVar(value=-1)
        ttk.Spinbox(control_frame, from_=-5, to=0, textvariable=self.mismatch_var, width=5).pack(side='left', padx=(0, 10))

        ttk.Label(control_frame, text="Gap Open:").pack(side='left')
        self.gap_open_var = tk.IntVar(value=-2)
        ttk.Spinbox(control_frame, from_=-5, to=0, textvariable=self.gap_open_var, width=5).pack(side='left', padx=(0, 10))

        self.align_btn = ttk.Button(control_frame, text="Run Genome Alignment",
                                    command=self._run_genome_alignment, state='disabled')
        self.align_btn.pack(side='right')

        viz_frame = ttk.Frame(main_frame)
        viz_frame.grid(row=3, column=0, columnspan=3, sticky='nsew', pady=(10, 0))
        viz_frame.columnconfigure(0, weight=1)
        viz_frame.columnconfigure(1, weight=1)
        viz_frame.rowconfigure(0, weight=1)

        self.notebook = ttk.Notebook(viz_frame)
        self.notebook.grid(row=0, column=0, columnspan=2, sticky='nsew')

        self.alignment_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.alignment_tab, text="Alignment")

        align_frame = ttk.Frame(self.alignment_tab)
        align_frame.pack(fill='both', expand=True, padx=5, pady=5)

        self.align_text1 = scrolledtext.ScrolledText(align_frame, height=15, width=60,
                                                     font=("Courier", 10), wrap='none')
        self.align_text1.pack(side='left', fill='both', expand=True)

        self.align_text2 = scrolledtext.ScrolledText(align_frame, height=15, width=60,
                                                     font=("Courier", 10), wrap='none')
        self.align_text2.pack(side='left', fill='both', expand=True)

        self.stats_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.stats_tab, text="Statistics")

        self.stats_text = scrolledtext.ScrolledText(self.stats_tab, height=20, width=80,
                                                    font=("Courier", 10))
        self.stats_text.pack(fill='both', expand=True, padx=5, pady=5)

        self.dotplot_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.dotplot_tab, text="Dot Plot")

        self.dotplot_canvas = tk.Canvas(self.dotplot_tab, bg='white')
        self.dotplot_canvas.pack(fill='both', expand=True, padx=5, pady=5)

        self.status_var = tk.StringVar(value="Ready. Please load genome files.")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief='sunken')
        status_bar.grid(row=4, column=0, columnspan=3, sticky='ew', pady=(10, 0))

    def _browse_file(self, path_var, seq_type):
        filename = filedialog.askopenfilename(
            title=f"Select {seq_type.capitalize()} FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
        )
        if filename:
            path_var.set(filename)
            self._check_files_loaded()

    def _check_files_loaded(self):
        if self.influenza_path.get() and self.covid_path.get():
            self.align_btn.config(state='normal')
            self.status_var.set("Files loaded. Click 'Run Genome Alignment' to start.")

    def _load_demo_sequences(self):
        influenza_sample = """
        ATGGAAGATTTTGTGCGACAATGCTTCAATCCGATGATCGTTCAAAAGAA
        AGGACATGACAATTGAAAAAATCATGGCGATCAACGCAAAAGTCAGCAAA
        TGCTCAAGATGAAATGACGCGTAGACGAGTTGTTGATGAAGCTGCACGAA
        TTTGCAAGATGAATACACAGTTCAATGCGATCGCGATCGTGGCGATGCTG
        GAGATGGATGCATTGCTCCGGATTGAAGATATGAAAGACGATGTTGCTCA
        """

        covid_sample = """
        ATGTTCGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAA
        TCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACAC
        GTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCA
        ACTCAAGGACTTTTCTGGGTTGACACCTAAACCTGCTTGTCCTAGGTGGT
        ATTAATTGTTACGACTATTGTATACCTTACAATAGTGTAACTTCTTCAAT
        """

        self.seq1 = ''.join(influenza_sample.split())
        self.seq2 = ''.join(covid_sample.split())

        self.align_btn.config(state='normal')
        self.status_var.set("Demo sequences loaded. Click 'Run Genome Alignment' to start.")

        messagebox.showinfo("Demo Loaded",
                            f"Loaded demo sequences:\n"
                            f"Influenza: {len(self.seq1)} bases\n"
                            f"COVID-19: {len(self.seq2)} bases")

    def _run_genome_alignment(self):
        try:
            self.align_btn.config(state='disabled')
            self.status_var.set("Loading sequences...")
            self.root.update()

            if hasattr(self, 'seq1') and self.seq1:
                pass
            else:
                if not os.path.exists(self.influenza_path.get()):
                    messagebox.showerror("Error", "Influenza FASTA file not found!")
                    return
                if not os.path.exists(self.covid_path.get()):
                    messagebox.showerror("Error", "COVID-19 FASTA file not found!")
                    return

                self.seq1 = read_fasta(self.influenza_path.get())
                self.seq2 = read_fasta(self.covid_path.get())

            if len(self.seq1) < 100 or len(self.seq2) < 100:
                messagebox.showerror("Error", "Sequences are too short for meaningful alignment!")
                return

            self.status_var.set(f"Sequences loaded: Influenza={len(self.seq1)}bp, COVID-19={len(self.seq2)}bp")
            self.root.update()

            display_seq1 = self.seq1[:10000] if len(self.seq1) > 10000 else self.seq1
            display_seq2 = self.seq2[:10000] if len(self.seq2) > 10000 else self.seq2

            self.status_var.set("Running hierarchical alignment...")
            self.root.update()

            align1, align2, score = hierarchical_alignment(
                display_seq1,
                display_seq2,
                chunk_size=2000
            )

            self._display_alignment(align1, align2, score)

            metrics = calculate_similarity_metrics(align1, align2)
            self._display_statistics(metrics)

            self._display_dot_plot(display_seq1[:500], display_seq2[:500])

            self.status_var.set(f"Alignment complete! Score: {score}")

        except Exception as e:
            messagebox.showerror("Error", f"Alignment failed: {str(e)}")
            self.status_var.set("Alignment failed")
        finally:
            self.align_btn.config(state='normal')

    def _display_alignment(self, align1, align2, score):
        self.align_text1.delete(1.0, tk.END)
        self.align_text2.delete(1.0, tk.END)

        if not align1 or not align2:
            self.align_text1.insert(1.0, "No significant alignment found.")
            return

        chunk_size = 60
        for i in range(0, min(len(align1), len(align2)), chunk_size):
            chunk1 = align1[i:i + chunk_size]
            chunk2 = align2[i:i + chunk_size]

            self.align_text1.insert(tk.END, chunk1 + "\n")
            self.align_text2.insert(tk.END, chunk2 + "\n")

        header = f"Alignment Score: {score}\n"
        header += "-" * 60 + "\n"
        self.align_text1.insert(1.0, header)
        self.align_text2.insert(1.0, header)

    def _display_statistics(self, metrics):
        self.stats_text.delete(1.0, tk.END)

        if not metrics:
            self.stats_text.insert(1.0, "No statistics available.")
            return

        stats_text = "ALIGNMENT STATISTICS\n"
        stats_text += "=" * 40 + "\n\n"

        stats_text += f"Sequence Identity: {metrics['identity']:.2f}%\n"
        stats_text += f"Sequence Similarity: {metrics['similarity']:.2f}%\n"
        stats_text += f"Gap Percentage: {metrics['gap_percentage']:.2f}%\n\n"

        stats_text += f"Matches: {metrics['matches']}\n"
        stats_text += f"Mismatches: {metrics['mismatches']}\n"
        stats_text += f"Gaps: {metrics['gaps']}\n"
        stats_text += f"Total Positions: {metrics['total_positions']}\n\n"

        stats_text += "GC CONTENT\n"
        stats_text += "-" * 20 + "\n"
        stats_text += f"Influenza: {metrics['gc_content_seq1']:.2f}%\n"
        stats_text += f"COVID-19: {metrics['gc_content_seq2']:.2f}%\n"
        stats_text += f"Difference: {metrics['gc_difference']:.2f}%\n\n"

        stats_text += "INTERPRETATION\n"
        stats_text += "-" * 20 + "\n"
        if metrics['identity'] > 70:
            stats_text += "High sequence identity detected\n"
        elif metrics['identity'] > 40:
            stats_text += "Moderate sequence identity detected\n"
        else:
            stats_text += "Low sequence identity detected\n"

        self.stats_text.insert(1.0, stats_text)

    def _display_dot_plot(self, seq1, seq2):
        self.dotplot_canvas.delete("all")

        matrix, seq1_display, seq2_display = create_dot_plot_matrix(
            seq1, seq2, window_size=15, threshold=0.6
        )

        canvas_width = self.dotplot_canvas.winfo_width() - 100
        canvas_height = self.dotplot_canvas.winfo_height() - 100

        if canvas_width < 100 or canvas_height < 100:
            canvas_width, canvas_height = 400, 400

        m, n = matrix.shape
        cell_width = canvas_width / n
        cell_height = canvas_height / m

        for i in range(m):
            for j in range(n):
                if matrix[i][j] > 0:
                    intensity = int(255 * matrix[i][j])
                    color = f'#{intensity:02x}{intensity:02x}{255:02x}'

                    x1 = 50 + j * cell_width
                    y1 = 50 + i * cell_height
                    x2 = x1 + cell_width
                    y2 = y1 + cell_height

                    self.dotplot_canvas.create_rectangle(
                        x1, y1, x2, y2,
                        fill=color, outline=color
                    )

        self.dotplot_canvas.create_line(50, 50, 50, 50 + canvas_height, fill='black', width=2)
        self.dotplot_canvas.create_line(50, 50 + canvas_height, 50 + canvas_width, 50 + canvas_height,
                                        fill='black', width=2)

        self.dotplot_canvas.create_text(25, 50 + canvas_height/2, text=self.seq1_name,
                                        angle=90, font=("Arial", 10, "bold"))
        self.dotplot_canvas.create_text(50 + canvas_width/2, 50 + canvas_height + 30,
                                        text=self.seq2_name, font=("Arial", 10, "bold"))

        self.dotplot_canvas.create_text(50 + canvas_width/2, 20,
                                        text="Sequence Similarity Dot Plot",
                                        font=("Arial", 12, "bold"))

        legend_x = 50 + canvas_width + 10
        legend_y = 100

        self.dotplot_canvas.create_text(legend_x, legend_y - 20, text="Similarity",
                                        font=("Arial", 9, "bold"))

        for i in range(5):
            y = legend_y + i * 20
            intensity = int(255 * (i / 4))
            color = f'#{intensity:02x}{intensity:02x}{255:02x}'

            self.dotplot_canvas.create_rectangle(legend_x, y, legend_x + 20, y + 15,
                                                 fill=color, outline='black')
            self.dotplot_canvas.create_text(legend_x + 40, y + 7,
                                            text=f"{i/4:.1f}", anchor='w')

def main():
    root = tk.Tk()
    app = GenomeAlignmentGUI(root)

    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)

    root.mainloop()

if __name__ == "__main__":
    main()