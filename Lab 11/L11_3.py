import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
import numpy as np
from collections import defaultdict, Counter
import math
import os
from math import log, exp


class SimilarityScorer:

    @staticmethod
    def equation1_jukes_cantor(align1, align2, match_score=2, mismatch_score=-1):
        L = len(align1)
        if L == 0:
            return 0

        matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')

        p_obs = matches / L

        if p_obs < 0.75:
            try:
                p_corrected = -0.75 * log(1.0 - (4.0/3.0) * p_obs)
            except (ValueError, ZeroDivisionError):
                p_corrected = p_obs
        else:
            p_corrected = 1.0

        similarity_score = min(100.0, p_corrected * 100)

        return similarity_score, {
            'raw_matches': matches,
            'observed_identity': p_obs * 100,
            'corrected_identity': p_corrected * 100,
            'jc_distance': -log(1 - (4/3) * p_obs) if p_obs < 0.75 else 0
        }

    @staticmethod
    def equation2_mutual_information(align1, align2, background_freq=None):
        nucleotides = ['A', 'C', 'G', 'T', '-']
        n = len(nucleotides)

        joint_counts = np.zeros((n, n))
        marginal1 = np.zeros(n)
        marginal2 = np.zeros(n)

        total_pairs = 0
        valid_pairs = 0

        for a, b in zip(align1, align2):
            if a in nucleotides and b in nucleotides:
                i = nucleotides.index(a)
                j = nucleotides.index(b)
                joint_counts[i, j] += 1

                if a != '-':
                    marginal1[i] += 1
                if b != '-':
                    marginal2[j] += 1

                total_pairs += 1
                if a != '-' and b != '-':
                    valid_pairs += 1

        if total_pairs == 0 or valid_pairs < 10:
            return 0, {'mutual_info': 0, 'alignment_length': len(align1)}

        joint_probs = joint_counts / total_pairs
        marginal1_probs = marginal1 / np.sum(marginal1) if np.sum(marginal1) > 0 else np.ones(n)/n
        marginal2_probs = marginal2 / np.sum(marginal2) if np.sum(marginal2) > 0 else np.ones(n)/n

        mi = 0.0
        for i in range(n):
            for j in range(n):
                if joint_probs[i, j] > 0 and marginal1_probs[i] > 0 and marginal2_probs[j] > 0:
                    mi += joint_probs[i, j] * log(joint_probs[i, j] /
                                                  (marginal1_probs[i] * marginal2_probs[j]))

        max_mi = log(n)
        normalized_mi = (mi / max_mi * 100) if max_mi > 0 else 0

        return min(100.0, normalized_mi), {
            'mutual_info': mi,
            'normalized_mi': normalized_mi,
            'max_possible_mi': max_mi,
            'alignment_pairs': valid_pairs
        }

    @staticmethod
    def equation3_weighted_motif(align1, align2, motif_weights=None,
                                 conservation_scores=None):
        L = len(align1)
        if L == 0:
            return 0, {}

        sub_matrix = {
            ('A', 'A'): 5, ('C', 'C'): 5, ('G', 'G'): 5, ('T', 'T'): 5,
            ('A', 'G'): -1, ('G', 'A'): -1,
            ('C', 'T'): -1, ('T', 'C'): -1,
            ('A', 'C'): -3, ('C', 'A'): -3,
            ('A', 'T'): -3, ('T', 'A'): -3,
            ('G', 'C'): -3, ('C', 'G'): -3,
            ('G', 'T'): -3, ('T', 'G'): -3,
        }

        if motif_weights is None:
            motif_weights = []
            center = L / 2
            sigma = L / 4
            for i in range(L):
                weight = exp(-0.5 * ((i - center) / sigma) ** 2)
                motif_weights.append(weight)
            motif_weights = [w * L / sum(motif_weights) for w in motif_weights]

        total_weight = 0
        weighted_score = 0

        for i, (a, b) in enumerate(zip(align1, align2)):
            if a == '-' and b == '-':
                continue

            weight = motif_weights[i] if i < len(motif_weights) else 1.0

            if (a, b) in sub_matrix:
                score = sub_matrix[(a, b)]
            elif a == '-' or b == '-':
                score = -4
            else:
                score = -3

            weighted_score += weight * score
            total_weight += weight

        if total_weight == 0:
            return 0, {}

        max_possible = 5 * total_weight
        min_possible = -4 * total_weight

        if max_possible > min_possible:
            normalized = ((weighted_score - min_possible) /
                          (max_possible - min_possible)) * 100
        else:
            normalized = 50

        return min(100.0, max(0.0, normalized)), {
            'raw_weighted_score': weighted_score,
            'total_weight': total_weight,
            'max_possible': max_possible,
            'min_possible': min_possible,
            'position_weights_used': True
        }

    @staticmethod
    def calculate_all_scores(align1, align2):
        scores = {}
        details = {}

        scores['jc_score'], details['jc_details'] = SimilarityScorer.equation1_jukes_cantor(
            align1, align2
        )

        scores['mi_score'], details['mi_details'] = SimilarityScorer.equation2_mutual_information(
            align1, align2
        )

        scores['wm_score'], details['wm_details'] = SimilarityScorer.equation3_weighted_motif(
            align1, align2
        )

        scores['combined_score'] = np.mean([scores['jc_score'],
                                            scores['mi_score'],
                                            scores['wm_score']])

        return scores, details


def hierarchical_alignment_with_scoring(seq1, seq2, chunk_size=2000):
    print(f"Aligning sequences: {len(seq1)}bp vs {len(seq2)}bp")

    sample_size = min(5000, len(seq1), len(seq2))
    sample1 = seq1[:sample_size]
    sample2 = seq2[:sample_size]

    seeds = find_kmer_matches(sample1, sample2, k=12)

    if not seeds:
        print("No seed matches found, using banded alignment...")
        align1, align2, score = banded_smith_waterman(
            sample1, sample2,
            bandwidth=100
        )

        if len(align1) > 0:
            scores, details = SimilarityScorer.calculate_all_scores(align1, align2)
            return align1, align2, scores, details
        else:
            return "", "", {}, {}

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

    all_alignments = []

    for cluster in clusters[:5]:
        pos1_min = min(s[0] for s in cluster)
        pos1_max = max(s[0] for s in cluster)
        pos2_min = min(s[1] for s in cluster)
        pos2_max = max(s[1] for s in cluster)

        start1 = max(0, pos1_min - 200)
        end1 = min(len(sample1), pos1_max + 800)
        start2 = max(0, pos2_min - 200)
        end2 = min(len(sample2), pos2_max + 800)

        region1 = sample1[start1:end1]
        region2 = sample2[start2:end2]

        align1, align2, score = banded_smith_waterman(
            region1, region2,
            bandwidth=100
        )

        if score > 30 and len(align1) > 50:
            all_alignments.append({
                'align1': align1,
                'align2': align2,
                'score': score,
                'pos1': start1,
                'pos2': start2
            })

    if all_alignments:
        all_alignments.sort(key=lambda x: x['score'], reverse=True)
        top_alignments = all_alignments[:min(3, len(all_alignments))]

        combined_align1 = "\n---\n".join([a['align1'] for a in top_alignments])
        combined_align2 = "\n---\n".join([a['align2'] for a in top_alignments])

        align1_for_scoring = combined_align1.replace('\n---\n', '')
        align2_for_scoring = combined_align2.replace('\n---\n', '')

        scores, details = SimilarityScorer.calculate_all_scores(
            align1_for_scoring, align2_for_scoring
        )

        return combined_align1, combined_align2, scores, details

    return "", "", {}, {}


class EnhancedGenomeAlignmentGUI:
    def __init__(self, root):
        self.root = root
        root.title("Genome Alignment with Multiple Similarity Metrics")
        root.geometry("1300x900")

        self.seq1 = ""
        self.seq2 = ""
        self.similarity_scores = {}
        self.score_details = {}

        self._setup_ui()

    def _setup_ui(self):
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.grid(row=0, column=0, sticky='nsew')

        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(3, weight=1)

        title = ttk.Label(main_frame, text="Genome Alignment: Multiple Similarity Metrics",
                          font=("Arial", 16, "bold"))
        title.grid(row=0, column=0, columnspan=2, pady=(0, 20))

        file_frame = ttk.LabelFrame(main_frame, text="Genome Input", padding=10)
        file_frame.grid(row=1, column=0, sticky='ew', pady=(0, 10), columnspan=2)
        file_frame.columnconfigure(1, weight=1)

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

        control_frame = ttk.Frame(main_frame)
        control_frame.grid(row=2, column=0, columnspan=2, pady=(0, 10), sticky='ew')

        self.align_btn = ttk.Button(control_frame, text="Run Alignment & Calculate Scores",
                                    command=self._run_alignment, state='disabled', width=30)
        self.align_btn.pack(side='left', padx=5)

        ttk.Button(control_frame, text="Load Demo Sequences",
                   command=self._load_demo_sequences).pack(side='left', padx=5)

        ttk.Button(control_frame, text="Clear Results",
                   command=self._clear_results).pack(side='right', padx=5)

        notebook_frame = ttk.Frame(main_frame)
        notebook_frame.grid(row=3, column=0, columnspan=2, sticky='nsew', pady=(10, 0))
        notebook_frame.columnconfigure(0, weight=1)
        notebook_frame.rowconfigure(0, weight=1)

        self.notebook = ttk.Notebook(notebook_frame)
        self.notebook.grid(row=0, column=0, sticky='nsew')

        self.alignment_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.alignment_tab, text="Alignment View")
        self._setup_alignment_tab()

        self.scoring_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.scoring_tab, text="Similarity Scores")
        self._setup_scoring_tab()

        self.analysis_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.analysis_tab, text="Detailed Analysis")
        self._setup_analysis_tab()

        self.visual_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.visual_tab, text="Visual Comparison")
        self._setup_visual_tab()

        self.status_var = tk.StringVar(value="Ready. Load genome files or use demo sequences.")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, relief='sunken')
        status_bar.grid(row=4, column=0, columnspan=2, sticky='ew', pady=(10, 0))

    def _setup_alignment_tab(self):
        align_frame = ttk.Frame(self.alignment_tab)
        align_frame.pack(fill='both', expand=True, padx=5, pady=5)

        text_frame = ttk.Frame(align_frame)
        text_frame.pack(fill='both', expand=True)

        left_frame = ttk.Frame(text_frame)
        left_frame.pack(side='left', fill='both', expand=True)
        ttk.Label(left_frame, text="Influenza Sequence", font=("Arial", 10, "bold")).pack()
        self.align_text1 = scrolledtext.ScrolledText(left_frame, height=20,
                                                     font=("Courier", 9), wrap='none')
        self.align_text1.pack(fill='both', expand=True)

        right_frame = ttk.Frame(text_frame)
        right_frame.pack(side='left', fill='both', expand=True)
        ttk.Label(right_frame, text="COVID-19 Sequence", font=("Arial", 10, "bold")).pack()
        self.align_text2 = scrolledtext.ScrolledText(right_frame, height=20,
                                                     font=("Courier", 9), wrap='none')
        self.align_text2.pack(fill='both', expand=True)

        viz_frame = ttk.Frame(align_frame)
        viz_frame.pack(fill='x', pady=(10, 0))
        ttk.Label(viz_frame, text="Match Visualization:", font=("Arial", 10, "bold")).pack(anchor='w')
        self.match_viz_text = scrolledtext.ScrolledText(viz_frame, height=3,
                                                        font=("Courier", 9), wrap='none')
        self.match_viz_text.pack(fill='x')

    def _setup_scoring_tab(self):
        scores_frame = ttk.Frame(self.scoring_tab)
        scores_frame.pack(fill='both', expand=True, padx=10, pady=10)

        panel_frame = ttk.Frame(scores_frame)
        panel_frame.pack(fill='both', expand=True)

        panel1 = ttk.LabelFrame(panel_frame, text="Equation 1: Jukes-Cantor Evolutionary Distance",
                                padding=10)
        panel1.grid(row=0, column=0, sticky='nsew', padx=5, pady=5)
        self.jc_score_var = tk.StringVar(value="Score: N/A")
        self.jc_details_var = tk.StringVar(value="No data")
        ttk.Label(panel1, textvariable=self.jc_score_var, font=("Arial", 12, "bold")).pack(anchor='w')
        ttk.Label(panel1, text="Details:", font=("Arial", 10, "bold")).pack(anchor='w', pady=(10, 0))
        ttk.Label(panel1, textvariable=self.jc_details_var, justify='left').pack(anchor='w', fill='x')

        panel2 = ttk.LabelFrame(panel_frame, text="Equation 2: Mutual Information Score",
                                padding=10)
        panel2.grid(row=0, column=1, sticky='nsew', padx=5, pady=5)
        self.mi_score_var = tk.StringVar(value="Score: N/A")
        self.mi_details_var = tk.StringVar(value="No data")
        ttk.Label(panel2, textvariable=self.mi_score_var, font=("Arial", 12, "bold")).pack(anchor='w')
        ttk.Label(panel2, text="Details:", font=("Arial", 10, "bold")).pack(anchor='w', pady=(10, 0))
        ttk.Label(panel2, textvariable=self.mi_details_var, justify='left').pack(anchor='w', fill='x')

        panel3 = ttk.LabelFrame(panel_frame, text="Equation 3: Weighted Position-Specific Score",
                                padding=10)
        panel3.grid(row=0, column=2, sticky='nsew', padx=5, pady=5)
        self.wm_score_var = tk.StringVar(value="Score: N/A")
        self.wm_details_var = tk.StringVar(value="No data")
        ttk.Label(panel3, textvariable=self.wm_score_var, font=("Arial", 12, "bold")).pack(anchor='w')
        ttk.Label(panel3, text="Details:", font=("Arial", 10, "bold")).pack(anchor='w', pady=(10, 0))
        ttk.Label(panel3, textvariable=self.wm_details_var, justify='left').pack(anchor='w', fill='x')

        combined_frame = ttk.LabelFrame(scores_frame, text="Combined Similarity Assessment",
                                        padding=10)
        combined_frame.pack(fill='x', pady=(10, 0))

        self.combined_score_var = tk.StringVar(value="Combined Score: N/A")
        self.assessment_var = tk.StringVar(value="Assessment: No data available")

        ttk.Label(combined_frame, textvariable=self.combined_score_var,
                  font=("Arial", 14, "bold")).pack(anchor='w')
        ttk.Label(combined_frame, textvariable=self.assessment_var,
                  font=("Arial", 11)).pack(anchor='w', pady=(5, 0))

        progress_frame = ttk.Frame(combined_frame)
        progress_frame.pack(fill='x', pady=(10, 0))

        ttk.Label(progress_frame, text="Score Comparison:").grid(row=0, column=0, sticky='w')

        self.jc_progress = ttk.Progressbar(progress_frame, length=200, mode='determinate')
        self.jc_progress.grid(row=0, column=1, padx=5)
        ttk.Label(progress_frame, text="JC").grid(row=0, column=2)

        self.mi_progress = ttk.Progressbar(progress_frame, length=200, mode='determinate')
        self.mi_progress.grid(row=0, column=3, padx=5)
        ttk.Label(progress_frame, text="MI").grid(row=0, column=4)

        self.wm_progress = ttk.Progressbar(progress_frame, length=200, mode='determinate')
        self.wm_progress.grid(row=0, column=5, padx=5)
        ttk.Label(progress_frame, text="WM").grid(row=0, column=6)

        panel_frame.columnconfigure(0, weight=1)
        panel_frame.columnconfigure(1, weight=1)
        panel_frame.columnconfigure(2, weight=1)

    def _setup_analysis_tab(self):
        analysis_frame = ttk.Frame(self.analysis_tab)
        analysis_frame.pack(fill='both', expand=True, padx=10, pady=10)

        ttk.Label(analysis_frame, text="Detailed Analysis Report",
                  font=("Arial", 12, "bold")).pack(anchor='w')

        self.analysis_text = scrolledtext.ScrolledText(analysis_frame, height=25,
                                                       font=("Courier", 10))
        self.analysis_text.pack(fill='both', expand=True, pady=(10, 0))

        ttk.Button(analysis_frame, text="Export Analysis Report",
                   command=self._export_analysis).pack(anchor='e', pady=(5, 0))

    def _setup_visual_tab(self):
        visual_frame = ttk.Frame(self.visual_tab)
        visual_frame.pack(fill='both', expand=True, padx=10, pady=10)

        self.visual_canvas = tk.Canvas(visual_frame, bg='white')
        self.visual_canvas.pack(fill='both', expand=True)

        legend_frame = ttk.Frame(visual_frame)
        legend_frame.pack(fill='x', pady=(10, 0))

        ttk.Label(legend_frame, text="Color Legend:").pack(side='left', padx=(0, 10))

        colors = [('#4CAF50', 'High Similarity (>70%)'),
                  ('#FFC107', 'Medium Similarity (40-70%)'),
                  ('#F44336', 'Low Similarity (<40%)'),
                  ('#2196F3', 'Gap/Insertion')]

        for color, label in colors:
            color_frame = ttk.Frame(legend_frame)
            color_frame.pack(side='left', padx=5)
            tk.Canvas(color_frame, width=20, height=20, bg=color, highlightthickness=0).pack(side='left')
            ttk.Label(color_frame, text=label).pack(side='left', padx=2)

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
            self.status_var.set("Files loaded. Ready to align.")

    def _load_demo_sequences(self):
        influenza_demo = """
        ATGGAAGATTTTGTGCGACAATGCTTCAATCCGATGATCGTTCAAAAGAA
        AGGACATGACAATTGAAAAAATCATGGCGATCAACGCAAAAGTCAGCAAA
        TGCTCAAGATGAAATGACGCGTAGACGAGTTGTTGATGAAGCTGCACGAA
        TTTGCAAGATGAATACACAGTTCAATGCGATCGCGATCGTGGCGATGCTG
        GAGATGGATGCATTGCTCCGGATTGAAGATATGAAAGACGATGTTGCTCA
        """

        covid_demo = """
        ATGTTCGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAA
        TCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACAC
        GTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCA
        ACTCAAGGACTTTTCTGGGTTGACACCTAAACCTGCTTGTCCTAGGTGGT
        ATTAATTGTTACGACTATTGTATACCTTACAATAGTGTAACTTCTTCAAT
        """

        self.seq1 = ''.join(influenza_demo.split())
        self.seq2 = ''.join(covid_demo.split())

        self.align_btn.config(state='normal')
        self.status_var.set("Demo sequences loaded. Ready to align.")

        messagebox.showinfo("Demo Loaded",
                            f"Loaded demo sequences:\n"
                            f"Influenza: {len(self.seq1)} bases\n"
                            f"COVID-19: {len(self.seq2)} bases")

    def _clear_results(self):
        self.align_text1.delete(1.0, tk.END)
        self.align_text2.delete(1.0, tk.END)
        self.match_viz_text.delete(1.0, tk.END)
        self.analysis_text.delete(1.0, tk.END)

        self.jc_score_var.set("Score: N/A")
        self.jc_details_var.set("No data")
        self.mi_score_var.set("Score: N/A")
        self.mi_details_var.set("No data")
        self.wm_score_var.set("Score: N/A")
        self.wm_details_var.set("No data")
        self.combined_score_var.set("Combined Score: N/A")
        self.assessment_var.set("Assessment: No data available")

        self.jc_progress['value'] = 0
        self.mi_progress['value'] = 0
        self.wm_progress['value'] = 0

        self.visual_canvas.delete("all")

        self.status_var.set("Results cleared. Ready for new analysis.")

    def _run_alignment(self):
        try:
            self.align_btn.config(state='disabled')
            self.status_var.set("Starting alignment and scoring...")
            self.root.update()

            if not hasattr(self, 'seq1') or not self.seq1:
                if not os.path.exists(self.influenza_path.get()):
                    messagebox.showerror("Error", "Influenza FASTA file not found!")
                    return
                if not os.path.exists(self.covid_path.get()):
                    messagebox.showerror("Error", "COVID-19 FASTA file not found!")
                    return

                self.seq1 = read_fasta(self.influenza_path.get())
                self.seq2 = read_fasta(self.covid_path.get())

            self.status_var.set("Performing hierarchical alignment...")
            self.root.update()

            align1, align2, scores, details = hierarchical_alignment_with_scoring(
                self.seq1, self.seq2
            )

            if not align1 or not align2:
                messagebox.showwarning("No Alignment",
                                       "No significant alignment found between sequences.")
                return

            self.similarity_scores = scores
            self.score_details = details

            self._display_alignment(align1, align2)
            self._display_scores()
            self._display_analysis()
            self._display_visual_comparison(align1, align2)

            self.status_var.set("Analysis complete! All similarity scores calculated.")

        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")
            self.status_var.set("Analysis failed")
        finally:
            self.align_btn.config(state='normal')

    def _display_alignment(self, align1, align2):
        self.align_text1.delete(1.0, tk.END)
        self.align_text2.delete(1.0, tk.END)
        self.match_viz_text.delete(1.0, tk.END)

        chunk_size = 60
        viz_lines = []

        for i in range(0, min(len(align1), len(align2)), chunk_size):
            chunk1 = align1[i:i + chunk_size]
            chunk2 = align2[i:i + chunk_size]

            self.align_text1.insert(tk.END, chunk1 + "\n")
            self.align_text2.insert(tk.END, chunk2 + "\n")

            viz_line = []
            for a, b in zip(chunk1, chunk2):
                if a == b and a != '-':
                    viz_line.append('|')
                elif a == '-' or b == '-':
                    viz_line.append(' ')
                elif (a in 'AG' and b in 'AG') or (a in 'CT' and b in 'CT'):
                    viz_line.append(':')
                else:
                    viz_line.append('.')
            viz_lines.append(''.join(viz_line))

        for viz_line in viz_lines:
            self.match_viz_text.insert(tk.END, viz_line + "\n")

    def _display_scores(self):
        jc_score = self.similarity_scores.get('jc_score', 0)
        jc_details = self.score_details.get('jc_details', {})
        self.jc_score_var.set(f"Jukes-Cantor Score: {jc_score:.1f}%")

        details_text = ""
        if jc_details:
            details_text = f"Raw Matches: {jc_details.get('raw_matches', 'N/A')}\n"
            details_text += f"Observed Identity: {jc_details.get('observed_identity', 0):.1f}%\n"
            details_text += f"Corrected Identity: {jc_details.get('corrected_identity', 0):.1f}%\n"
            details_text += f"JC Distance: {jc_details.get('jc_distance', 0):.3f}"
        self.jc_details_var.set(details_text)

        mi_score = self.similarity_scores.get('mi_score', 0)
        mi_details = self.score_details.get('mi_details', {})
        self.mi_score_var.set(f"Mutual Information: {mi_score:.1f}%")

        details_text = ""
        if mi_details:
            details_text = f"Raw MI: {mi_details.get('mutual_info', 0):.3f}\n"
            details_text += f"Normalized MI: {mi_details.get('normalized_mi', 0):.1f}%\n"
            details_text += f"Max Possible MI: {mi_details.get('max_possible_mi', 0):.3f}\n"
            details_text += f"Alignment Pairs: {mi_details.get('alignment_pairs', 0)}"
        self.mi_details_var.set(details_text)

        wm_score = self.similarity_scores.get('wm_score', 0)
        wm_details = self.score_details.get('wm_details', {})
        self.wm_score_var.set(f"Weighted Motif Score: {wm_score:.1f}%")

        details_text = ""
        if wm_details:
            details_text = f"Weighted Score: {wm_details.get('raw_weighted_score', 0):.1f}\n"
            details_text += f"Total Weight: {wm_details.get('total_weight', 0):.1f}\n"
            details_text += f"Position Weights Used: {'Yes' if wm_details.get('position_weights_used', False) else 'No'}"
        self.wm_details_var.set(details_text)

        combined = self.similarity_scores.get('combined_score', 0)
        self.combined_score_var.set(f"Combined Similarity Score: {combined:.1f}%")

        if combined >= 70:
            assessment = "HIGH SIMILARITY: Sequences show significant evolutionary relationship"
        elif combined >= 40:
            assessment = "MODERATE SIMILARITY: Sequences share some conserved regions"
        else:
            assessment = "LOW SIMILARITY: Limited sequence conservation detected"
        self.assessment_var.set(f"Assessment: {assessment}")

        self.jc_progress['value'] = jc_score
        self.mi_progress['value'] = mi_score
        self.wm_progress['value'] = wm_score

    def _display_analysis(self):
        self.analysis_text.delete(1.0, tk.END)

        report = "GENOME SIMILARITY ANALYSIS REPORT\n"
        report += "=" * 50 + "\n\n"

        report += "SCORING EQUATIONS USED:\n"
        report += "1. Jukes-Cantor Evolutionary Distance: Corrects for multiple substitutions\n"
        report += "2. Mutual Information: Measures statistical dependence between positions\n"
        report += "3. Weighted Position-Specific: Considers biological importance of regions\n\n"

        report += "RESULTS SUMMARY:\n"
        report += "-" * 30 + "\n"

        for eq_name, score in [("Jukes-Cantor", self.similarity_scores.get('jc_score', 0)),
                               ("Mutual Information", self.similarity_scores.get('mi_score', 0)),
                               ("Weighted Motif", self.similarity_scores.get('wm_score', 0))]:
            report += f"{eq_name}: {score:.1f}%\n"

        report += f"\nCombined Score: {self.similarity_scores.get('combined_score', 0):.1f}%\n\n"

        report += "INTERPRETATION:\n"
        report += "-" * 30 + "\n"

        combined = self.similarity_scores.get('combined_score', 0)
        if combined >= 80:
            report += "◆ Very high similarity suggesting close evolutionary relationship\n"
            report += "◆ Possible homologous regions or conserved functional domains\n"
        elif combined >= 60:
            report += "◆ Significant similarity with detectable evolutionary relationship\n"
            report += "◆ Some conserved regions likely with functional importance\n"
        elif combined >= 40:
            report += "◆ Moderate similarity indicating distant evolutionary relationship\n"
            report += "◆ Limited conserved regions, possibly convergent evolution\n"
        else:
            report += "◆ Low similarity suggesting independent evolutionary origins\n"
            report += "◆ Any matches likely due to chance or very short motifs\n"

        report += "\nSCORE AGREEMENT ANALYSIS:\n"
        report += "-" * 30 + "\n"

        scores = [self.similarity_scores.get('jc_score', 0),
                  self.similarity_scores.get('mi_score', 0),
                  self.similarity_scores.get('wm_score', 0)]

        score_range = max(scores) - min(scores)
        if score_range < 20:
            report += "◆ Scores show good agreement across different metrics\n"
            report += "◆ Result is robust to different similarity definitions\n"
        elif score_range < 40:
            report += "◆ Moderate agreement between different scoring methods\n"
            report += "◆ Similarity pattern depends on measurement approach\n"
        else:
            report += "◆ Poor agreement between different scoring methods\n"
            report += "◆ Similarity assessment highly metric-dependent\n"

        self.analysis_text.insert(1.0, report)

    def _display_visual_comparison(self, align1, align2):
        self.visual_canvas.delete("all")

        canvas_width = self.visual_canvas.winfo_width() - 40
        canvas_height = self.visual_canvas.winfo_height() - 40

        if canvas_width < 100 or canvas_height < 100:
            canvas_width, canvas_height = 600, 400

        L = len(align1)
        if L == 0:
            return

        cell_width = min(5, canvas_width / L)
        cell_height = 20

        y_position = 30

        for seq_idx, seq in enumerate([align1, align2]):
            y = y_position + seq_idx * (cell_height + 5)

            self.visual_canvas.create_text(10, y + cell_height/2,
                                           text="Influenza" if seq_idx == 0 else "COVID-19",
                                           anchor='w', font=("Arial", 9, "bold"))

            x = 100
            for i, base in enumerate(seq):
                if i < len(align2):
                    other_base = align2[i] if seq_idx == 0 else align1[i]

                    if base == '-' or other_base == '-':
                        color = '#2196F3'
                    elif base == other_base:
                        match_score = self.similarity_scores.get('combined_score', 0)
                        if match_score > 70:
                            color = '#4CAF50'
                        elif match_score > 40:
                            color = '#FFC107'
                        else:
                            color = '#F44336'
                    elif (base in 'AG' and other_base in 'AG') or (base in 'CT' and other_base in 'CT'):
                        color = '#9C27B0'
                    else:
                        color = '#FF5722'
                else:
                    color = '#9E9E9E'

                self.visual_canvas.create_rectangle(
                    x, y, x + cell_width, y + cell_height,
                    fill=color, outline=color
                )

                if cell_width > 8:
                    self.visual_canvas.create_text(
                        x + cell_width/2, y + cell_height/2,
                        text=base, font=("Courier", 8)
                    )

                x += cell_width

        y_position += 2 * (cell_height + 5) + 20

        scores = [
            ("JC", self.similarity_scores.get('jc_score', 0)),
            ("MI", self.similarity_scores.get('mi_score', 0)),
            ("WM", self.similarity_scores.get('wm_score', 0))
        ]

        bar_width = 100
        bar_spacing = 150

        for i, (label, score) in enumerate(scores):
            x_start = 50 + i * bar_spacing

            self.visual_canvas.create_rectangle(
                x_start, y_position, x_start + bar_width, y_position + 20,
                fill='#E0E0E0', outline='#BDBDBD'
            )

            fill_width = (score / 100) * bar_width
            color = '#4CAF50' if score > 70 else '#FFC107' if score > 40 else '#F44336'

            self.visual_canvas.create_rectangle(
                x_start, y_position, x_start + fill_width, y_position + 20,
                fill=color, outline=color
            )

            self.visual_canvas.create_text(
                x_start + bar_width/2, y_position - 10,
                text=label, font=("Arial", 10, "bold")
            )

            self.visual_canvas.create_text(
                x_start + bar_width/2, y_position + 10,
                text=f"{score:.1f}%", font=("Arial", 9)
            )

        self.visual_canvas.create_text(
            canvas_width/2, 15,
            text="Sequence Similarity Visualization",
            font=("Arial", 12, "bold")
        )

    def _export_analysis(self):
        filename = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                with open(filename, 'w') as f:
                    f.write(self.analysis_text.get(1.0, tk.END))
                messagebox.showinfo("Export Successful",
                                    f"Analysis report saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Export Failed", f"Could not save file: {str(e)}")


def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip().upper()
    return sequence

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

def banded_smith_waterman(seq1, seq2, match=2, mismatch=-1, gap_open=-2,
                          gap_extend=-1, bandwidth=50):
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

def main():
    root = tk.Tk()
    app = EnhancedGenomeAlignmentGUI(root)

    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)

    root.mainloop()

if __name__ == "__main__":
    main()