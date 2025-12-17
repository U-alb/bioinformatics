import math
import re
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


class GenomeMotifScanner:
    def __init__(self, motif_sequences=None):
        self.nucleotides = ['A', 'C', 'G', 'T']
        self.background_freq = 0.25
        self.motif_length = 9
        self.motif_sequences = motif_sequences or []

        if motif_sequences:
            self._build_motif_model()

    def _build_motif_model(self):
        print("Building motif model...")
        self.count_matrix = self.create_count_matrix(self.motif_sequences)
        self.weight_matrix = self.create_weight_matrix(self.count_matrix)
        self.rel_freq_matrix = self.create_relative_frequencies(self.weight_matrix)
        self.log_matrix = self.create_log_likelihood_matrix(self.rel_freq_matrix)
        print("Motif model built successfully!\n")

    def create_count_matrix(self, sequences):
        count_matrix = {nt: [0] * self.motif_length for nt in self.nucleotides}

        for i in range(self.motif_length):
            for seq in sequences:
                if i < len(seq):
                    nt = seq[i]
                    if nt in count_matrix:
                        count_matrix[nt][i] += 1

        return count_matrix

    def create_weight_matrix(self, count_matrix, pseudo_count=1):
        weight_matrix = {nt: [0] * self.motif_length for nt in self.nucleotides}
        num_sequences = sum(count_matrix[self.nucleotides[0]][:])

        for nt in self.nucleotides:
            for i in range(self.motif_length):
                weight_matrix[nt][i] = (count_matrix[nt][i] + pseudo_count) / (num_sequences + len(self.nucleotides) * pseudo_count)

        return weight_matrix

    def create_relative_frequencies(self, weight_matrix):
        return weight_matrix

    def create_log_likelihood_matrix(self, relative_freq_matrix):
        log_matrix = {nt: [0] * self.motif_length for nt in self.nucleotides}

        for nt in self.nucleotides:
            for i in range(self.motif_length):
                freq = relative_freq_matrix[nt][i]
                log_matrix[nt][i] = math.log(freq / self.background_freq)

        return log_matrix

    def read_fasta_file(self, filepath: str) -> Dict[str, str]:
        genomes = {}
        current_id = None
        current_seq = []

        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        genomes[current_id] = ''.join(current_seq)

                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.upper())

        if current_id:
            genomes[current_id] = ''.join(current_seq)

        print(f"Read {len(genomes)} genomes from FASTA file")
        return genomes

    def calculate_window_score(self, window: str, log_matrix: Dict) -> float:
        if len(window) != self.motif_length:
            return float('-inf')

        score = 0
        for i, nt in enumerate(window):
            if nt in log_matrix:
                score += log_matrix[nt][i]
            else:
                score += math.log(self.background_freq / self.background_freq)

        return score

    def scan_genome(self, genome_id: str, sequence: str, log_matrix: Dict) -> Tuple[List[int], List[float]]:
        positions = []
        scores = []

        for i in range(len(sequence) - self.motif_length + 1):
            window = sequence[i:i + self.motif_length]
            score = self.calculate_window_score(window, log_matrix)
            positions.append(i)
            scores.append(score)

        return positions, scores

    def find_top_motif_locations(self, positions: List[int], scores: List[float], top_n: int = 10) -> List[Tuple[int, float, str]]:
        if not scores:
            return []

        score_positions = list(zip(positions, scores))
        score_positions.sort(key=lambda x: x[1], reverse=True)

        return score_positions[:top_n]

    def plot_genome_signal(self, genome_id: str, positions: List[int], scores: List[float],
                           sequence: str = None, save_path: str = None):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), gridspec_kw={'height_ratios': [3, 1]})

        ax1.plot(positions, scores, 'b-', alpha=0.7, linewidth=0.5)
        ax1.fill_between(positions, 0, scores, alpha=0.3, color='blue')

        mean_score = np.mean(scores)
        std_score = np.std(scores)
        threshold = mean_score + std_score

        ax1.axhline(y=threshold, color='r', linestyle='--', alpha=0.7,
                    label=f'Threshold (μ+σ) = {threshold:.2f}')
        ax1.axhline(y=mean_score, color='g', linestyle='--', alpha=0.5,
                    label=f'Mean = {mean_score:.2f}')

        top_locations = self.find_top_motif_locations(positions, scores, top_n=10)
        top_positions = [pos for pos, score in top_locations]
        top_scores = [score for pos, score in top_locations]

        ax1.scatter(top_positions, top_scores, color='red', s=50, zorder=5,
                    label=f'Top {len(top_locations)} sites')

        ax1.set_xlabel('Genome Position (bp)', fontsize=12)
        ax1.set_ylabel('Motif Score (log-likelihood)', fontsize=12)
        ax1.set_title(f'Motif Signal Scan - {genome_id}', fontsize=14, fontweight='bold')
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)

        ax2.hist(scores, bins=50, alpha=0.7, color='blue', edgecolor='black')
        ax2.axvline(x=threshold, color='r', linestyle='--', alpha=0.7,
                    label=f'Threshold = {threshold:.2f}')
        ax2.set_xlabel('Score Distribution', fontsize=12)
        ax2.set_ylabel('Frequency', fontsize=12)
        ax2.set_title('Distribution of Motif Scores', fontsize=12)
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Plot saved to {save_path}")

        plt.show()

        if top_locations and sequence:
            print(f"\nTop {len(top_locations)} motif locations for {genome_id}:")
            print("-" * 80)
            print(f"{'Rank':<5} {'Position':<10} {'Score':<10} {'Sequence':<20} {'Context'}")
            print("-" * 80)

            for rank, (pos, score) in enumerate(top_locations, 1):
                motif_seq = sequence[pos:pos + self.motif_length]
                start_context = max(0, pos - 20)
                end_context = min(len(sequence), pos + self.motif_length + 20)
                context = sequence[start_context:end_context]

                motif_start_in_context = pos - start_context
                highlighted_context = (
                        context[:motif_start_in_context] +
                        "[" + context[motif_start_in_context:motif_start_in_context + self.motif_length] + "]" +
                        context[motif_start_in_context + self.motif_length:]
                )

                print(f"{rank:<5} {pos:<10} {score:<10.3f} {motif_seq:<20} {highlighted_context}")

    def analyze_genomes_from_fasta(self, fasta_filepath: str, output_dir: str = "output"):
        import os

        os.makedirs(output_dir, exist_ok=True)

        print(f"Reading genomes from {fasta_filepath}...")
        genomes = self.read_fasta_file(fasta_filepath)

        results = {}

        for genome_id, sequence in genomes.items():
            print(f"\n{'='*80}")
            print(f"Analyzing genome: {genome_id}")
            print(f"Sequence length: {len(sequence)} bp")
            print(f"{'='*80}")

            positions, scores = self.scan_genome(genome_id, sequence, self.log_matrix)

            mean_score = np.mean(scores)
            max_score = np.max(scores)
            threshold = mean_score + np.std(scores)

            sites_above_threshold = sum(1 for s in scores if s > threshold)

            print(f"Scan completed: {len(scores)} windows analyzed")
            print(f"Mean score: {mean_score:.3f}")
            print(f"Max score: {max_score:.3f}")
            print(f"Threshold (μ+σ): {threshold:.3f}")
            print(f"Potential motif sites above threshold: {sites_above_threshold}")

            plot_filename = f"{output_dir}/{genome_id}_motif_scan.png"
            plot_filename = re.sub(r'[^\w\-_\. ]', '_', plot_filename)

            self.plot_genome_signal(genome_id, positions, scores, sequence, plot_filename)

            results[genome_id] = {
                'positions': positions,
                'scores': scores,
                'mean_score': mean_score,
                'max_score': max_score,
                'threshold': threshold,
                'sites_above_threshold': sites_above_threshold,
                'top_locations': self.find_top_motif_locations(positions, scores, top_n=10)
            }

        self.generate_summary_report(results, output_dir)

        return results

    def generate_summary_report(self, results: Dict, output_dir: str):
        report_path = f"{output_dir}/motif_scan_summary.txt"

        with open(report_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("MOTIF SCAN SUMMARY REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"Motif Model Information:\n")
            f.write(f"- Motif length: {self.motif_length} bp\n")
            f.write(f"- Background frequency: {self.background_freq}\n")
            f.write(f"- Number of training sequences: {len(self.motif_sequences)}\n\n")

            f.write("Genome Analysis Results:\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"{'Genome ID':<30} {'Length':<10} {'Mean Score':<12} {'Max Score':<12} {'Sites > Thresh':<15}\n")
            f.write("-" * 80 + "\n")

            for genome_id, data in results.items():
                f.write(f"{genome_id:<30} {len(data['positions']):<10} {data['mean_score']:<12.3f} "
                        f"{data['max_score']:<12.3f} {data['sites_above_threshold']:<15}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("DETAILED ANALYSIS OF TOP MOTIF SITES\n")
            f.write("=" * 80 + "\n\n")

            for genome_id, data in results.items():
                f.write(f"\n{genome_id} - Top 5 Motif Locations:\n")
                f.write("-" * 80 + "\n")

                for rank, (pos, score) in enumerate(data['top_locations'][:5], 1):
                    f.write(f"  {rank}. Position: {pos}, Score: {score:.3f}\n")

        print(f"\nSummary report saved to {report_path}")

        self.plot_summary_comparison(results, output_dir)

    def plot_summary_comparison(self, results: Dict, output_dir: str):
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        genome_ids = list(results.keys())
        mean_scores = [results[gid]['mean_score'] for gid in genome_ids]
        max_scores = [results[gid]['max_score'] for gid in genome_ids]
        sites_above_thresh = [results[gid]['sites_above_threshold'] for gid in genome_ids]

        axes[0, 0].bar(range(len(genome_ids)), mean_scores, color='skyblue', edgecolor='black')
        axes[0, 0].set_xticks(range(len(genome_ids)))
        axes[0, 0].set_xticklabels(genome_ids, rotation=45, ha='right', fontsize=8)
        axes[0, 0].set_ylabel('Mean Score')
        axes[0, 0].set_title('Mean Motif Scores Across Genomes')
        axes[0, 0].grid(True, alpha=0.3, axis='y')

        axes[0, 1].bar(range(len(genome_ids)), max_scores, color='lightcoral', edgecolor='black')
        axes[0, 1].set_xticks(range(len(genome_ids)))
        axes[0, 1].set_xticklabels(genome_ids, rotation=45, ha='right', fontsize=8)
        axes[0, 1].set_ylabel('Max Score')
        axes[0, 1].set_title('Maximum Motif Scores Across Genomes')
        axes[0, 1].grid(True, alpha=0.3, axis='y')

        axes[1, 0].bar(range(len(genome_ids)), sites_above_thresh, color='lightgreen', edgecolor='black')
        axes[1, 0].set_xticks(range(len(genome_ids)))
        axes[1, 0].set_xticklabels(genome_ids, rotation=45, ha='right', fontsize=8)
        axes[1, 0].set_ylabel('Number of Sites')
        axes[1, 0].set_title('Potential Motif Sites (Above μ+σ Threshold)')
        axes[1, 0].grid(True, alpha=0.3, axis='y')

        all_scores = [results[gid]['scores'] for gid in genome_ids]
        axes[1, 1].boxplot(all_scores, labels=[gid[:15] + '...' for gid in genome_ids])
        axes[1, 1].set_xticklabels([gid[:15] + '...' for gid in genome_ids], rotation=45, ha='right', fontsize=8)
        axes[1, 1].set_ylabel('Score Distribution')
        axes[1, 1].set_title('Distribution of Motif Scores Across Genomes')
        axes[1, 1].grid(True, alpha=0.3, axis='y')

        plt.suptitle('Comparative Analysis of Motif Signals in Influenza Genomes', fontsize=16, fontweight='bold')
        plt.tight_layout()

        summary_plot_path = f"{output_dir}/summary_comparison.png"
        plt.savefig(summary_plot_path, dpi=150, bbox_inches='tight')
        plt.show()

        print(f"Summary comparison plot saved to {summary_plot_path}")

def main():

    motif_sequences = [
        "GAGGTAAAC",
        "TCCGTAAGT",
        "CAGGTTGGA",
        "ACAGTCAGT",
        "TAGGTCATT",
        "TAGGTACTG",
        "ATGGTAACT",
        "CAGGTATAC",
        "TGTGTGAGT",
        "AAGGTAAAT"
    ]

    print("=" * 80)
    print("GENOME MOTIF SCANNER")
    print("Scanning influenza genomes for potential exon-intron borders")
    print("=" * 80)
    print()

    scanner = GenomeMotifScanner(motif_sequences)

    fasta_filepath = input("Enter the path to your FASTA file (or press Enter to use example): ")

    if not fasta_filepath.strip():
        print("\n No file provided")

    try:
        results = scanner.analyze_genomes_from_fasta(fasta_filepath, output_dir="motif_scan_results")

        print("\n" + "=" * 80)
        print("OVERALL CONCLUSION")
        print("=" * 80)

        strongest_genome = max(results.items(), key=lambda x: x[1]['sites_above_threshold'])
        weakest_genome = min(results.items(), key=lambda x: x[1]['sites_above_threshold'])

        print(f"\nGenome with strongest motif signals: {strongest_genome[0]}")
        print(f"  - {strongest_genome[1]['sites_above_threshold']} potential motif sites")
        print(f"  - Maximum score: {strongest_genome[1]['max_score']:.3f}")

        print(f"\nGenome with weakest motif signals: {weakest_genome[0]}")
        print(f"  - {weakest_genome[1]['sites_above_threshold']} potential motif sites")
        print(f"  - Maximum score: {weakest_genome[1]['max_score']:.3f}")

        total_sites = sum(data['sites_above_threshold'] for data in results.values())
        print(f"\nTotal potential motif sites across all genomes: {total_sites}")

        if total_sites > 0:
            print("\n SIGNIFICANT SIGNALS DETECTED: Some genomes show potential exon-intron border motifs.")
            print("   Check the individual genome plots and summary report for details.")
        else:
            print("\n WEAK SIGNALS: Few potential motif sites detected across genomes.")
            print("   The motif model may need refinement or the genomes may have different characteristics.")

    except FileNotFoundError:
        print(f"\n ERROR: FASTA file not found at '{fasta_filepath}'")
        print("Please check the file path and try again.")
    except Exception as e:
        print(f"\n ERROR: {str(e)}")
        print("Please ensure your FASTA file is properly formatted.")

if __name__ == "__main__":
    try:
        import matplotlib
    except ImportError:
        print("Matplotlib is required for visualization.")
        print("Install it using: pip install matplotlib")
        import sys
        sys.exit(1)

    main()