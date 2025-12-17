import math

import numpy as np

class MotifFinder:
    def __init__(self):
        self.nucleotides = ['A', 'C', 'G', 'T']
        self.background_freq = 0.25  # P(A)=P(T)=P(C)=P(G)=0.25
        self.motif_length = 9

    def create_count_matrix(self, sequences):
        count_matrix = {nt: [0] * self.motif_length for nt in self.nucleotides}

        for i in range(self.motif_length):
            for seq in sequences:
                nt = seq[i]
                if nt in count_matrix:
                    count_matrix[nt][i] += 1

        return count_matrix

    def create_weight_matrix(self, count_matrix, pseudo_count=1):
        weight_matrix = {nt: [0] * self.motif_length for nt in self.nucleotides}
        num_sequences = sum(count_matrix[self.nucleotides[0]][:])

        for nt in self.nucleotides:
            for i in range(self.motif_length):
                weight_matrix[nt][i] = (count_matrix[nt][i] + pseudo_count) / (
                            num_sequences + len(self.nucleotides) * pseudo_count)

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

    def calculate_sequence_score(self, sequence, log_matrix):
        score = 0
        if len(sequence) < self.motif_length:
            return 0

        for i in range(self.motif_length):
            nt = sequence[i]
            if nt in log_matrix:
                score += log_matrix[nt][i]

        return score

    def scan_sequence(self, sequence, log_matrix):
        scores = []
        positions = []

        for i in range(len(sequence) - self.motif_length + 1):
            window = sequence[i:i + self.motif_length]
            score = self.calculate_sequence_score(window, log_matrix)
            scores.append(score)
            positions.append(i)

        return positions, scores

    def find_high_scoring_windows(self, sequence, log_matrix, threshold=0):
        positions, scores = self.scan_sequence(sequence, log_matrix)
        high_scoring = []

        for pos, score in zip(positions, scores):
            if score > threshold:
                high_scoring.append((pos, score, sequence[pos:pos + self.motif_length]))

        return high_scoring


def main():
    finder = MotifFinder()

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

    print("=== DNA Motif Finder Application ===\n")
    print("1. Count Matrix:")
    count_matrix = finder.create_count_matrix(motif_sequences)
    for nt in finder.nucleotides:
        print(f"   {nt}: {count_matrix[nt]}")

    print("\n2. Weight Matrix (with pseudocounts):")
    weight_matrix = finder.create_weight_matrix(count_matrix, pseudo_count=1)
    for nt in finder.nucleotides:
        print(f"   {nt}: [{', '.join([f'{x:.3f}' for x in weight_matrix[nt]])}]")

    print("\n3. Relative Frequencies Matrix:")
    rel_freq_matrix = finder.create_relative_frequencies(weight_matrix)
    for nt in finder.nucleotides:
        print(f"   {nt}: [{', '.join([f'{x:.3f}' for x in rel_freq_matrix[nt]])}]")

    print("\n4. Log-Likelihoods Matrix (natural log):")
    log_matrix = finder.create_log_likelihood_matrix(rel_freq_matrix)
    for nt in finder.nucleotides:
        print(f"   {nt}: [{', '.join([f'{x:.3f}' for x in log_matrix[nt]])}]")

    S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"
    print(f"\n5. Analyzing sequence S:")
    print(f"   S = '{S}'")
    print(f"   Length: {len(S)}")

    positions, scores = finder.scan_sequence(S, log_matrix)

    print(f"\n   Scores for each sliding window (position:score):")
    for pos, score in zip(positions, scores):
        window = S[pos:pos + finder.motif_length]
        print(f"   Position {pos:2d}: {window} - Score: {score:.3f}")

    print(f"\n   High-scoring windows (potential motif matches):")
    high_scoring = finder.find_high_scoring_windows(S, log_matrix, threshold=0)

    if high_scoring:
        for pos, score, window in high_scoring:
            print(f"   Position {pos:2d}: {window} - Score: {score:.3f}")

        print(f"\n   Analysis:")
        print(f"   Found {len(high_scoring)} potential motif sites with positive scores.")
        print(f"   Positive scores indicate that the sequence is more likely to match")
        print(f"   the motif model than a random sequence.")

        best_pos, best_score, best_window = max(high_scoring, key=lambda x: x[1])
        print(f"\n   Best candidate (position {best_pos}): {best_window}")
        print(f"   Score: {best_score:.3f}")

        print(f"\n   Known motifs have scores:")
        for i, motif in enumerate(motif_sequences, 1):
            motif_score = finder.calculate_sequence_score(motif, log_matrix)
            print(f"   Motif {i}: {motif} - Score: {motif_score:.3f}")

        avg_motif_score = np.mean([finder.calculate_sequence_score(m, log_matrix) for m in motif_sequences])
        print(f"\n   Average score of known motifs: {avg_motif_score:.3f}")

        if best_score > 0:
            print(f"\n SIGNAL DETECTED: The sequence S contains regions with positive")
            print(f"   log-likelihood scores, indicating potential exon-intron borders.")
            print(f"   The highest scoring window (position {best_pos}) has a score of {best_score:.3f},")
            print(f"   which suggests it may be a functional site.")
        else:
            print(f"\n WEAK SIGNAL: All scores are non-positive. The sequence S may not")
            print(f"   contain strong exon-intron border signals.")
    else:
        print(f"\n   NO SIGNAL DETECTED: No windows with positive scores found.")
        print(f"   The sequence S likely does not contain exon-intron borders.")


if __name__ == "__main__":
    main()
