def normalize_transitions(counts):
    smoothed_counts = [[count + 1 for count in row] for row in counts]

    probabilities = [[0] * 4 for _ in range(4)]

    for i in range(4):
        row_total = sum(smoothed_counts[i])
        for j in range(4):
            probabilities[i][j] = smoothed_counts[i][j] / row_total

    return probabilities


def log2(x):
    import math
    return math.log2(x) if x > 0 else float('-inf')


def create_log_likelihood_matrix(pos_probs, neg_probs):
    log_matrix = [[0] * 4 for _ in range(4)]

    for i in range(4):
        for j in range(4):
            if pos_probs[i][j] > 0 and neg_probs[i][j] > 0:
                log_matrix[i][j] = log2(pos_probs[i][j] / neg_probs[i][j])
            else:
                eps = 1e-10
                ratio = (pos_probs[i][j] + eps) / (neg_probs[i][j] + eps)
                log_matrix[i][j] = log2(ratio)

    return log_matrix


class CpGIslandDetector:
    def __init__(self):
        self.nucleotides = ['A', 'C', 'G', 'T']
        self.nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    def count_transitions(self, sequence):
        counts = [[0] * 4 for _ in range(4)]

        for i in range(len(sequence) - 1):
            from_nuc = sequence[i]
            to_nuc = sequence[i + 1]
            from_idx = self.nuc_to_idx[from_nuc]
            to_idx = self.nuc_to_idx[to_nuc]
            counts[from_idx][to_idx] += 1

        return counts

    def calculate_sequence_score(self, sequence, log_matrix):
        score = 0

        for i in range(len(sequence) - 1):
            from_nuc = sequence[i]
            to_nuc = sequence[i + 1]
            from_idx = self.nuc_to_idx[from_nuc]
            to_idx = self.nuc_to_idx[to_nuc]
            score += log_matrix[from_idx][to_idx]

        return score

    def display_matrix(self, matrix, title):
        print(f"\n{title}:")
        print("    A       C       G       T")
        for i, row_label in enumerate(self.nucleotides):
            row_str = f"{row_label} "
            for val in matrix[i]:
                if isinstance(val, float):
                    row_str += f"{val:7.3f} "
                else:
                    row_str += f"{val:7} "
            print(row_str)

    def run_analysis(self):
        S1 = "ATCGATTCGATATCATACACGTAT"
        S2 = "CTCGACTAGTATGAAGTCCACGCTTG"
        S_test = "CAGGTTGGAAACGTAA"

        print("=" * 60)
        print("CpG Island Detection using Log-Likelihood Matrix")
        print("=" * 60)

        print("\n1. Counting transitions for CpG+ model (from S1):")
        pos_counts = self.count_transitions(S1)
        print(f"Sequence S1: {S1}")
        self.display_matrix(pos_counts, "Transition Counts (CpG+)")

        print("\n2. Counting transitions for CpG- model (from S2):")
        neg_counts = self.count_transitions(S2)
        print(f"Sequence S2: {S2}")
        self.display_matrix(neg_counts, "Transition Counts (CpG-)")

        print("\n3. Normalizing counts to get transition probabilities:")
        pos_probs = normalize_transitions(pos_counts)
        self.display_matrix(pos_probs, "Transition Probabilities (CpG+)")

        neg_probs = normalize_transitions(neg_counts)
        self.display_matrix(neg_probs, "Transition Probabilities (CpG-)")

        print("\n4. Creating log-likelihood matrix (log2(P+/P-)):")
        log_matrix = create_log_likelihood_matrix(pos_probs, neg_probs)
        self.display_matrix(log_matrix, "Log-Likelihood Matrix")

        print("\n5. Testing sequence S:")
        print(f"Test sequence S: {S_test}")

        score = self.calculate_sequence_score(S_test, log_matrix)
        print(f"Total log-likelihood score: {score:.3f}")

        if score > 0:
            print(f"Decision: S BELONGS to a CpG island (score > 0)")
        else:
            print(f"Decision: S DOES NOT BELONG to a CpG island (score ≤ 0)")

        print("\n" + "=" * 60)
        print("Additional Analysis:")
        print("=" * 60)

        score_S1 = self.calculate_sequence_score(S1, log_matrix)
        score_S2 = self.calculate_sequence_score(S2, log_matrix)

        print(f"Score for S1 (known CpG island): {score_S1:.3f}")
        print(f"Score for S2 (known non-island): {score_S2:.3f}")
        print(f"Score for S (test sequence): {score:.3f}")

        print("\n" + "=" * 60)
        print("Comparison with Provided Matrix from Image:")
        print("=" * 60)

        provided_matrix = [
            [0, -0.415, 0, 2.17],  # A->A, A->C, A->G, A->T
            [1.485, 0, 1.07, 0],  # C->A, C->C, C->G, C->T
            [1, 0, 0, 0],  # G->A, G->C, G->G, G->T
            [0.392, 0.392, 0, -0.193]  # T->A, T->C, T->G, T->T
        ]

        print("Using provided log-likelihood matrix from image:")
        score_provided = self.calculate_sequence_score(S_test, provided_matrix)
        print(f"Total score using provided matrix: {score_provided:.3f}")

        if score_provided > 0:
            print(f"Decision using provided matrix: S BELONGS to a CpG island")
        else:
            print(f"Decision using provided matrix: S DOES NOT BELONG to a CpG island")

        return {
            'pos_probs': pos_probs,
            'neg_probs': neg_probs,
            'log_matrix': log_matrix,
            'score': score,
            'score_S1': score_S1,
            'score_S2': score_S2,
            'score_provided': score_provided
        }


def main():
    detector = CpGIslandDetector()
    results = detector.run_analysis()

    print("\n" + "=" * 60)
    print("Summary:")
    print("=" * 60)
    print(f"Test sequence score (from our model): {results['score']:.3f}")
    print(f"Test sequence score (from provided matrix): {results['score_provided']:.3f}")

    print("\n" + "=" * 60)
    print("CONCLUSION:")
    print("=" * 60)

    print("Based on OUR model built from S1 and S2, the test sequence")
    if results['score'] > 0:
        print("likely BELONGS to a CpG island (positive log-likelihood score).")
    else:
        print("likely DOES NOT BELONG to a CpG island (non-positive score).")

    print(f"\nReference: Known CpG island (S1) score: {results['score_S1']:.3f}")
    print(f"Reference: Known non-island (S2) score: {results['score_S2']:.3f}")


if __name__ == "__main__":
    main()
