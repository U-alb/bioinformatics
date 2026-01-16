import numpy as np
import json
import random

def generate_dna_sequence(length=50):
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def calculate_transition_matrix(sequence):
    bases = ['A', 'C', 'G', 'T']
    base_to_index = {base: i for i, base in enumerate(bases)}

    transition_counts = np.zeros((4, 4))

    for i in range(len(sequence) - 1):
        current_base = sequence[i]
        next_base = sequence[i + 1]
        current_idx = base_to_index[current_base]
        next_idx = base_to_index[next_base]
        transition_counts[current_idx][next_idx] += 1

    transition_matrix = np.zeros((4, 4))

    for i in range(4):
        row_sum = np.sum(transition_counts[i])
        if row_sum > 0:
            transition_matrix[i] = transition_counts[i] / row_sum
        else:
            transition_matrix[i] = np.ones(4) / 4

    return transition_matrix, bases

def create_initial_vector(sequence, bases):
    first_base = sequence[0]
    base_to_index = {base: i for i, base in enumerate(bases)}
    initial_vector = np.zeros(4)
    initial_vector[base_to_index[first_base]] = 1
    return initial_vector

def save_to_json(matrix, vector, bases, filename="dna_transition.json"):
    data = {
        "transition_matrix": matrix.tolist(),
        "initial_vector": vector.tolist(),
        "base_order": bases,
        "description": "DNA transition probabilities (A,C,G,T order)"
    }

    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)

def load_from_json(filename="dna_transition.json"):
    with open(filename, 'r') as f:
        data = json.load(f)

    matrix = np.array(data["transition_matrix"])
    vector = np.array(data["initial_vector"])
    bases = data["base_order"]

    return matrix, vector, bases

class DNAPredictor:
    def __init__(self, matrix, initial_vector, bases):
        self.matrix = np.array(matrix)
        self.state = np.array(initial_vector)
        self.bases = bases

    def predict_steps(self, num_steps=5):
        predictions = [self.state.copy()]

        for _ in range(num_steps):
            self.state = self.matrix @ self.state
            predictions.append(self.state.copy())

        return predictions

    def get_base_probabilities(self, state_vector):
        return {base: prob for base, prob in zip(self.bases, state_vector)}

    def print_predictions(self, num_steps=5):
        predictions = self.predict_steps(num_steps)

        print("DNA Sequence Predictions")
        print("Bases in order:", self.bases)
        print()

        for i, state in enumerate(predictions):
            if i == 0:
                print(f"Step {i} (Initial):")
            else:
                print(f"Step {i}:")

            probs = self.get_base_probabilities(state)
            for base, prob in probs.items():
                print(f"  {base}: {prob:.4f}")

            print()

def main():
    dna_sequence = generate_dna_sequence(50)
    print(f"Generated DNA sequence (50 bases):")
    print(dna_sequence)
    print()

    transition_matrix, bases = calculate_transition_matrix(dna_sequence)
    initial_vector = create_initial_vector(dna_sequence, bases)

    print("Transition Matrix (A,C,G,T order):")
    for i, row in enumerate(transition_matrix):
        print(f"{bases[i]}: {[round(x, 4) for x in row]}")
    print()

    print(f"Initial Vector (first base: '{dna_sequence[0]}'):")
    print([round(x, 4) for x in initial_vector])
    print()

    save_to_json(transition_matrix, initial_vector, bases)
    print("Saved to dna_transition.json")
    print()

    loaded_matrix, loaded_vector, loaded_bases = load_from_json()

    predictor = DNAPredictor(loaded_matrix, loaded_vector, loaded_bases)
    predictor.print_predictions(5)

    final_state = predictor.predict_steps(5)[-1]
    most_likely_base = loaded_bases[np.argmax(final_state)]
    print(f"After 5 steps, most likely base: {most_likely_base}")

    with open("dna_original_sequence.txt", "w") as f:
        f.write(dna_sequence)

if __name__ == "__main__":
    main()