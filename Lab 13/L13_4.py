import json
import numpy as np
import random


class DNAGenerator:
    def __init__(self, json_file="dna_transition.json"):
        with open(json_file, 'r') as f:
            data = json.load(f)

        self.matrix = np.array(data["transition_matrix"])
        self.initial_vector = np.array(data["initial_vector"])
        self.bases = data["base_order"]
        self.base_index = {base: i for i, base in enumerate(self.bases)}
        self.current_state = self.initial_vector.copy()

    def reset(self):
        self.current_state = self.initial_vector.copy()

    def get_next_base(self):
        probabilities = list(self.current_state)
        next_base = random.choices(self.bases, weights=probabilities)[0]
        next_index = self.base_index[next_base]
        self.current_state = self.matrix[next_index]
        return next_base

    def generate_sequence(self, length=50):
        sequence = []
        self.reset()

        first_index = np.argmax(self.initial_vector)
        first_base = self.bases[first_index]
        sequence.append(first_base)
        self.current_state = self.matrix[first_index]

        for _ in range(length - 1):
            base = self.get_next_base()
            sequence.append(base)

        return ''.join(sequence)

    def generate_markov_sequence(self, length=50):
        sequence = []
        self.reset()

        first_index = np.argmax(self.initial_vector)
        first_base = self.bases[first_index]
        sequence.append(first_base)
        current_index = first_index

        for _ in range(length - 1):
            next_probs = self.matrix[current_index]
            next_base = random.choices(self.bases, weights=next_probs)[0]
            sequence.append(next_base)
            current_index = self.base_index[next_base]

        return ''.join(sequence)


class TextGenerator:
    def __init__(self, json_file="word_transitions.json"):
        with open(json_file, 'r') as f:
            data = json.load(f)

        self.matrix = np.array(data["transition_matrix"])
        self.initial_vector = np.array(data["initial_vector"])
        self.symbols = data["symbols"]
        self.symbol_to_word = data["symbol_to_word"]
        self.symbol_index = {symbol: i for i, symbol in enumerate(self.symbols)}
        self.current_state = self.initial_vector.copy()

    def reset(self):
        self.current_state = self.initial_vector.copy()

    def get_next_symbol(self):
        probabilities = list(self.current_state)
        next_symbol = random.choices(self.symbols, weights=probabilities)[0]
        next_index = self.symbol_index[next_symbol]
        self.current_state = self.matrix[next_index]
        return next_symbol

    def generate_text(self, num_words=20):
        words = []
        self.reset()

        first_index = np.argmax(self.initial_vector)
        first_symbol = self.symbols[first_index]
        words.append(self.symbol_to_word[first_symbol])
        self.current_state = self.matrix[first_index]

        for _ in range(num_words - 1):
            symbol = self.get_next_symbol()
            words.append(self.symbol_to_word[symbol])

        return ' '.join(words)

    def generate_markov_text(self, num_words=20):
        words = []
        self.reset()

        first_index = np.argmax(self.initial_vector)
        first_symbol = self.symbols[first_index]
        words.append(self.symbol_to_word[first_symbol])
        current_index = first_index

        for _ in range(num_words - 1):
            next_probs = self.matrix[current_index]
            next_symbol = random.choices(self.symbols, weights=next_probs)[0]
            words.append(self.symbol_to_word[next_symbol])
            current_index = self.symbol_index[next_symbol]

        return ' '.join(words)

    def generate_with_beam_search(self, num_words=20, beam_width=3):
        beams = []

        first_index = np.argmax(self.initial_vector)
        first_symbol = self.symbols[first_index]
        first_word = self.symbol_to_word[first_symbol]

        initial_beam = {
            'sequence': [first_word],
            'symbols': [first_symbol],
            'probability': 1.0,
            'current_index': first_index
        }
        beams = [initial_beam]

        for step in range(1, num_words):
            new_beams = []

            for beam in beams:
                current_idx = beam['current_index']
                probs = self.matrix[current_idx]

                top_indices = np.argsort(probs)[-beam_width:][::-1]

                for idx in top_indices:
                    if probs[idx] > 0:
                        new_symbol = self.symbols[idx]
                        new_word = self.symbol_to_word[new_symbol]
                        new_prob = beam['probability'] * probs[idx]

                        new_sequence = beam['sequence'] + [new_word]
                        new_symbol_seq = beam['symbols'] + [new_symbol]

                        new_beams.append({
                            'sequence': new_sequence,
                            'symbols': new_symbol_seq,
                            'probability': new_prob,
                            'current_index': idx
                        })

            new_beams.sort(key=lambda x: x['probability'], reverse=True)
            beams = new_beams[:beam_width]

        best_beam = max(beams, key=lambda x: x['probability'])
        return ' '.join(best_beam['sequence'])


class TrainingDataCreator:
    @staticmethod
    def create_dna_training_data():
        bases = ['A', 'C', 'G', 'T']
        sequence = ''.join(random.choices(bases, k=50))

        base_index = {base: i for i, base in enumerate(bases)}
        transition_counts = np.zeros((4, 4))

        for i in range(len(sequence) - 1):
            current_idx = base_index[sequence[i]]
            next_idx = base_index[sequence[i + 1]]
            transition_counts[current_idx][next_idx] += 1

        transition_matrix = np.zeros((4, 4))
        for i in range(4):
            row_sum = np.sum(transition_counts[i])
            if row_sum > 0:
                transition_matrix[i] = transition_counts[i] / row_sum
            else:
                transition_matrix[i] = np.ones(4) / 4

        initial_vector = np.zeros(4)
        initial_vector[base_index[sequence[0]]] = 1

        data = {
            "transition_matrix": transition_matrix.tolist(),
            "initial_vector": initial_vector.tolist(),
            "base_order": bases,
            "original_sequence": sequence
        }

        with open("dna_training.json", 'w') as f:
            json.dump(data, f, indent=2)

        return sequence

    @staticmethod
    def create_text_training_data():
        words = ["the", "cat", "sat", "on", "mat", "dog", "ran", "fast",
                 "quick", "brown", "fox", "jumps", "over", "lazy", "moon",
                 "sun", "shines", "bright", "stars", "twinkle", "night"]

        text_words = []
        for _ in range(30):
            text_words.append(random.choice(words))

        text = ' '.join(text_words)

        word_to_symbol = {}
        symbol_to_word = {}
        symbol_counter = 1

        for word in text_words:
            if word not in word_to_symbol:
                symbol = f"W{symbol_counter}"
                word_to_symbol[word] = symbol
                symbol_to_word[symbol] = word
                symbol_counter += 1

        symbols = list(word_to_symbol.values())
        symbol_index = {symbol: idx for idx, symbol in enumerate(symbols)}

        n = len(symbols)
        transition_counts = np.zeros((n, n))

        for i in range(len(text_words) - 1):
            current_symbol = word_to_symbol[text_words[i]]
            next_symbol = word_to_symbol[text_words[i + 1]]
            current_idx = symbol_index[current_symbol]
            next_idx = symbol_index[next_symbol]
            transition_counts[current_idx][next_idx] += 1

        transition_matrix = np.zeros((n, n))
        for i in range(n):
            row_sum = np.sum(transition_counts[i])
            if row_sum > 0:
                transition_matrix[i] = transition_counts[i] / row_sum
            else:
                transition_matrix[i] = np.ones(n) / n

        first_symbol = word_to_symbol[text_words[0]]
        initial_vector = np.zeros(n)
        initial_vector[symbol_index[first_symbol]] = 1

        data = {
            "transition_matrix": transition_matrix.tolist(),
            "initial_vector": initial_vector.tolist(),
            "symbols": symbols,
            "symbol_to_word": symbol_to_word,
            "original_text": text
        }

        with open("text_training.json", 'w') as f:
            json.dump(data, f, indent=2)

        return text


def main():
    print("=" * 60)
    print("SEQUENCE GENERATION ENGINE")
    print("=" * 60)

    while True:
        print("\n1. Generate DNA Sequences")
        print("2. Generate English Text")
        print("3. Create Training Data")
        print("4. Exit")

        choice = input("\nSelect option (1-4): ").strip()

        if choice == '1':
            try:
                dna_gen = DNAGenerator("dna_training.json")

                print("\nDNA Sequence Generation")
                print("-" * 40)

                length = int(input("Enter sequence length: "))
                num_sequences = int(input("Enter number of sequences to generate: "))

                print("\nGenerated DNA Sequences:")
                for i in range(num_sequences):
                    seq = dna_gen.generate_markov_sequence(length)
                    print(f"{i + 1}: {seq}")

                with open("dna_training.json", 'r') as f:
                    training_data = json.load(f)
                print(f"\nOriginal training sequence: {training_data['original_sequence']}")

            except FileNotFoundError:
                print("DNA training file not found. Create training data first.")

        elif choice == '2':
            try:
                text_gen = TextGenerator("text_training.json")

                print("\nText Generation")
                print("-" * 40)

                num_words = int(input("Enter number of words: "))
                num_texts = int(input("Enter number of texts to generate: "))

                print("\n1. Random sampling")
                print("2. Beam search (higher quality)")
                method = input("Select generation method: ").strip()

                print("\nGenerated Texts:")
                for i in range(num_texts):
                    if method == '2':
                        text = text_gen.generate_with_beam_search(num_words)
                    else:
                        text = text_gen.generate_markov_text(num_words)
                    print(f"{i + 1}: {text}")

                with open("text_training.json", 'r') as f:
                    training_data = json.load(f)
                print(f"\nOriginal training text: {training_data['original_text']}")

            except FileNotFoundError:
                print("Text training file not found. Create training data first.")

        elif choice == '3':
            print("\nCreating Training Data")
            print("-" * 40)

            dna_seq = TrainingDataCreator.create_dna_training_data()
            print(f"Created DNA training data")
            print(f"DNA sequence: {dna_seq}")

            text_seq = TrainingDataCreator.create_text_training_data()
            print(f"\nCreated Text training data")
            print(f"Text: {text_seq}")

            print("\nTraining data saved to:")
            print("- dna_training.json")
            print("- text_training.json")

        elif choice == '4':
            print("Exiting...")
            break

        else:
            print("Invalid choice. Please try again.")


if __name__ == "__main__":
    main()
