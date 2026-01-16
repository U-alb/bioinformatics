import json
import random
import re
from collections import defaultdict, Counter

def generate_english_text():
    words = ["the", "quick", "brown", "fox", "jumps", "over", "lazy", "dog",
             "sun", "shines", "bright", "morning", "birds", "sing", "sweet",
             "air", "fresh", "cool", "wind", "blows", "gentle", "breeze",
             "river", "flows", "calm", "peaceful", "nature", "beautiful",
             "world", "full", "wonder", "mystery", "life", "simple", "joy",
             "happiness", "love", "kindness", "peace", "harmony", "balance"]

    text_words = []
    current_length = 0
    target_length = 300

    while current_length < target_length:
        word = random.choice(words)
        if current_length + len(word) + (1 if text_words else 0) <= target_length:
            text_words.append(word)
            current_length += len(word) + (1 if text_words else 0)
        else:
            remaining = target_length - current_length
            if remaining >= 3:
                short_word = random.choice([w for w in words if len(w) <= remaining])
                text_words.append(short_word)
                current_length += len(short_word) + (1 if text_words else 0)
            break

    text = " ".join(text_words)
    return text[:300]

def process_text_to_words(text):
    text = re.sub(r'[^\w\s]', '', text.lower())
    words = text.split()
    return words

def create_word_mapping(words):
    word_to_symbol = {}
    symbol_to_word = {}
    symbol_counter = 1

    for word in words:
        if word not in word_to_symbol:
            symbol = f"W{symbol_counter}"
            word_to_symbol[word] = symbol
            symbol_to_word[symbol] = word
            symbol_counter += 1

    return word_to_symbol, symbol_to_word

def calculate_word_transitions(words, word_to_symbol):
    symbols = list(word_to_symbol.values())
    symbol_index = {symbol: idx for idx, symbol in enumerate(symbols)}

    n = len(symbols)
    transition_counts = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(len(words) - 1):
        current_word = words[i]
        next_word = words[i + 1]

        current_symbol = word_to_symbol[current_word]
        next_symbol = word_to_symbol[next_word]

        current_idx = symbol_index[current_symbol]
        next_idx = symbol_index[next_symbol]

        transition_counts[current_idx][next_idx] += 1

    transition_matrix = []
    for i in range(n):
        row_sum = sum(transition_counts[i])
        if row_sum > 0:
            normalized_row = [count / row_sum for count in transition_counts[i]]
        else:
            normalized_row = [1/n for _ in range(n)]
        transition_matrix.append(normalized_row)

    return transition_matrix, symbols, symbol_index

def create_initial_vector(words, symbols, symbol_index):
    first_word = words[0]
    word_to_symbol, _ = create_word_mapping([first_word])
    first_symbol = word_to_symbol[first_word]

    initial_vector = [0 for _ in range(len(symbols))]
    if first_symbol in symbol_index:
        initial_vector[symbol_index[first_symbol]] = 1

    return initial_vector

def save_transition_data(transition_matrix, initial_vector, symbols, symbol_to_word, filename="word_transitions.json"):
    data = {
        "transition_matrix": transition_matrix,
        "initial_vector": initial_vector,
        "symbols": symbols,
        "symbol_to_word": symbol_to_word,
        "num_states": len(symbols)
    }

    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)

class WordTransitionPredictor:
    def __init__(self, matrix, initial_vector, symbols, symbol_to_word):
        self.matrix = matrix
        self.state = initial_vector
        self.symbols = symbols
        self.symbol_to_word = symbol_to_word
        self.symbol_index = {symbol: idx for idx, symbol in enumerate(symbols)}

    def predict_steps(self, num_steps=5):
        predictions = [self.state.copy()]

        for _ in range(num_steps):
            next_state = [0 for _ in range(len(self.symbols))]

            for i in range(len(self.symbols)):
                for j in range(len(self.symbols)):
                    next_state[j] += self.state[i] * self.matrix[i][j]

            self.state = next_state
            predictions.append(next_state.copy())

        return predictions

    def get_word_probabilities(self, state_vector):
        probabilities = {}
        for symbol, prob in zip(self.symbols, state_vector):
            word = self.symbol_to_word[symbol]
            probabilities[word] = prob
        return probabilities

def main():
    text = generate_english_text()
    print("Generated English Text (300 letters):")
    print(text)
    print()

    words = process_text_to_words(text)
    print(f"Total words: {len(words)}")
    print(f"Unique words: {len(set(words))}")
    print(f"Words list: {words}")
    print()

    word_to_symbol, symbol_to_word = create_word_mapping(words)
    print("Word to Symbol Mapping:")
    for word, symbol in word_to_symbol.items():
        print(f"  {word}: {symbol}")
    print()

    transition_matrix, symbols, symbol_index = calculate_word_transitions(words, word_to_symbol)
    initial_vector = create_initial_vector(words, symbols, symbol_index)

    print("Transition Matrix Size:", len(transition_matrix), "x", len(transition_matrix[0]))
    print("First few rows of transition matrix:")
    for i in range(min(5, len(transition_matrix))):
        print(f"  {symbols[i]}: {[round(p, 3) for p in transition_matrix[i][:5]]}...")
    print()

    save_transition_data(transition_matrix, initial_vector, symbols, symbol_to_word)
    print("Saved to word_transitions.json")
    print()

    with open("word_transitions.json", 'r') as f:
        loaded_data = json.load(f)

    predictor = WordTransitionPredictor(
        loaded_data["transition_matrix"],
        loaded_data["initial_vector"],
        loaded_data["symbols"],
        loaded_data["symbol_to_word"]
    )

    predictions = predictor.predict_steps(5)
    print("Word Probability Predictions for 5 steps:")
    print("-" * 40)

    for step, state in enumerate(predictions):
        print(f"\nStep {step}:")
        word_probs = predictor.get_word_probabilities(state)

        sorted_probs = sorted(word_probs.items(), key=lambda x: x[1], reverse=True)[:5]

        for word, prob in sorted_probs:
            if prob > 0.001:
                print(f"  {word}: {prob:.4f}")

        if step < len(predictions) - 1:
            most_likely_word = max(word_probs.items(), key=lambda x: x[1])[0]
            print(f"  Most likely next word: {most_likely_word}")

    with open("original_text.txt", "w") as f:
        f.write(text)

if __name__ == "__main__":
    main()