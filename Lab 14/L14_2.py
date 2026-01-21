import numpy as np
import matplotlib.pyplot as plt
import re
from collections import defaultdict

class PoetryStyleAnalyzer:
    def __init__(self):
        self.chars = 'abcdefghijklmnopqrstuvwxyzăâîșț,.!?;- '
        self.char_to_idx = {char: i for i, char in enumerate(self.chars)}
        self.idx_to_char = {i: char for char, i in self.char_to_idx.items()}
        self.num_chars = len(self.chars)

        self.poet_models = {}
        self.log_likelihood_matrix = None

    def preprocess_text(self, text):
        text = text.lower()
        text = ''.join([c for c in text if c in self.char_to_idx])
        text = re.sub(r'\s+', ' ', text)
        return text.strip()

    def build_transition_matrix(self, text, smoothing=0.1):
        counts = np.ones((self.num_chars, self.num_chars)) * smoothing

        for i in range(len(text) - 1):
            from_char = text[i]
            to_char = text[i + 1]
            if from_char in self.char_to_idx and to_char in self.char_to_idx:
                from_idx = self.char_to_idx[from_char]
                to_idx = self.char_to_idx[to_char]
                counts[from_idx][to_idx] += 1

        row_sums = counts.sum(axis=1, keepdims=True)
        probabilities = counts / row_sums

        return probabilities

    def train_poet_model(self, poet_name, poems):
        combined_text = ""
        for poem in poems:
            combined_text += self.preprocess_text(poem) + " "

        print(f"Training {poet_name} model with {len(combined_text)} characters")
        transition_matrix = self.build_transition_matrix(combined_text)
        self.poet_models[poet_name] = transition_matrix
        return transition_matrix

    def create_log_likelihood_matrix(self):
        if 'Eminescu' not in self.poet_models or 'Stanescu' not in self.poet_models:
            raise ValueError("Both poet models must be trained first")

        eminescu_probs = self.poet_models['Eminescu']
        stanescu_probs = self.poet_models['Stanescu']

        epsilon = 1e-10
        eminescu_probs = np.maximum(eminescu_probs, epsilon)
        stanescu_probs = np.maximum(stanescu_probs, epsilon)

        log_likelihood = np.log2(eminescu_probs / stanescu_probs)
        self.log_likelihood_matrix = log_likelihood

        return log_likelihood

    def calculate_text_score(self, text, window_size=100, step_size=50):
        text = self.preprocess_text(text)
        text_length = len(text)

        if text_length < window_size:
            windows = [text]
            positions = [0]
        else:
            windows = []
            positions = []
            for start in range(0, text_length - window_size + 1, step_size):
                windows.append(text[start:start + window_size])
                positions.append(start)

        scores = []

        for window in windows:
            window_score = 0
            valid_transitions = 0

            for i in range(len(window) - 1):
                from_char = window[i]
                to_char = window[i + 1]

                if from_char in self.char_to_idx and to_char in self.char_to_idx:
                    from_idx = self.char_to_idx[from_char]
                    to_idx = self.char_to_idx[to_char]
                    window_score += self.log_likelihood_matrix[from_idx][to_idx]
                    valid_transitions += 1

            if valid_transitions > 0:
                scores.append(window_score / valid_transitions)
            else:
                scores.append(0)

        return positions, scores, windows

    def visualize_results(self, positions, scores, text_length, threshold=0):
        plt.figure(figsize=(12, 6))

        plt.plot(positions, scores, 'b-', linewidth=2, label='Style Score')
        plt.axhline(y=threshold, color='r', linestyle='--', label='Decision Threshold')
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)

        for i in range(len(positions) - 1):
            color = 'lightgreen' if scores[i] > threshold else 'lightcoral'
            alpha = 0.3 if scores[i] > threshold else 0.2
            plt.axvspan(positions[i], positions[i+1], alpha=alpha, color=color)

        plt.xlabel('Position in Text (characters)', fontsize=12)
        plt.ylabel('Style Score (log-likelihood ratio)', fontsize=12)
        plt.title('Poetry Style Analysis: Eminescu vs Stănescu', fontsize=14, fontweight='bold')

        plt.legend()

        plt.text(0.02, 0.98, 'Score > 0: More like Eminescu\nScore < 0: More like Stănescu',
                 transform=plt.gca().transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        return plt

eminescu_poems = [
    """Trecut-au anii ca nori lungi pe şesuri
    Şi niciodată n-or să vie iară,
    Căci nu mă-ncântă azi cum mă mişcară
    Poveşti şi doine, ghicitori, eresuri,.
    
    Ce fruntea-mi de copil o-nseninară,
    Abia-nţelese, pline de-nţelesuri -
    Cu-a tale umbre azi în van mă-mpesuri,
    O, ceas al tainei, asfinţit de sară.
    
    Să smulg un sunet din trecutul vieţii,
    Să fac, o, suflet, ca din nou să tremuri
    Cu mâna mea în van pe liră lunec;.
    
    Pierdut e totu-n zarea tinereţii
    Şi mută-i gura dulce-a altor vremuri,
    Iar timpul creşte-n urma mea.. mă-ntunec!""",

    """A fost odata ca-n povesti,
    A fost ca niciodata,
    Din rude mari imparatesti,
    O prea frumoasa fata.
    
    Si era una la părinti
    Si mindra-n toate cele,
    Cum e Fecioara intre sfinti
    Si luna intre stele."""
]

stanescu_poems = [
    """Noi ştim că unu ori unu fac unu,
    dar un inorog ori o pară
    nu ştim cât face.
    Ştim că cinci fără patru fac unu,
    dar un nor fără o corabie
    nu ştim cât face.
    Ştim, noi ştim că opt
    împărţit la opt fac unu,
    dar un munte împărţit la o capră
    nu ştim cât face.""",

    """Ea stă plictisită şi foarte frumoasă
    părul ei negru este supărat
    mâna ei luminoasă
    demult m-a uitat, -
    demult s-a uitat şi pe sine
    cum atârnă pe ceafa scaunului.
    
    Eu mă înec în lumine
    şi scrişnesc în crugul anului.
    Îi arăt dinţii din gură,
    dar ea ştie că eu nu râd,
    dulcea luminii faptură
    mie, pe mine mă înfăţişează pe când
    ea stă plictisită şi foarte frumoasa
    şi eu numai pentru ea trîăiesc
    în lumea fioroasă
    de sub ceresc.""",

    """Ea devenise încetul cu încetul cuvânt,
    fuioare de suflet de vânt,
    delfin în ghearele sprâncenelor mele,
    piatră stârnind în apă inele,
    stea înlăuntrul genunchiului meu,
    cer înlăuntrul umărului meu,
    eu înlăuntrul eului meu."""
    ]

test_text = """
Trecut-au anii ca nori lungi pe şesuri,
A fost odata ca-n povesti, a fost ca niciodata;
Şi niciodată n-or să vie iară,
Din rude mari imparatesti, o prea frumoasa fata.

Ce fruntea-mi de copil o-nseninară,
Si era una la părinti şi mindra-n toate cele,
Abia-nţelese, pline de-nţelesuri —
Cum e Fecioara intre sfinti şi luna intre stele.

Să smulg un sunet din trecutul vieţii,
Noi ştim că unu ori unu fac unu,
Să fac, o, suflet, ca din nou să tremuri
dar un inorog ori o pară nu ştim cât face.

Pierdut e totu-n zarea tinereţii,
Ştim că cinci fără patru fac unu,
Şi mută-i gura dulce-a altor vremuri—
dar un nor fără o corabie nu ştim cât face.

Ea stă plictisită şi foarte frumoasă,
Ce fruntea-mi de copil o-nseninară;  
părul ei negru este supărat,
Abia-nţelese, pline de-nţelesuri.

Eu mă înec în lumine şi scrişnesc în crugul anului;
Cu mâna mea în van pe liră lunec;
Îi arăt dinţii din gură, dar ea ştie că eu nu râd,
dulcea luminii faptură mie, pe mine mă înfăţişează.

Ea devenise încetul cu încetul cuvânt,
fuioare de suflet de vânt,
delfin în ghearele sprâncenelor mele,
piatră stârnind în apă inele.

Să smulg un sunet din trecutul vieţii,
stea înlăuntrul genunchiului meu,
cer înlăuntrul umărului meu,
eu înlăuntrul eului meu.

Iar timpul creşte-n urma mea... mă-ntunec!
Ştim, noi ştim că opt împărţit la opt fac unu,
dar un munte împărţit la o capră nu ştim cât face —
O, ceas al tainei, asfinţit de sară.
"""

def main():
    print("=" * 60)
    print("POETRY STYLE ANALYSIS: Eminescu vs Stănescu")
    print("=" * 60)

    analyzer = PoetryStyleAnalyzer()

    print("\n1. Training poet models...")
    analyzer.train_poet_model('Eminescu', eminescu_poems)
    analyzer.train_poet_model('Stanescu', stanescu_poems)

    print("\n2. Creating log-likelihood matrix...")
    log_matrix = analyzer.create_log_likelihood_matrix()
    print(f"Log-likelihood matrix shape: {log_matrix.shape}")

    print("\nSample log-likelihood values (Eminescu vs Stănescu):")
    sample_chars = ['a', 'e', 'o', 't', ' ']
    for from_char in sample_chars:
        for to_char in sample_chars:
            if from_char in analyzer.char_to_idx and to_char in analyzer.char_to_idx:
                from_idx = analyzer.char_to_idx[from_char]
                to_idx = analyzer.char_to_idx[to_char]
                value = log_matrix[from_idx][to_idx]
                print(f"  {from_char}->{to_char}: {value:.3f}")

    print(f"\n3. Analyzing test text (length: {len(test_text)} characters)...")
    print(f"Test text preview: {test_text[:200]}...")

    positions, scores, windows = analyzer.calculate_text_score(
        test_text, window_size=50, step_size=25
    )

    print(f"\n4. Analysis Results:")
    print(f"   Number of windows analyzed: {len(scores)}")
    print(f"   Window size: 50 characters")
    print(f"   Step size: 25 characters")

    eminescu_windows = sum(1 for score in scores if score > 0)
    stanescu_windows = len(scores) - eminescu_windows

    print(f"\n   Classification:")
    print(f"   Windows classified as Eminescu-like: {eminescu_windows} ({eminescu_windows/len(scores)*100:.1f}%)")
    print(f"   Windows classified as Stănescu-like: {stanescu_windows} ({stanescu_windows/len(scores)*100:.1f}%)")

    print(f"\n   Example window classifications:")
    for i in range(min(3, len(windows))):
        poet = "Eminescu" if scores[i] > 0 else "Stănescu"
        print(f"   Window starting at position {positions[i]}: {poet} (score: {scores[i]:.3f})")
        print(f"     Text: '{windows[i][:30]}...'")

    print(f"\n5. Generating visualization...")
    plt = analyzer.visualize_results(positions, scores, len(test_text))

    plt.savefig('poetry_style_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()

    print(f"\n" + "=" * 60)
    print("STYLISTIC PATTERNS ANALYSIS")
    print("=" * 60)

    comparison_chars = [('e', 'm'), ('a', ' '), ('o', 'n'), (' ', 't')]

    print("\nTransition probability comparison:")
    print("Char Pair | Eminescu Prob | Stănescu Prob | Difference")
    print("-" * 50)

    for from_char, to_char in comparison_chars:
        if from_char in analyzer.char_to_idx and to_char in analyzer.char_to_idx:
            from_idx = analyzer.char_to_idx[from_char]
            to_idx = analyzer.char_to_idx[to_char]

            em_prob = analyzer.poet_models['Eminescu'][from_idx][to_idx]
            st_prob = analyzer.poet_models['Stanescu'][from_idx][to_idx]
            diff = em_prob - st_prob

            print(f"{from_char}->{to_char}  | {em_prob:.4f}       | {st_prob:.4f}       | {diff:+.4f}")

    return analyzer, positions, scores

if __name__ == "__main__":
    analyzer, positions, scores = main()