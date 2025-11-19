import re

# Sample DNA Sequence
dna_sequence = """ATGCTAGCTAGCTTAGCTAGCTGATGCGTAGCTAGCGTACGTAGCTGCTAGCTAGCTAGC
ATCGATGCACTACGTAGCTGATGC
TACGATCGTACGTAGCTAGCTGTCG
TAGCTAGCGTATCGATGTTAGCTAG
CTGATCGTAGCTAGTAGCTAGCCCCCA
"""

transposable_elements = [
    "ATCGATGCACTACGTAGCTGATGC",
    "TACGATCGTACGTAGCTAGCTGTCG",
    "TAGCTAGCGTATCGATGTTAGCTAG"
]

def find_transposable_elements(dna_sequence, elements):
    positions = {}
    for element in elements:
        pattern = re.escape(element)
        matches = list(re.finditer(pattern, dna_sequence))
        positions[element] = [(match.start(), match.end()) for match in matches]

    return positions

# Locate transposable elements
positions = find_transposable_elements(dna_sequence, transposable_elements)

# Output results
for element, pos in positions.items():
    for start, end in pos:
        print(f"Transposable Element: {element} | Start: {start} | End: {end}")
