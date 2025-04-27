import csv  # For reading base mass values from a CSV file
from core.modifications import get_modification

# Global dictionary to store monoisotopic mass of each nucleotide base
BASE_MASSES = {}

def load_base_masses(path="data/base_masses.csv"):
    global BASE_MASSES  # Required to update the global BASE_MASSES dictionary
    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)  # Read CSV rows as dictionaries
        for row in reader:
            BASE_MASSES[row["Base"]] = float(row["Mass"])  # Store each base's mass as float
        

def calculate_momo_iso_mass(seq):
    if not BASE_MASSES:
        # Load the base mass values if they haven't been loaded yet
        load_base_masses()
    return sum(BASE_MASSES.get(base, 0) for base in seq)

# def calculate_momo_iso_mass(seq):
#     if not BASE_MASSES:
#         load_base_masses()  # Load mass values if not already loaded
#     from core.modifications import get_modification, MODIFICATIONS
#     from core.sequence_utils import tokenize_sequence

#     tokens = tokenize_sequence(seq, MODIFICATIONS)
#     total = 0
#     print("mass_tokens:", tokens)
#     for token in tokens:
#             if token in BASE_MASSES:
#                 total += BASE_MASSES[token]
#             else:
#                 mod = get_modification(token)
#                 if mod:
#                     total += mod["mass"]
#                 else:
#                     raise ValueError(f"Unknown base/modification: {token}")
#     return total