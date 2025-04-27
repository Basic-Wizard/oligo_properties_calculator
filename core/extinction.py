import csv  # Import CSV module to read extinction coefficients from a CSV file
from core.modifications import get_modification, MODIFICATIONS
#from core.sequence_utils import tokenize_sequence


# Global dictionaries to store extinction coefficients for single bases and base pairs
EXTINCTION_BASES = {}
EXTINCTION_PAIRS = {}

def load_extinction_coeffs(path="data/extinction_coeffs.csv"):
    global EXTINCTION_BASES, EXTINCTION_PAIRS  # Needed to modify global dictionaries
    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)  # Reads the CSV file into a dictionary per row
        for row in reader:
            # Store extinction coefficient based on the type (Base or Pair)
            if row["Type"] == "Base":
                EXTINCTION_BASES[row["Sequence"]] = int(row["Value"])
            elif row["Type"] == "Pair":
                EXTINCTION_PAIRS[row["Sequence"]] = int(row["Value"])

def calculate_ext(seq):
    # Ensure coefficients are loaded before calculation
    if not EXTINCTION_BASES or not EXTINCTION_PAIRS:
        load_extinction_coeffs()

    total = 0  # Running total of the extinction coefficient

    # Loop through sequence and sum extinction values of each base pair
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]  # Get overlapping 2-base pair
        total += EXTINCTION_PAIRS.get(pair, 0)  # Add its coefficient (0 if not found)

    # Add single-base contributions from the 5' and 3' ends
    if seq:
        total += EXTINCTION_BASES.get(seq[0], 0)  # First base
        total += EXTINCTION_BASES.get(seq[-1], 0)  # Last base

    return total  # Final extinction coefficient for the full sequence


    
#     if not MODIFICATIONS:
#         from core.modifications import load_modifications
#         load_modifications()

#     # Sum extinction values of adjacent base pairs in the sequence
#     for i in range(len(tokens) - 1):
#         pair = seq[i:i+2]
#         a, b = tokens[i], tokens[i + 1]
#         if a in EXTINCTION_BASES and b in EXTINCTION_BASES:
#             pair = a + b
#             total += EXTINCTION_PAIRS.get(pair, 0)

#     if tokens:
#         first = tokens[0]
#         last = tokens[-1]
    
#         total += get_ext_coeff_for_token(first)
#         if last != first:  # Avoid double-counting if length 1
#             total += get_ext_coeff_for_token(last)
    
#     for token in tokens[1:-1]:
#         if token not in EXTINCTION_BASES:
#             mod = get_modification(token)
#             if mod:
#                 total += mod["ext"]

#     return total

# def get_ext_coeff_for_token(token):
#     if token in EXTINCTION_BASES:
#         return EXTINCTION_BASES[token]
#     mod = get_modification(token)
#     return mod["ext"] if mod else 0
