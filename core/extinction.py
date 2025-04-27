import csv  # Import CSV module to read extinction coefficients from a CSV file
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
