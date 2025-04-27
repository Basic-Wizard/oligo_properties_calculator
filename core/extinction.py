import csv  # For reading extinction coefficients from a CSV file

# Global dictionaries to store extinction coefficients for single bases and base pairs
EXTINCTION_BASES = {}
EXTINCTION_PAIRS = {}


def load_extinction_coeffs(path="data/extinction_coeffs.csv"):
    """
    Load extinction coefficients from a CSV file into global dictionaries.

    Args:
        path (str): Path to the extinction coefficients CSV file.
    """
    global EXTINCTION_BASES, EXTINCTION_PAIRS

    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            seq_type = row.get("Type")
            sequence = row.get("Sequence")
            value = row.get("Value")
            
            if seq_type and sequence and value:
                if seq_type == "Base":
                    EXTINCTION_BASES[sequence.upper()] = int(value)
                elif seq_type == "Pair":
                    EXTINCTION_PAIRS[sequence.upper()] = int(value)


def calculate_ext(seq):
    """
    Calculate the extinction coefficient of a DNA sequence using nearest-neighbor and terminal contributions.

    Args:
        seq (str): DNA sequence.

    Returns:
        int: Total extinction coefficient.
    """
    if not EXTINCTION_BASES or not EXTINCTION_PAIRS:
        load_extinction_coeffs()

    seq = seq.upper()
    total = 0

    # Sum extinction values for adjacent base pairs
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        total += EXTINCTION_PAIRS.get(pair, 0)

    # Add contributions from the first and last single bases
    if seq:
        total += EXTINCTION_BASES.get(seq[0], 0)
        total += EXTINCTION_BASES.get(seq[-1], 0)

    return total

