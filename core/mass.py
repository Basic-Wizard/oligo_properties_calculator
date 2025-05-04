import csv  # For reading base mass values from a CSV file

# Global dictionary to store monoisotopic mass of each nucleotide base
BASE_MASSES = {}


def load_base_masses(path="data/base_masses.csv"):
    """
    Load base monoisotopic mass values from a CSV file into the global BASE_MASSES dictionary.

    Args:
        path (str): Path to the base mass CSV file.
    """
    global BASE_MASSES
    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            base = row.get("Base")
            mass = row.get("Mass")
            if base and mass:
                BASE_MASSES[base.upper()] = float(mass)


def calculate_momo_iso_mass(seq):
    """
    Calculate the monoisotopic mass of a DNA sequence.

    Args:
        seq (str): DNA sequence composed of bases A, C, G, T.

    Returns:
        float: Total monoisotopic mass of the sequence.
    """
    if not BASE_MASSES:
        load_base_masses()

    # Sum the mass of each base; unknown bases default to 0 mass
    total_mass = sum(BASE_MASSES.get(base.upper(), 0) for base in seq)
    phosphates = (len(seq) - 1) * 61.9558
    total_mass += phosphates
    return total_mass
