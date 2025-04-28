import csv  # For reading nearest neighbor melting parameters
import math  # For temperature calculations

# Global dictionaries to store nearest-neighbor enthalpy and entropy values
NN_DELTA_H = {}
NN_DELTA_S = {}


def load_melt_values(path="data/melting_values.csv"):
    """
    Load nearest-neighbor enthalpy (Delta H) and entropy (Delta S) values from a CSV file.

    Args:
        path (str): Path to the melting values CSV file.
    """
    global NN_DELTA_H, NN_DELTA_S

    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interaction = row.get("Interaction", "").strip()
            dna_delta_h = row.get("dna_delta_h", "")
            dna_delta_s = row.get("dna_delta_s", "")

            if interaction and dna_delta_h and dna_delta_s:
                NN_DELTA_H[interaction] = float(dna_delta_h)
                NN_DELTA_S[interaction] = float(dna_delta_s)



def calculate_tm(seq, concentration=500e-9):
    """
    Calculate the melting temperature (Tm) of a DNA sequence using the nearest-neighbor method.

    Args:
        seq (str): DNA sequence.
        concentration (float): Strand concentration in mol/L (default 500 nM).

    Returns:
        float: Melting temperature (Tm) in degrees Celsius.
    """
    if not NN_DELTA_H or not NN_DELTA_S:
        load_melt_values()

    seq = seq.upper()

    total_h = 0  # Total enthalpy (kcal/mol)
    total_s = 0  # Total entropy (cal/mol/K)

    # Add initiation contributions if available
    total_h += NN_DELTA_H.get("Initiation", 0)
    total_s += NN_DELTA_S.get("Initiation", 0)

    # Sum contributions from each neighboring pair
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        pair_key = f"{pair[0]}{pair[1]}/{complement(pair[0])}{complement(pair[1])}"
    
        # Try direct key, else reverse
        if pair_key not in NN_DELTA_H:
            pair_key = f"{pair[1]}{pair[0]}/{complement(pair[1])}{complement(pair[0])}"
           
        total_h += NN_DELTA_H.get(pair_key, 0)
        total_s += NN_DELTA_S.get(pair_key, 0)

    # Apply the nearest-neighbor formula
    R = 1.987  # Gas constant in cal/(K*mol)
    effective_concentration = concentration / 4  # Duplex correction
    denominator = total_s + (R * math.log(effective_concentration))
    numerator = (total_h * 1000)

    tm_kelvin = numerator / denominator
    tm_celsius = tm_kelvin - 273.15


    return tm_celsius



def complement(base):
    """
    Return the complement base for DNA (A<->T, G<->C).

    Args:
        base (str): A single DNA base.

    Returns:
        str: Complementary base.
    """
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return comp.get(base.upper(), base)
