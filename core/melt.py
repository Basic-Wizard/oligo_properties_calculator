import csv
import math

# Global dictionaries
NN_DELTA_H = {}
NN_DELTA_S = {}

# Load nearest neighbor parameters
def load_melt_values(path="data/melting_values.csv"):
    global NN_DELTA_H, NN_DELTA_S
    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            interaction = row["Interaction"].strip()
            dna_delta_h = float(row["dna_delta_h"])
            dna_delta_s = float(row["dna_delta_s"])
            NN_DELTA_H[interaction] = dna_delta_h
            NN_DELTA_S[interaction] = dna_delta_s

# Calculate Tm for a given DNA sequence
def calculate_tm(seq, concentration=500e-9):
    if not NN_DELTA_H or not NN_DELTA_S:
        load_melt_values()

    seq = seq.upper()

    total_h = 0  # Total enthalpy (kcal/mol)
    total_s = 0  # Total entropy (cal/mol/K)

    # Add initiation parameters manually
    total_h += NN_DELTA_H.get("Initiation", 0)
    total_s += NN_DELTA_S.get("Initiation", 0)

    # Sum for each neighboring pair
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        pair_key = f"{pair[0]}{pair[1]}/{complement(pair[1])}{complement(pair[0])}"
        if pair_key not in NN_DELTA_H:
            pair_key = f"{pair[1]}{pair[0]}/{complement(pair[0])}{complement(pair[1])}"
        
        total_h += NN_DELTA_H.get(pair_key, 0)
        total_s += NN_DELTA_S.get(pair_key, 0)

    # Apply nearest neighbor formula for Tm
    R = 1.987  # cal/(K*mol)
    tm_kelvin = (total_h * 1000) / (total_s + (R * math.log(concentration)))
    tm_celsius = tm_kelvin - 273.15

    return tm_celsius

# Helper to get the complement base
def complement(base):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return comp.get(base, base)

if __name__ == "__main__":
    test_seq = "ACGTACGT"
    tm = calculate_tm(test_seq)
    print(f"Melting temperature (Tm) of {test_seq}: {tm:.2f} °C")