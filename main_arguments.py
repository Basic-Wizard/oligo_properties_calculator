import argparse
import re
from collections import Counter

from core.sequence_utils import validate_seq,break_seq, calculate_gc,count_bases, calculate_formula
from core.io import save_to_csv,load_sequences_from_csv, save_batch_results_to_csv
from core.extinction import calculate_ext
from core.mass import calculate_momo_iso_mass
from core.modifications import load_modifications, get_modification

def parse_args():
    parser = argparse.ArgumentParser(description="Split a CSV into multiple single-row CSV files.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input CSV file.")
    parser.add_argument("-o", "--output", required=True, help="Directory to save individual CSV files.")
    return parser.parse_args()

args = parse_args()
i = args.input
o = args.output



def parse_formula(formula):
    pattern = re.compile(r"([A-Z][a-z]*)(\d*)")
    counter = Counter()
    for (element, count) in re.findall(pattern, formula):
        counter[element] += int(count) if count else 1
    return counter


def run_csv(input):
    sequences = load_sequences_from_csv(input)

    for seq in sequences:
        idx = seq[0]
        seq = seq[1]
        tokens = break_seq(seq)

        total_mass = 0
        total_extinction = 0
        total_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        total_formula = Counter()

        for token in tokens:
            if validate_seq(token):
                # Normal base sequence
                bases = count_bases(token)
                formula = calculate_formula(token)
                ext = calculate_ext(token)
                mass = calculate_momo_iso_mass(token)

                # Accumulate totals
                for base in "ACGT":
                    total_counts[base] += bases.get(base, 0)

                # Parse formula like C10H13N5O4P → Counter({'C':10,'H':13,...})
                token_formula = parse_formula(formula)
                total_formula += token_formula

                total_extinction += ext
                total_mass += mass
            else:
                # Modification
                load_modifications
                mod = get_modification(token)
                if mod:
                    print(mod["mass"],mod["ext"])
                    total_mass += mod["mass"]
                    total_extinction += mod["ext"]

                    # If you know the chemical formula for mods, you could add that too here!
                    # (Would require loading formula strings for mods, not just mass/ext)
                else:
                    print(f"Unknown modification: {token}")

        # After processing all tokens, now output the final total properties
        final_formula = ''.join(f"{elem}{total_formula[elem]}" for elem in sorted(total_formula))

        # Create result dictionary
        results = {
            "Base Counts": total_counts,
            "Chemical Formula": final_formula,
            "Extinction Coefficient": total_extinction,
            "Monoisotopic Mass": total_mass
        }

        print(f"\nResults for sequence {seq}:")
        print(results)
        save_to_csv(o, seq, results, idx)
        
        # You can call save_to_csv() here if you want to write to a file
        # save_to_csv(filepath, seq, results)

if __name__ == "__main__":
    run_csv(i)