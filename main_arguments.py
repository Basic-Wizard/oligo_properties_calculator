import argparse
import re
import os
from collections import Counter

from core.sequence_utils import validate_seq, break_seq, calculate_gc, count_bases, calculate_formula, format_formula
from core.io import save_to_csv, load_sequences_from_csv
from core.extinction import calculate_ext
from core.mass import calculate_momo_iso_mass
from core.modifications import load_modifications, get_modification
from core.melt import load_melt_values, calculate_tm


def parse_args():
    """
    Parse command-line arguments for input file, output directory, and verbosity.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Oligo Properties Calculator: Process sequences from a CSV file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input CSV file.")
    parser.add_argument("-o", "--output", required=True, help="Directory to save individual CSV results.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print detailed output to terminal.")
    return parser.parse_args()


def preload_data():
    """
    Load necessary base properties and modifications into memory.
    """
    load_melt_values()
    load_modifications()


def run_csv(input_path, output_dir, verbose=False):
    """
    Main function to process sequences from an input CSV and save their properties.

    Args:
        input_path (str): Path to the input CSV file.
        output_dir (str): Directory to save output CSV files.
        verbose (bool): Whether to print results to the terminal.
    """
    sequences = load_sequences_from_csv(input_path)

    for seq_id, sequence in sequences:
        tokens = break_seq(sequence)

        total_mass = 0
        total_extinction = 0
        total_tm = 0
        total_counts = Counter({"A": 0, "C": 0, "G": 0, "T": 0})
        total_formula = Counter()

        for token in tokens:
            if validate_seq(token):
                # Normal DNA segment
                token = token.upper()
                total_mass += calculate_momo_iso_mass(token)
                total_extinction += calculate_ext(token)
                total_tm += calculate_tm(token)
                total_counts.update(count_bases(token))
                total_formula.update(calculate_formula(token))
            else:
                # Assume modification
                mod = get_modification(token)
                if mod:
                    total_mass += mod["mass"]
                    total_extinction += mod["ext"]
                    # Future enhancement: add formula if available
                else:
                    print(f"Warning: Unknown modification '{token}' in sequence '{seq_id}'.")

        # Compile final output
        results = {
            "Sequence ID": seq_id,
            "Base Counts": dict(total_counts),
            "Chemical Formula": format_formula(total_formula),
            "Extinction Coefficient": total_extinction,
            "Monoisotopic Mass": total_mass,
            "Melting Temp": total_tm
        }

        if verbose:
            print(f"\nResults for sequence {seq_id}:")
            for key, value in results.items():
                print(f"{key}: {value}")

        save_to_csv(output_dir, sequence, results, seq_id)


if __name__ == "__main__":
    args = parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    preload_data()
    run_csv(args.input, args.output, verbose=args.verbose)
