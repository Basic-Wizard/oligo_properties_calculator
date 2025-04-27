from collections import Counter  # Used for counting atoms and base occurrences
from core.modifications import MODIFICATIONS, get_modification
import re


def validate_seq(seq):
    """
    Validate that the sequence only contains valid DNA bases (case-insensitive).
    """
    return all(base in "AGTCagtc" for base in seq)


def break_seq(seq):
    """
    Break the sequence into tokens, splitting around brackets [] used for modifications.
    Filters out any empty tokens that may arise from improper bracket use.

    Example:
    AC[FAM]GT -> ['AC', 'FAM', 'GT']
    """
    tokens = [t for t in re.split(r'[\[\]]', seq) if t]
    return tokens


def calculate_gc(seq):
    """
    Calculate the GC content (percentage of G and C bases) of a DNA sequence.
    """
    gc_count = seq.upper().count("G") + seq.upper().count("C")
    return (gc_count / len(seq)) * 100 if seq else 0


def count_bases(seq):
    """
    Count the occurrences of each DNA base (A, G, C, T) in the sequence.

    Returns a dictionary with counts.
    """
    return {base: seq.upper().count(base) for base in "AGCT"}


def calculate_formula(seq):
    """
    Calculate the chemical formula for a DNA sequence based on base composition.

    Returns a Counter object representing the total atom counts.
    """
    # Base compositions (monophosphate forms)
    base_composition = {
        'A': Counter({'C': 10, 'H': 13, 'N': 5, 'O': 4, 'P': 1}),
        'T': Counter({'C': 10, 'H': 14, 'N': 2, 'O': 5, 'P': 1}),
        'C': Counter({'C': 9, 'H': 14, 'N': 3, 'O': 5, 'P': 1}),
        'G': Counter({'C': 10, 'H': 14, 'N': 5, 'O': 5, 'P': 1}),
    }

    total = Counter()  # Running total for the entire sequence
    for base in seq.upper():
        if base in base_composition:
            total += base_composition[base]
        else:
            raise ValueError(f"Invalid base '{base}' in sequence.")

    return total  # Return Counter for further manipulation


def format_formula(counter_obj):
    """
    Format a Counter object representing a chemical formula into a readable string.

    Example:
    Counter({'C': 20, 'H': 27, 'N': 7, 'O': 9, 'P': 2}) -> 'C20H27N7O9P2'
    """
    return ''.join(f"{elem}{counter_obj[elem]}" for elem in sorted(counter_obj))
