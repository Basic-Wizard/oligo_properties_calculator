from collections import Counter  # Used for counting atoms and base occurrences
from core.modifications import MODIFICATIONS, get_modification
import re

def validate_seq(seq):
    return all(base in "AGTC" for base in seq)

    # tokens = tokenize_sequence(seq, MODIFICATIONS)
    # for token in tokens:
    #     if token not in "AGTC" and not get_modification(token):
    #         print(f"Unknown token in sequence: {token}")
    #         return False
    # return True

def break_seq(seq):
    tokens = re.split('[\[\]]',seq)
    return tokens

def calculate_gc(seq):
    gc_count = seq.count("G") + seq.count("C")  # Count G and C bases
    return (gc_count / len(seq)) * 100 if seq else 0  # Calculate GC% only if seq is non-empty

def count_bases(seq):
    # Count each base in "AGCT" and return as dictionary
    return {base: seq.count(base) for base in "AGCT"}

# def rev_comp(seq):
#     complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  # DNA base complements
#     comp_strand = ''.join(complement[base] for base in seq)  # Get complement strand
#     return comp_strand[::-1]  # Reverse to get reverse complement

def calculate_formula(seq):
    base_composition = {
        'A': Counter({'C': 10, 'H': 13, 'N': 5, 'O': 4, 'P': 1}),
        'T': Counter({'C': 10, 'H': 14, 'N': 2, 'O': 5, 'P': 1}),
        'C': Counter({'C': 9, 'H': 14, 'N': 3, 'O': 5, 'P': 1}),
        'G': Counter({'C': 10, 'H': 14, 'N': 5, 'O': 5, 'P': 1}),
    }

    total = Counter()  # Accumulate total atomic composition
    for base in seq.upper():
        if base in base_composition:
            total += base_composition[base]  # Add base's atom counts
        else:
            raise ValueError(f"invalid base '{base}' in sequence.")  # Alert for invalid base

    # Create molecular formula string in sorted atomic order
    formula = ''.join(f"{elem}{total[elem]}" for elem in sorted(total))
    
    return formula

