import csv  # Module for reading/writing CSV files
import os    # Module for filesystem path management


def save_to_csv(output_dir, seq, results, index):
    """
    Save results for a single sequence to a CSV file.

    Args:
        output_dir (str): Directory where output file will be saved.
        seq (str): Original input sequence.
        results (dict): Calculated properties (mass, extinction, GC%, etc.).
        index (str): Sequence identifier for naming the output file.
    """
    # Create a safe filename based on sequence index
    filename = f"sequence_{index}.csv"
    filepath = os.path.join(output_dir, filename)

    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Property", "Value"])  # Write header
        writer.writerow(["Sequence", seq])       # Save original sequence
        for key, value in results.items():
            writer.writerow([key, value])         # Save result entries


def load_sequences_from_csv(filepath):
    """
    Load sequences from a CSV file with 'seq_id' and 'sequence' columns.
    If 'seq_id' is missing, auto-generate a label.
    If 'sequence' is missing, fill with a placeholder.

    Args:
        filepath (str): Path to the input CSV file.

    Returns:
        list of (seq_id, sequence) tuples.
    """
    sequences = []
    unlabeled_counter = 1

    with open(filepath, newline="") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]  # Normalize headers

        for row in reader:
            seq_id = row.get("seq_id")
            sequence = row.get("sequence")

            # Handle missing seq_id
            if not seq_id or seq_id.strip() == "":
                seq_id = f"unlabeled_sequence_{unlabeled_counter}"
                unlabeled_counter += 1
            else:
                seq_id = seq_id.strip()

            # Handle missing sequence
            if not sequence or sequence.strip() == "":
                sequence = "NO_INPUT_SEQUENCE"
            else:
                sequence = sequence.strip()

            sequences.append((seq_id, sequence))

    return sequences


def save_batch_results_to_csv(filepath, batch_results):
    """
    Save multiple sequence analysis results into one CSV file.

    Args:
        filepath (str): Output file path.
        batch_results (list): List of dictionaries containing sequence results.
    """
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        # Write headers for all properties
        headers = [
            "Sequence", "Reverse Complement", "Chemical Formula", "Length",
            "A", "C", "G", "T", "Monoisotopic Mass",
            "GC Content", "Extinction Coefficient"
        ]
        writer.writerow(headers)

        for result in batch_results:
            row = [
                result.get("Sequence", ""),
                result.get("Reverse Complement", ""),
                result.get("Chemical Formula", ""),
                result.get("Length", 0),
                result.get("Counts", {}).get("A", 0),
                result.get("Counts", {}).get("C", 0),
                result.get("Counts", {}).get("G", 0),
                result.get("Counts", {}).get("T", 0),
                result.get("Monoisotopic Mass", 0),
                result.get("GC Content", 0),
                result.get("Extinction Coefficient", 0)
            ]
            writer.writerow(row)
