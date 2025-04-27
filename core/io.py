import csv  # CSV module for writing results to a file
import os

def save_to_csv(output_dir, seq, results, index):
    # Sanitize filename (optional, based on sequence or just use index)
    filename = f"sequence_{index}.csv"
    filepath = os.path.join(output_dir, filename)

    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Property", "Value"])  # Header
        writer.writerow(["Sequence", seq])       # First row: the actual sequence
        for key, value in results.items():
            writer.writerow([key, value])         # Then each property


def load_sequences_from_csv(filepath):
    sequences = []
    unlabeled_counter = 1

    with open(filepath, newline="") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]

        for row in reader:
            seq_id = row.get("seq_id")
            sequence = row.get("sequence")

            # Fill missing seq_id
            if not seq_id or seq_id.strip() == "":
                seq_id = f"unlabeled_sequence_{unlabeled_counter}"
                unlabeled_counter += 1
            else:
                seq_id = seq_id.strip()

            # Fill missing sequence
            if not sequence or sequence.strip() == "":
                sequence = "NO_INPUT_SEQUENCE"
            else:
                sequence = sequence.strip()

            sequences.append((seq_id, sequence))

    return sequences


def save_batch_results_to_csv(filepath, batch_results):
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        # Headers
        headers = [
            "Sequence", "Reverse Complement", "Chemical Formula", "Length",
            "A", "C", "G", "T", "Monoisotopic Mass",
            "GC Content", "Extinction Coefficient"
        ]
        writer.writerow(headers)
        for result in batch_results:
            row = [
                result["Sequence"],
                result["Reverse Complement"],
                result["Chemical Formula"],
                result["Length"],
                result["Counts"].get("A", 0),
                result["Counts"].get("C", 0),
                result["Counts"].get("G", 0),
                result["Counts"].get("T", 0),
                result["Monoisotopic Mass"],
                result["GC Content"],
                result["Extinction Coefficient"]
            ]
            writer.writerow(row)