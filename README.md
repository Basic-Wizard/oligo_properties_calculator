# Oligonucleotide Properties Calculator

A Python tool to process DNA sequences and compute important properties such as:

- Base counts
- Monoisotopic mass
- Extinction coefficient
- Melting temperature (Tm)
- Chemical formula

Each sequence is analyzed and saved individually as a CSV file.

---

## Project Structure

```
core/
  extinction.py
  io.py
  mass.py
  melt.py
  modifications.py
  sequence_utils.py

main_arguments.py
README.md
```

---

## Installation

1. **Clone the repository:**

```bash
git clone https://github.com/your-username/oligo_properties_calculator.git
cd oligo_properties_calculator
```

2. **Install required dependencies:**

This tool only uses Python's standard library (no external packages needed).

Ensure you are using Python 3.7+.

---

## Usage

```bash
python main_arguments.py -i <input_csv> -o <output_directory> [-v]
```

### Arguments:

- `-i`, `--input` **(required)**: Path to the input CSV file containing sequences.
- `-o`, `--output` **(required)**: Path to the output directory where results will be saved.
- `-v`, `--verbose` **(optional)**: Print detailed results to the terminal.

### Example:

```bash
python main_arguments.py -i data/test_sequences.csv -o results/ -v
```

If the output directory does not exist, it will be created automatically.

---

## Input CSV Format

The input CSV must have at least the following columns:

| seq_id    | sequence  |
| --------- | --------- |
| Example_1 | ACGTACGT  |
| Example_2 | AC[FAM]GT |

- `seq_id`: A unique label for each sequence (optional, autogenerated if missing)
- `sequence`: The DNA sequence (modifications should be enclosed in square brackets, e.g., `[FAM]`)

If a sequence ID or sequence is missing, the tool will auto-generate a placeholder.

---

## Output

Each sequence will produce its own CSV file under the specified output directory.

Each result file contains:

| Property               | Value                      |
| ---------------------- | -------------------------- |
| Sequence               | Original sequence          |
| Reverse Compliment     | Reverse Compliment         |
| Sequence ID            | ID from input or generated |
| Base Counts            | Number of A, C, G, T bases |
| Chemical Formula       | Aggregate chemical formula |
| Extinction Coefficient | Absorbance measure         |
| Monoisotopic Mass      | Total mass in Daltons      |
| Melting Temp           | Calculated Tm (Celsius)    |

---

## Supported Modifications

The tool can recognize chemical modifications based on the `data/modifications.csv` file.

- Unknown modifications will be reported as warnings.

To add new supported modifications:

- Edit the `modifications.csv` and reload.

---

## Notes

- Only DNA bases (A, C, G, T) and bracketed modifications are supported.
- All calculations are based on common nearest-neighbor thermodynamic models and extinction coefficients.

---

## License

MIT License. See `LICENSE` file for details.

---

## Future Improvements

- Support for RNA sequences and mixed chirality
- Detailed chemical formula parsing for modifications
- Batch summary CSV output

---

## Contact

For questions or contributions, please contact:

**Your Name**  
[GitHub Profile](https://github.com/your-username)

---

Happy calculating! 🚀
