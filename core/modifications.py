import csv

MODIFICATIONS = {}

def load_modifications(path="data/modifications.csv"):
    global MODIFICATIONS
    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                MODIFICATIONS[row["Synoligo Code"]] = {
                    "name": row["Description"],
                    "mass": float(row["Formula Weight"]),
                    "ext": int(row["Extinction Coefficient"]),
                    "positions": {
                        "5'": row["5′ Mod"].strip().lower() == "yes",
                        "internal": row["Internal Mod"].strip().lower() == "yes",
                        "3'": row["3′ Mod"].strip().lower() == "yes"
                    }
                }
            except ValueError as e:
                print(f"Error parsing row {row}: {e}")

def get_modification(symbol):
    return MODIFICATIONS.get(symbol)