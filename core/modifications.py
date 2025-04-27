import csv

MODIFICATIONS = {}


def normalize_token(token):
    # Replace curly quotes with straight quotes
    token = token.strip().upper()
    token = token.replace("’", "'").replace("′", "'")
    token = token.replace("‘", "'").replace("`", "'")
    return token

def load_modifications(path="data/modifications.csv"):
    global MODIFICATIONS
    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                code = normalize_token(row["Synoligo Code"])
                MODIFICATIONS[code] = {
                    "name": row["Description"],
                    "mass": float(row["Formula Weight"]),
                    "ext": int(row["Extinction Coefficient"]),
                }
            except ValueError as e:
                print(f"Error parsing row {row}: {e}")


def get_modification(symbol):
    symbol = normalize_token(symbol)
    return MODIFICATIONS.get(symbol)