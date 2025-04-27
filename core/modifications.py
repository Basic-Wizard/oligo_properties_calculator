import csv  # For reading modification data from a CSV file

# Global dictionary to store modification data
MODIFICATIONS = {}


def normalize_token(token):
    """
    Normalize a modification token by stripping whitespace,
    converting to uppercase, and fixing curly/special quotes.

    Args:
        token (str): Raw token string.

    Returns:
        str: Normalized token.
    """
    token = token.strip().upper()
    token = token.replace("’", "'").replace("′", "'")
    token = token.replace("‘", "'").replace("`", "'")
    return token


def load_modifications(path="data/modifications.csv"):
    """
    Load chemical modifications from a CSV file into the global MODIFICATIONS dictionary.

    Args:
        path (str): Path to the modifications CSV file.
    """
    global MODIFICATIONS

    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                code = normalize_token(row.get("Synoligo Code", ""))
                description = row.get("Description", "")
                mass = row.get("Formula Weight", "")
                ext = row.get("Extinction Coefficient", "")

                if code and mass and ext:
                    MODIFICATIONS[code] = {
                        "name": description,
                        "mass": float(mass),
                        "ext": int(ext)
                    }
            except ValueError as e:
                print(f"Error parsing row {row}: {e}")


def get_modification(symbol):
    """
    Retrieve the modification details for a given symbol.

    Args:
        symbol (str): Modification token.

    Returns:
        dict or None: Modification details (mass, extinction) if found.
    """
    symbol = normalize_token(symbol)
    return MODIFICATIONS.get(symbol)