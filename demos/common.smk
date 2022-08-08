from pathlib import Path

# Directory / file constants
SRC_DIR = Path("src")
DATA_DIR = Path("data")
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

# Helper functions
def str2bool(v):
  return v is not None and v.lower() in ("yes", "true", "t", "1")