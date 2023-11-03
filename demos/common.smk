from pathlib import Path
import platform

# Check if this is running on O2
IS_O2 = (platform.system() == "Linux")

# Directory / file constants
SRC_DIR = Path("src")
DATA_DIR = Path("data" if not IS_O2 else "/n/data1/hms/dbmi/gehlenborg/lab/vitessce-python-demos")
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

# Helper functions
def str2bool(v):
  return v is not None and v.lower() in ("yes", "true", "t", "1")

def flatten(l):
  return [item for sublist in l for item in sublist]

def is_aws(output_path):
  return not output_path.endswith('.ome.zarr')