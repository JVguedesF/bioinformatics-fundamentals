import sys
import warnings
from pathlib import Path
from rich.console import Console
from Bio import BiopythonWarning

def init_pipeline(script_file: str) -> Console:
    warnings.simplefilter('ignore', BiopythonWarning)
    
    src_path = str(Path(script_file).resolve().parent)
    if src_path not in sys.path:
        sys.path.append(src_path)
        
    return Console()