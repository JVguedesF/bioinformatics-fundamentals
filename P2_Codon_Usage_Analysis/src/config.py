import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict


ROOT_REPO = Path(__file__).resolve().parents[2]
if str(ROOT_REPO) not in sys.path:
    sys.path.append(str(ROOT_REPO))

from utils.base_config import BaseAppConfig

CURRENT_PROJECT_DIR = Path(__file__).resolve().parent.parent

@dataclass
class AppConfig(BaseAppConfig):
    MIN_ORF_LENGTH: int = 100

    DATASETS: List[Dict[str, str]] = field(default_factory=lambda: [
        {"id": "NC_001416.1", "name": "lambda_phage.fasta", "type": "fasta"},
        {"id": "NC_000913.3", "name": "e_coli_k12.fasta", "type": "fasta"},
        {"id": "NM_000207.3", "name": "human_insulin.fasta", "type": "fasta"}
    ])

settings = AppConfig(PROJECT_ROOT=CURRENT_PROJECT_DIR)