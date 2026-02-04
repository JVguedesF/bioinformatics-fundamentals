from pathlib import Path
import os
from dataclasses import dataclass, field
from typing import List, Dict
from dotenv import load_dotenv

load_dotenv()


@dataclass
class AppConfig:
    BASE_PATH: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = BASE_PATH / "data" / "sequences"
    RESULTS_DIR: Path = BASE_PATH / "results"

    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: ["*.fasta", "*.fa"])

    MIN_ORF_LENGTH: int = 100

    ENTREZ_EMAIL: str = field(default_factory=lambda: os.getenv("ENTREZ_EMAIL"))

    DATASETS: List[Dict[str, str]] = field(default_factory=lambda: [
        {"id": "NC_001416.1", "name": "lambda_phage.fasta", "type": "fasta"},
        {"id": "NC_000913.3", "name": "e_coli_k12.fasta", "type": "fasta"},
        {"id": "NM_000207.3", "name": "human_insulin.fasta", "type": "fasta"}
    ])


settings = AppConfig()