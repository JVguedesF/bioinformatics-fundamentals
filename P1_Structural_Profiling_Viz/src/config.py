from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List

@dataclass
class AppConfig:
    """Centraliza configurações da aplicação."""

    BASE_PATH: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = BASE_PATH / "data" / "sequences"
    RESULTS_DIR: Path = BASE_PATH / "results"

    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: ["*.fasta", "*.fa"])

    MIN_PROTEIN_LENGTH: int = 5
    PREVIEW_LENGTH: int = 60

    TYPE_COLORS: Dict[str, str] = field(default_factory=lambda: {
        'DNA': 'blue',
        'RNA': 'magenta',
        'PROTEIN': 'green',
        'UNKNOWN': 'red'
    })

settings = AppConfig()