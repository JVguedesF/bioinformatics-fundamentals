from pathlib import Path
import os
from dataclasses import dataclass, field
from typing import List, Dict
from dotenv import load_dotenv

load_dotenv()


@dataclass
class BaseAppConfig:
    """
    Configuração Base compartilhada por todos os projetos.
    Define caminhos padrões e carrega variáveis de ambiente.
    """
    PROJECT_ROOT: Path

    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: ["*.fasta", "*.fa"])
    ENTREZ_EMAIL: str = field(default_factory=lambda: os.environ.get("ENTREZ_EMAIL"))

    DATASETS: List[Dict[str, str]] = field(default_factory=list)

    def __post_init__(self):
        self.DATA_DIR = self.PROJECT_ROOT / "data" / "sequences"
        self.RESULTS_DIR = self.PROJECT_ROOT / "results"

        self.DATA_DIR.mkdir(parents=True, exist_ok=True)
        self.RESULTS_DIR.mkdir(parents=True, exist_ok=True)