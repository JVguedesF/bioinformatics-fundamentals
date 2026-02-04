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
        {"id": "NG_017013.2", "name": "tp53_human_genomic_refseqgene.fasta", "type": "fasta"},
        {"id": "NM_000546.6", "name": "tp53_human_transcript_variant1.fasta", "type": "fasta"},
        {"id": "NM_001126112.2", "name": "tp53_human_transcript_variant2.fasta", "type": "fasta"},
        {"id": "NG_007557.1", "name": "trp53_mouse_genomic_ortholog.fasta", "type": "fasta"}
    ])

settings = AppConfig(PROJECT_ROOT=CURRENT_PROJECT_DIR)