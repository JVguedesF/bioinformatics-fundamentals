import sys
from pathlib import Path

ROOT_REPO = Path(__file__).resolve().parents[2]
if str(ROOT_REPO) not in sys.path:
    sys.path.append(str(ROOT_REPO))

from config import settings
from utils.bio_downloader import download_entrez_data

def download_datasets(output_dir: Path):
    download_entrez_data(
        email=settings.ENTREZ_EMAIL,
        datasets=settings.DATASETS,
        output_dir=output_dir
    )