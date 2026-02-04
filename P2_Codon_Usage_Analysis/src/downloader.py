from Bio import Entrez
from pathlib import Path
from rich.console import Console
from config import settings
import sys

Entrez.email = settings.ENTREZ_EMAIL
console = Console()

def download_datasets(output_dir: Path):
    if not Entrez.email:
        console.print("[red]Erro: ENTREZ_EMAIL não configurado no arquivo .env![/]")
        console.print("[yellow]Crie um arquivo .env na raiz do projeto com: ENTREZ_EMAIL=seu@email.com[/]")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    for data in settings.DATASETS:
        file_path = output_dir / data["name"]

        if not file_path.exists():
            console.print(f"[yellow]Downloading {data['name']} ({data['id']})...[/]")
            try:
                with Entrez.efetch(db="nucleotide", id=data["id"], rettype=data["type"], retmode="text") as handle:
                    seq_data = handle.read()
                    with open(file_path, "w") as f:
                        f.write(seq_data)
                console.print(f"[green]Successfully saved {data['name']}[/]")
            except Exception as e:
                console.print(f"[red]Failed to download {data['id']}: {e}[/]")