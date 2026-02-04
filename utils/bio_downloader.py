from Bio import Entrez
from pathlib import Path
from rich.console import Console
from typing import List, Dict
import sys

console = Console()


def download_entrez_data(email: str, datasets: List[Dict], output_dir: Path):
    """
    Função genérica para baixar dados do NCBI/Entrez.
    Recebe as configurações como argumentos, desacoplando do arquivo settings global.
    """
    if not email:
        console.print("[red]Erro: Email do Entrez não fornecido![/]")
        console.print("[yellow]Verifique seu arquivo .env e as configurações.[/]")
        sys.exit(1)

    Entrez.email = email
    output_dir.mkdir(parents=True, exist_ok=True)

    for data in datasets:
        file_path = output_dir / data["name"]
        if file_path.exists():
            console.print(f"[cyan]Arquivo já existe: {data['name']} (pulando download)[/]")
            continue

        console.print(f"[yellow]Downloading {data['name']} ({data['id']})...[/]")

        try:
            r_type = data.get("type", "fasta")
            r_mode = "text"

            with Entrez.efetch(db="nucleotide", id=data["id"], rettype=r_type, retmode=r_mode) as handle:
                seq_data = handle.read()
                with open(file_path, "w") as f:
                    f.write(seq_data)

            console.print(f"[green]Successfully saved {data['name']}[/]")

        except Exception as e:
            console.print(f"[red]Failed to download {data['id']}: {e}[/]")