import sys
import argparse
from pathlib import Path
from datetime import datetime

from Bio import SeqIO
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from utils.startup import init_pipeline

console = init_pipeline(__file__)

from analyzer import CentralDogmaAnalyzer
from view import StructuralView
from config import settings
from utils.exceptions import BioPipelineError

def parse_arguments():
    """Configura e processa argumentos da linha de comando."""
    parser = argparse.ArgumentParser(description="BioPipeline v4.0 - Molecular Analysis Pipeline")
    parser.add_argument("--input", "-i", type=Path, help="Diretório de entrada (arquivos FASTA)")
    parser.add_argument("--output", "-o", type=Path, help="Diretório de saída para relatórios")
    return parser.parse_args()

def main():
    """
    Main pipeline execution.
    """
    args = parse_arguments()
    if args.input:
        settings.DATA_DIR = args.input
    if args.output:
        settings.RESULTS_DIR = args.output

    if not settings.DATA_DIR.exists():
        console.print(f"[red]Error: Directory not found: {settings.DATA_DIR}[/]")
        return

    fasta_files = []
    for ext in settings.FILE_EXTENSIONS:
        fasta_files.extend(settings.DATA_DIR.glob(ext))
    fasta_files = sorted(fasta_files)

    if not fasta_files:
        console.print("[yellow]No FASTA files found in data/sequences[/]")
        return

    console.print()
    console.rule("[bold]BioPipeline v4.0 - Molecular Analysis Pipeline[/bold]")
    console.print(f"[dim]Analysis started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}[/dim]")
    console.print()

    results = []
    analyzer = CentralDogmaAnalyzer()

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console
    ) as progress:

        task = progress.add_task("Processing sequences...", total=len(fasta_files))

        for i, fasta_path in enumerate(fasta_files, 1):
            progress.update(task, description=f"Analyzing {fasta_path.name}")

            try:
                record = SeqIO.read(fasta_path, "fasta")
                result = analyzer.process_sequence(record, fasta_path.name)
                # Usa a View para criar a visualização
                header, tree = StructuralView.create_visualization(result, str(record.seq))
                
                results.append(result)

                console.rule(f"[bold]Sequence {i}/{len(fasta_files)}[/bold]")
                console.print(header)
                console.print(tree)
                console.print()
            except BioPipelineError as err:
                console.print(f"[orange1]Analysis Skipped for {fasta_path.name}: {err}[/]")
            except Exception as err:
                console.print(f"[red]Unexpected Error in {fasta_path.name}: {err}[/]")

            progress.advance(task)

    console.print()
    if results:
        console.print(StructuralView.create_summary_table(results))
        console.print()
        console.print(StructuralView.generate_statistics_panel(results))
        StructuralView.export_results(results, settings.RESULTS_DIR)

    console.print()
    console.rule("[bold green]Analysis Complete[/bold green]")


if __name__ == "__main__":
    main()