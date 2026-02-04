import sys
import argparse
from pathlib import Path
from datetime import datetime

from Bio import SeqIO
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from utils.startup import init_pipeline

console = init_pipeline(__file__)

from analyzer import SequenceAnalyzer
from downloader import download_datasets
from view import CodonView
from config import settings
from utils.exceptions import BioPipelineError


def parse_arguments():
    parser = argparse.ArgumentParser(description="GeneCodePro v1.0 - Genetic Code Analysis")
    parser.add_argument("--input", "-i", type=Path, help="Diretório de entrada (arquivos FASTA)")
    parser.add_argument("--output", "-o", type=Path, help="Diretório de saída para relatórios")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.input:
        settings.DATA_DIR = args.input
    if args.output:
        settings.RESULTS_DIR = args.output

    console.print()
    console.rule("[bold magenta]GeneCodePro v1.0 - Genetic Code Analysis[/bold magenta]")
    console.print(f"[dim]Started: {datetime.now().isoformat()}[/dim]")
    console.print()

    with console.status("[bold green]Syncing datasets...[/]"):
        download_datasets(settings.DATA_DIR)

    files = []
    for ext in settings.FILE_EXTENSIONS:
        files.extend(settings.DATA_DIR.glob(ext))

    pipeline_results = []

    with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
    ) as progress:

        task = progress.add_task("Analyzing Genomes...", total=len(files))

        for fasta_path in files:
            progress.update(task, description=f"Processing {fasta_path.name}")

            dataset_info = next((d for d in settings.DATASETS if d["name"] == fasta_path.name), None)
            table_id = dataset_info.get("table", 1) if dataset_info else 1

            try:
                record = SeqIO.read(fasta_path, "fasta")
                analyzer = SequenceAnalyzer(record.id, record.seq, table_id)

                orfs = analyzer.find_orfs(min_len_aa=settings.MIN_ORF_LENGTH)
                metrics = analyzer.analyze_codon_usage()

                header, tree = CodonView.create_analysis_view(
                    filename=fasta_path.name,
                    record_id=record.id,
                    length=len(record.seq),
                    table_id=table_id,
                    orfs=orfs,
                    metrics=metrics
                )

                console.print(header)
                console.print(tree)
                console.print()

                gc_pct = (metrics.gc_content / (metrics.total_codons * 3)) * 100 if metrics.total_codons > 0 else 0
                pipeline_results.append({
                    "id": record.id,
                    "file": fasta_path.name,
                    "table": table_id,
                    "length": len(record.seq),
                    "orf_count": len(orfs),
                    "gc_percent": gc_pct,
                    "top_codons": metrics.top_codons
                })

            except BioPipelineError as err:
                console.print(f"[orange1]Skipping {fasta_path.name}: {err}[/]")
            except Exception as err:
                console.print(f"[red]Unexpected Error in {fasta_path.name}: {err}[/]")

            progress.advance(task)

    console.rule("[bold]Generating Reports[/bold]")
    if pipeline_results:
        CodonView.export_results(pipeline_results, settings.RESULTS_DIR)

    console.print()
    console.rule("[bold green]Pipeline Complete[/bold green]")


if __name__ == "__main__":
    main()