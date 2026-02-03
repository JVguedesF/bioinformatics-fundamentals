from typing import List, Tuple
from pathlib import Path
from dataclasses import asdict
from rich.panel import Panel
from rich.tree import Tree
from rich.text import Text
from rich.table import Table
from rich.console import Console

from analyzer import AnalysisResult
from utils.visualizer import BioVisualizer
from utils.reporter import BioReporter
from config import settings

console = Console()

class StructuralView:
    """Gerencia toda a visualização e relatórios do Projeto 1."""

    @staticmethod
    def create_visualization(result: AnalysisResult, sequence_preview: str) -> Tuple[Panel, Tree]:
        color = settings.TYPE_COLORS.get(result.molecule_type, 'white')

        header = BioVisualizer.create_header_panel(
            result.filename, result.sequence_id, result.molecule_type, result.length, color
        )

        branches = []
        seq_text = Text("Sequence: ").append(BioVisualizer.format_sequence_preview(sequence_preview))
        branches.append((None, [seq_text]))

        if result.gc_content is not None:
            items = [
                f"GC Content: {result.gc_content:.2f}%",
                f"Melting Temperature: {result.melting_temp:.2f}°C"
            ]
            branches.append(("[blue]Genomic Analysis[/blue]", items))

        if result.mfe is not None:
            items = [f"Minimum Free Energy: {result.mfe:.2f} kcal/mol"]
            if result.secondary_structure:
                items.append(f"Secondary Structure: {result.secondary_structure[:40]}...")
            branches.append(("[magenta]Transcriptomic Analysis[/magenta]", items))

        if result.molecular_weight is not None:
            items = [
                f"Molecular Weight: {result.molecular_weight:.2f} Da",
                f"Isoelectric Point: {result.isoelectric_point:.2f}"
            ]
            if result.stability_index:
                items.append(f"Instability Index: {result.stability_index}")
            branches.append(("[green]Proteomic Analysis[/green]", items))

        tree = BioVisualizer.create_result_tree("Central Dogma Analysis", branches)
        return header, tree

    @staticmethod
    def create_summary_table(results: List[AnalysisResult]) -> Table:
        table = Table(title="Comparative Analysis Summary", border_style="bright_black", header_style="bold cyan", show_lines=True)
        table.add_column("File", style="white", no_wrap=True)
        table.add_column("Type", justify="center")
        table.add_column("Length", justify="right")
        table.add_column("GC% / pI", justify="right")
        table.add_column("Tm / MFE", justify="right")

        for res in results:
            if res.molecule_type == 'DNA':
                m1, m2, style = f"{res.gc_content:.1f}%", f"{res.melting_temp:.1f}°C", "blue"
            elif res.molecule_type == 'RNA':
                m1, m2, style = "N/A", f"{res.mfe:.1f} kcal/mol", "magenta"
            else:
                m1, m2, style = f"{res.isoelectric_point:.2f}", f"{res.molecular_weight:.0f} Da", "green"

            table.add_row(res.filename, f"[{style}]{res.molecule_type}[/{style}]", str(res.length), m1, m2)
        return table

    @staticmethod
    def generate_statistics_panel(results: List[AnalysisResult]) -> Panel:
        total = len(results)
        counts = {t: sum(1 for r in results if r.molecule_type == t) for t in ['DNA', 'RNA', 'PROTEIN']}
        avg_len = sum(r.length for r in results) / total if total > 0 else 0
        
        text = (f"[bold]Dataset Overview[/bold]\n\nTotal: {total}\n"
                f"DNA: {counts['DNA']}\nRNA: {counts['RNA']}\nProtein: {counts['PROTEIN']}\n"
                f"Avg Length: {avg_len:.0f}")
        return Panel(text, title="Summary Statistics", border_style="cyan")

    @staticmethod
    def export_results(results: List[AnalysisResult], output_dir: Path):
        output_dir.mkdir(exist_ok=True)
        counts = {t: sum(1 for r in results if r.molecule_type == t) for t in ['DNA', 'RNA', 'PROTEIN']}
        summary = {"Total Sequences": len(results), **{f"{k} Count": v for k, v in counts.items()}}
        
        headers = ["File", "ID", "Type", "Length", "GC% / pI", "Tm / MFE"]
        rows = []
        for res in results:
            if res.molecule_type == 'DNA':
                m1, m2 = f"{res.gc_content:.2f}", f"{res.melting_temp:.2f}"
            elif res.molecule_type == 'RNA':
                m1, m2 = "-", f"{res.mfe:.2f}"
            else:
                m1, m2 = f"{res.isoelectric_point:.2f}", f"{res.molecular_weight:.0f}"
            rows.append([res.filename, res.sequence_id, res.molecule_type, str(res.length), m1, m2])

        console.print(f"\n[green]Results exported to:[/green]")
        data_dicts = [asdict(r) for r in results]
        BioReporter.export_json(data_dicts, output_dir)
        BioReporter.export_csv(data_dicts, output_dir)
        tex = BioReporter.generate_latex("Structural Profiling Report", summary, headers, rows, output_dir)
        console.print(f"  LaTeX: {tex.name}")