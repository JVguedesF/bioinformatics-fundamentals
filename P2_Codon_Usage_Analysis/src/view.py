from typing import List, Tuple, Union
from pathlib import Path
from rich.panel import Panel
from rich.tree import Tree
from rich.console import Console, RenderableType

from analyzer import ORFResult, CodonMetrics
from utils.visualizer import BioVisualizer
from utils.reporter import BioReporter

console = Console()

class CodonView:
    """Gerencia visualização e relatórios do Projeto 2."""

    @staticmethod
    def create_analysis_view(filename: str, record_id: str, length: int, table_id: int,
                            orfs: List[ORFResult], metrics: CodonMetrics) -> Tuple[Panel, Tree]:
        org_type = "BACTERIA (Table 11)" if table_id == 11 else "STANDARD (Table 1)"
        color = "cyan" if table_id == 11 else "green"

        header = BioVisualizer.create_header_panel(filename, record_id, org_type, length, color)
        branches = []
        
        # Genomic Metrics
        gc_pct = (metrics.gc_content / (metrics.total_codons * 3)) * 100 if metrics.total_codons > 0 else 0
        branches.append((f"[{color}]Genomic Metrics[/{color}]", [
            f"GC Content: {gc_pct:.2f}%",
            f"Total Codons: {metrics.total_codons}"
        ]))

        # ORF Prediction
        orf_items: List[Union[str, RenderableType]] = [f"Total ORFs detected: {len(orfs)} (>100aa)"]
        if orfs:
            best = orfs[0]
            best_tree = Tree(f"[bold]Longest Candidate:[/bold] {best.length_aa} aa")
            best_tree.add(f"Frame: {best.frame} | Strand: {best.strand}")
            best_tree.add(f"Preview: [dim]{best.protein_seq[:25]}...[/dim]")
            orf_items.append(best_tree)
        branches.append(("[yellow]ORF Prediction[/yellow]", orf_items))

        # Codon Bias
        top_codons = [f"[bold]{c}[/bold]: {n} occurrences" for c, n in metrics.top_codons[:3]]
        branches.append(("[magenta]Codon Bias Profile[/magenta]", top_codons))

        tree = BioVisualizer.create_result_tree("Genetic Code & ORF Analysis", branches)
        return header, tree

    @staticmethod
    def export_results(results: List[dict], output_dir: Path):
        output_dir.mkdir(parents=True, exist_ok=True)
        total_seqs = len(results)
        avg_gc = sum(r['gc_percent'] for r in results) / total_seqs if total_seqs else 0
        
        summary = {
            "Total Genomes Analyzed": total_seqs,
            "Average GC Content": f"{avg_gc:.2f}%",
            "Total ORFs Detected": sum(r['orf_count'] for r in results)
        }
        headers = ["Organism/ID", "Table", "Length", "GC%", "ORFs", "Top Codons"]
        rows = [[r['id'], str(r['table']), str(r['length']), f"{r['gc_percent']:.2f}%", str(r['orf_count']), ", ".join([x[0] for x in r['top_codons'][:3]])] for r in results]

        console.print(f"\n[green]Reports Generated:[/green]")
        BioReporter.export_json(results, output_dir, prefix="codon_usage")
        BioReporter.export_csv(results, output_dir, prefix="codon_usage")
        tex = BioReporter.generate_latex("GeneCodePro Analysis Report", summary, headers, rows, output_dir)
        console.print(f"  LaTeX: {tex.name}")