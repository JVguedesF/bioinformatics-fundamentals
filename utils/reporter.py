from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any
import json
import csv
from rich.console import Console

console = Console()

class BioReporter:
    """
    Gerador de relatórios genérico para todos os projetos de bioinformática.
    Gera PDF (LaTeX) e JSON padronizados.
    """
    
    @staticmethod
    def generate_latex(
        title: str,
        summary_data: Dict[str, Any],
        headers: List[str],
        rows: List[List[str]],
        output_dir: Path,
        filename_prefix: str = "report"
    ) -> Path:
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        tex_path = output_dir / f"{filename_prefix}_{timestamp}.tex"

        summary_text = "\n".join([f"\\item \\textbf{{{k}}}: {v}" for k, v in summary_data.items()])

        col_format = "l" * len(headers) # Alinha tudo à esquerda
        header_row = " & ".join([f"\\textbf{{{h}}}" for h in headers])

        table_content = ""
        for row in rows:
            safe_row = [str(item).replace('_', r'\_').replace('%', r'\%') for item in row]
            table_content += " & ".join(safe_row) + r" \\" + "\n"

        latex_template = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{geometry}
\geometry{margin=2.5cm}

\title{""" + title + r"""}
\author{Bioinformatica Iniciante Pipeline}
\date{\today}

\begin{document}
\maketitle

\section{Summary}
\begin{itemize}
""" + summary_text + r"""
\end{itemize}

\section{Detailed Results}
\begin{longtable}{""" + col_format + r"""}
\toprule
""" + header_row + r""" \\
\midrule
""" + table_content + r"""
\bottomrule
\end{longtable}
\end{document}
"""
        
        with open(tex_path, 'w', encoding='utf-8') as f:
            f.write(latex_template)
            
        return tex_path

    @staticmethod
    def export_json(data: List[Dict], output_dir: Path, prefix: str = "data"):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        json_path = output_dir / f"{prefix}_{timestamp}.json"
        
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
            
        console.print(f"  JSON:  {json_path.name}")
        return json_path

    @staticmethod
    def export_csv(data: List[Dict], output_dir: Path, prefix: str = "data"):
        if not data:
            return None

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_path = output_dir / f"{prefix}_{timestamp}.csv"

        fieldnames = list(data[0].keys())
        
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
            
        console.print(f"  CSV:   {csv_path.name}")
        return csv_path