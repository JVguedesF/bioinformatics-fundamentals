from dataclasses import dataclass
from typing import List, Dict, Tuple
from collections import Counter
from Bio.Seq import Seq
from Bio.Data import CodonTable

from utils.exceptions import SequenceAnalysisError, EmptySequenceError


@dataclass
class ORFResult:
    """Container for identified Open Reading Frame data."""
    id: str
    frame: int
    strand: int
    length_aa: int
    length_bp: int
    protein_seq: str
    start_pos: int  # Posição em nucleotídeos (bp)


@dataclass
class CodonMetrics:
    """Container for Codon Usage statistics."""
    total_codons: int
    gc_content: float
    codon_counts: Dict[str, int]
    amino_acid_counts: Dict[str, int]
    top_codons: List[Tuple[str, int]]


class SequenceAnalyzer:
    def __init__(self, sequence_id: str, sequence: Seq, table_id: int = 1):
        if not sequence:
            raise EmptySequenceError(f"Sequence {sequence_id} is empty.")

        self.sequence_id = sequence_id
        self.sequence = sequence
        self.table_id = table_id
        self.table = CodonTable.unambiguous_dna_by_id[table_id]

    def find_orfs(self, min_len_aa: int = 100) -> List[ORFResult]:
        """
        Scan sequence for ORFs in all 6 frames.
        Returns list of ORFResult objects sorted by length.
        """
        try:
            results = []

            for strand, nuc_seq in [(1, self.sequence), (-1, self.sequence.reverse_complement())]:
                for frame in range(3):
                    trans = nuc_seq[frame:].translate(table=self.table_id)
                    proteins = trans.split("*")

                    current_aa_pos = 0
                    for protein in proteins:
                        if len(protein) >= min_len_aa:
                            # Converte a posição do aminoácido para a posição do nucleotídeo
                            nt_start_pos = (current_aa_pos * 3) + frame

                            results.append(ORFResult(
                                id=f"{self.sequence_id}_fr{frame}_st{strand}",
                                frame=frame + 1,
                                strand=strand,
                                length_aa=len(protein),
                                length_bp=len(protein) * 3,
                                protein_seq=str(protein),
                                start_pos=nt_start_pos
                            ))
                        # +1 para compensar o códon de parada (*) removido no split
                        current_aa_pos += len(protein) + 1

            return sorted(results, key=lambda x: x.length_aa, reverse=True)
        except Exception as e:
            raise SequenceAnalysisError(f"Error finding ORFs: {e}")

    def analyze_codon_usage(self) -> CodonMetrics:
        """
        Calculate Codon Usage Bias and GC content.
        """
        try:
            trim = len(self.sequence) // 3 * 3
            coding_seq = self.sequence[:trim]

            codons = [str(coding_seq[i:i + 3]) for i in range(0, len(coding_seq), 3)]
            codon_counts = Counter(codons)
            aa_counts = Counter(Seq(coding_seq).translate(table=self.table_id))

            return CodonMetrics(
                total_codons=sum(codon_counts.values()),
                gc_content=self.sequence.upper().count("G") + self.sequence.upper().count("C"),
                codon_counts=dict(codon_counts),
                amino_acid_counts=dict(aa_counts),
                top_codons=codon_counts.most_common(5)
            )
        except Exception as e:
            raise SequenceAnalysisError(f"Error analyzing codon usage: {e}")