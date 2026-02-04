from dataclasses import dataclass
from typing import Optional
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data.CodonTable import TranslationError

from config import settings
from utils.exceptions import SequenceAnalysisError, EmptySequenceError

try:
    import RNA
    HAS_VIENNA_RNA = True
except ImportError:
    HAS_VIENNA_RNA = False

DNA_BASES = set("ATGCN")
RNA_BASES = set("AUGCN")

@dataclass
class AnalysisResult:
    filename: str
    sequence_id: str
    molecule_type: str
    length: int
    gc_content: Optional[float] = None
    melting_temp: Optional[float] = None
    mfe: Optional[float] = None
    isoelectric_point: Optional[float] = None
    molecular_weight: Optional[float] = None
    stability_index: Optional[float] = None
    timestamp: str = None
    secondary_structure: Optional[str] = None

    def __post_init__(self):
        if self.timestamp is None:
            self.timestamp = datetime.now().isoformat()

def calculate_dna_metrics(seq):
    length = len(seq)
    gc = gc_fraction(seq) * 100
    tm = mt.Tm_NN(seq, strict=False)
    return {
        "length": length,
        "gc_content": gc,
        "tm": tm
    }

def get_folding_data(sequence_str):
    if HAS_VIENNA_RNA:
        structure, mfe = RNA.fold(sequence_str)
        return structure, mfe
    return None, None

def analyze_protein(sequence_str):
    clean_seq = sequence_str.rstrip('*')
    if not clean_seq:
        return None
    try:
        X = ProteinAnalysis(clean_seq)
        return {
            "length": len(clean_seq),
            "molecular_weight": X.molecular_weight(),
            "isoelectric_point": X.isoelectric_point(),
            "gravy": X.gravy(),
            "instability_index": X.instability_index()
        }
    except ValueError:
        return None

class CentralDogmaAnalyzer:
    @staticmethod
    def detect_molecule_type(sequence: str) -> str:
        unique_chars = set(sequence.upper()) - set("\n\r\t ")

        if len(unique_chars - DNA_BASES - RNA_BASES) > 0:
            return 'PROTEIN'
        if 'T' in unique_chars and 'U' not in unique_chars:
            return 'DNA'
        if 'U' in unique_chars and 'T' not in unique_chars:
            return 'RNA'
        return 'DNA'

    def process_sequence(self, record: SeqRecord, filename: str) -> AnalysisResult:
        if not record.seq:
            raise EmptySequenceError(f"Sequence in {filename} is empty.")

        mol_type = self.detect_molecule_type(str(record.seq)[:1000])
        
        result = AnalysisResult(
            filename=filename,
            sequence_id=record.id,
            molecule_type=mol_type,
            length=len(record.seq)
        )

        rna_seq = None
        protein_seq = record.seq

        try:
            if mol_type == 'DNA':
                dna_metrics = calculate_dna_metrics(record.seq)
                result.gc_content = dna_metrics.get('gc_content')
                result.melting_temp = dna_metrics.get('tm')
                rna_seq = record.seq.transcribe()
            elif mol_type == 'RNA':
                rna_seq = record.seq

            if rna_seq and HAS_VIENNA_RNA:
                structure, mfe = get_folding_data(str(rna_seq))
                if structure:
                    result.mfe = mfe
                    result.secondary_structure = structure
                try:
                    protein_seq = rna_seq.translate(to_stop=True)
                except (ValueError, TranslationError):
                    protein_seq = ""

            if len(protein_seq) >= settings.MIN_PROTEIN_LENGTH:
                prot_data = analyze_protein(str(protein_seq))
                if prot_data:
                    result.isoelectric_point = prot_data.get('isoelectric_point')
                    result.molecular_weight = prot_data.get('molecular_weight')
                    result.stability_index = prot_data.get('instability_index')
        except Exception as e:
            raise SequenceAnalysisError(f"Failed to analyze {mol_type} sequence: {str(e)}")

        return result