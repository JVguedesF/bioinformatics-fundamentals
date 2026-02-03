class BioPipelineError(Exception):
    """Classe base para exceções do pipeline."""
    pass

class SequenceAnalysisError(BioPipelineError):
    """Erro ao processar a biologia da sequência."""
    pass

class EmptySequenceError(SequenceAnalysisError):
    """Levantado quando a sequência está vazia ou inválida."""
    pass