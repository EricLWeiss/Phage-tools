"""
Phage Display NGS Analysis Pipeline

A comprehensive toolkit for analyzing Next-Generation Sequencing data from
phage display experiments, including:
- Read quality assessment
- Mapping reads to input libraries
- Enrichment analysis of unmapped sequences
"""

__version__ = "0.1.0"

from .quality_assessment import QualityAssessor, QualityReport
from .library_mapping import LibraryMapper, MappingResult
from .enrichment_analysis import EnrichmentAnalyzer, UnmappedAnalysis
from .utils import read_fastq, read_fasta, SequenceRecord

__all__ = [
    "QualityAssessor",
    "QualityReport",
    "LibraryMapper",
    "MappingResult",
    "EnrichmentAnalyzer",
    "UnmappedAnalysis",
    "read_fastq",
    "read_fasta",
    "SequenceRecord",
]
