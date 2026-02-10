"""
Utility functions and data classes for phage display NGS analysis.
"""

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional
import re


@dataclass
class SequenceRecord:
    """Represents a single sequence record from FASTQ/FASTA."""
    id: str
    sequence: str
    quality: Optional[str] = None  # None for FASTA
    description: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def quality_scores(self) -> Optional[list[int]]:
        """Convert Phred+33 quality string to integer scores."""
        if self.quality is None:
            return None
        return [ord(c) - 33 for c in self.quality]

    @property
    def gc_content(self) -> float:
        """Calculate GC content as a fraction."""
        if not self.sequence:
            return 0.0
        gc_count = sum(1 for base in self.sequence.upper() if base in 'GC')
        return gc_count / len(self.sequence)

    def __hash__(self):
        return hash(self.sequence)


@dataclass
class LibraryMember:
    """Represents a member of the phage display library."""
    id: str
    sequence: str
    peptide: Optional[str] = None
    description: str = ""

    @classmethod
    def from_fasta_record(cls, record: SequenceRecord) -> "LibraryMember":
        return cls(
            id=record.id,
            sequence=record.sequence,
            description=record.description
        )


@dataclass
class ReadCount:
    """Counts for a specific sequence across samples."""
    sequence: str
    library_member_id: Optional[str] = None
    counts: dict[str, int] = field(default_factory=dict)

    @property
    def total_count(self) -> int:
        return sum(self.counts.values())

    def add_count(self, sample_name: str, count: int = 1):
        self.counts[sample_name] = self.counts.get(sample_name, 0) + count


def open_file(filepath: Path | str):
    """Open a file, handling gzip compression automatically."""
    filepath = Path(filepath)
    if filepath.suffix == '.gz':
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def read_fastq(filepath: Path | str) -> Iterator[SequenceRecord]:
    """
    Read sequences from a FASTQ file.

    Handles both plain text and gzip-compressed files.

    Args:
        filepath: Path to FASTQ file (.fastq or .fastq.gz)

    Yields:
        SequenceRecord objects for each read
    """
    with open_file(filepath) as f:
        while True:
            # Read header line
            header = f.readline().strip()
            if not header:
                break
            if not header.startswith('@'):
                raise ValueError(f"Invalid FASTQ format: expected '@', got '{header[0]}'")

            # Parse header
            parts = header[1:].split(None, 1)
            seq_id = parts[0]
            description = parts[1] if len(parts) > 1 else ""

            # Read sequence
            sequence = f.readline().strip()

            # Read separator line
            separator = f.readline().strip()
            if not separator.startswith('+'):
                raise ValueError(f"Invalid FASTQ format: expected '+', got '{separator}'")

            # Read quality
            quality = f.readline().strip()

            yield SequenceRecord(
                id=seq_id,
                sequence=sequence,
                quality=quality,
                description=description
            )


def read_fasta(filepath: Path | str) -> Iterator[SequenceRecord]:
    """
    Read sequences from a FASTA file.

    Handles both plain text and gzip-compressed files.
    Multi-line sequences are supported.

    Args:
        filepath: Path to FASTA file (.fasta, .fa, or .gz variants)

    Yields:
        SequenceRecord objects for each sequence
    """
    with open_file(filepath) as f:
        current_id = None
        current_desc = ""
        current_seq = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Yield previous record if exists
                if current_id is not None:
                    yield SequenceRecord(
                        id=current_id,
                        sequence=''.join(current_seq),
                        description=current_desc
                    )

                # Parse new header
                parts = line[1:].split(None, 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq = []
            else:
                current_seq.append(line)

        # Yield last record
        if current_id is not None:
            yield SequenceRecord(
                id=current_id,
                sequence=''.join(current_seq),
                description=current_desc
            )


def translate_dna(sequence: str, frame: int = 0) -> str:
    """
    Translate DNA sequence to amino acids.

    Args:
        sequence: DNA sequence
        frame: Reading frame (0, 1, or 2)

    Returns:
        Amino acid sequence
    """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    sequence = sequence.upper().replace('U', 'T')
    peptide = []

    for i in range(frame, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        aa = codon_table.get(codon, 'X')  # X for unknown codons
        peptide.append(aa)

    return ''.join(peptide)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))


def hamming_distance(seq1: str, seq2: str) -> int:
    """Calculate Hamming distance between two sequences of equal length."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def extract_variable_region(
    sequence: str,
    upstream_flank: str,
    downstream_flank: str,
    allow_mismatches: int = 1
) -> Optional[str]:
    """
    Extract variable region between conserved flanking sequences.

    Args:
        sequence: Full read sequence
        upstream_flank: Conserved sequence upstream of variable region
        downstream_flank: Conserved sequence downstream of variable region
        allow_mismatches: Number of mismatches allowed in flanking regions

    Returns:
        Extracted variable region or None if flanks not found
    """
    sequence = sequence.upper()
    upstream_flank = upstream_flank.upper()
    downstream_flank = downstream_flank.upper()

    # Find upstream flank with allowed mismatches
    upstream_pos = None
    for i in range(len(sequence) - len(upstream_flank) + 1):
        substr = sequence[i:i + len(upstream_flank)]
        if hamming_distance(substr, upstream_flank) <= allow_mismatches:
            upstream_pos = i + len(upstream_flank)
            break

    if upstream_pos is None:
        return None

    # Find downstream flank with allowed mismatches
    downstream_pos = None
    for i in range(upstream_pos, len(sequence) - len(downstream_flank) + 1):
        substr = sequence[i:i + len(downstream_flank)]
        if hamming_distance(substr, downstream_flank) <= allow_mismatches:
            downstream_pos = i
            break

    if downstream_pos is None:
        return None

    return sequence[upstream_pos:downstream_pos]


def calculate_enrichment_ratio(
    count_selected: int,
    count_input: int,
    total_selected: int,
    total_input: int,
    pseudocount: float = 0.5
) -> float:
    """
    Calculate enrichment ratio between selected and input samples.

    Uses pseudocounts to handle zero counts.

    Args:
        count_selected: Counts in selected/enriched sample
        count_input: Counts in input/naive sample
        total_selected: Total reads in selected sample
        total_input: Total reads in input sample
        pseudocount: Pseudocount for zero-count handling

    Returns:
        Log2 enrichment ratio
    """
    import math

    # Normalize by total counts
    freq_selected = (count_selected + pseudocount) / (total_selected + pseudocount)
    freq_input = (count_input + pseudocount) / (total_input + pseudocount)

    return math.log2(freq_selected / freq_input)
