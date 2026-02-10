"""
Enrichment analysis for unmapped sequences in phage display NGS data.

Analyzes sequences that are enriched but don't map to the input library,
identifying potential sources such as:
- Point mutations of library members
- Frameshifts and deletions
- Recombination/chimeras
- Contamination
- Novel sequences
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from collections import Counter, defaultdict
import json
import re

import numpy as np

from .utils import (
    SequenceRecord,
    translate_dna,
    reverse_complement,
    hamming_distance,
)


@dataclass
class SequenceCluster:
    """A cluster of similar sequences."""
    id: int
    representative: str
    members: list[str] = field(default_factory=list)
    count: int = 0
    consensus: Optional[str] = None
    annotation: str = ""

    @property
    def size(self) -> int:
        return len(self.members)

    def __hash__(self):
        return hash(self.id)


@dataclass
class VariantCall:
    """A variant relative to a reference sequence."""
    reference_id: str
    reference_seq: str
    variant_seq: str
    variant_type: str  # 'substitution', 'insertion', 'deletion', 'frameshift'
    position: int
    ref_base: str
    alt_base: str
    count: int = 0


@dataclass
class UnmappedAnalysis:
    """
    Results from analyzing unmapped sequences.
    """
    sample_name: str
    total_unmapped: int = 0
    unique_unmapped: int = 0

    # Sequence clusters
    clusters: list[SequenceCluster] = field(default_factory=list)

    # Potential library variants
    library_variants: list[VariantCall] = field(default_factory=list)

    # Sequence composition
    length_distribution: dict[int, int] = field(default_factory=dict)
    gc_distribution: dict[int, int] = field(default_factory=dict)

    # Motif analysis
    common_motifs: list[tuple[str, int]] = field(default_factory=list)
    amino_acid_frequency: dict[str, int] = field(default_factory=dict)

    # Potential sources
    likely_mutations: int = 0
    likely_chimeras: int = 0
    likely_contamination: int = 0
    novel_sequences: int = 0

    # Top sequences
    top_sequences: list[tuple[str, int]] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "sample_name": self.sample_name,
            "total_unmapped": self.total_unmapped,
            "unique_unmapped": self.unique_unmapped,
            "cluster_count": len(self.clusters),
            "length_distribution": self.length_distribution,
            "gc_distribution": self.gc_distribution,
            "common_motifs": self.common_motifs[:20],
            "source_breakdown": {
                "likely_mutations": self.likely_mutations,
                "likely_chimeras": self.likely_chimeras,
                "likely_contamination": self.likely_contamination,
                "novel_sequences": self.novel_sequences,
            },
            "top_sequences": self.top_sequences[:50],
            "library_variants": [
                {
                    "reference_id": v.reference_id,
                    "variant_type": v.variant_type,
                    "position": v.position,
                    "ref": v.ref_base,
                    "alt": v.alt_base,
                    "count": v.count,
                }
                for v in self.library_variants[:100]
            ],
        }

    def save_json(self, filepath: Path | str):
        """Save analysis to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def summary(self) -> str:
        """Generate text summary of analysis."""
        lines = [
            f"Unmapped Sequence Analysis: {self.sample_name}",
            "=" * 50,
            f"Total unmapped reads: {self.total_unmapped:,}",
            f"Unique sequences: {self.unique_unmapped:,}",
            f"Sequence clusters: {len(self.clusters)}",
            "",
            "Potential Sources:",
            f"  - Library mutations: {self.likely_mutations:,}",
            f"  - Chimeric sequences: {self.likely_chimeras:,}",
            f"  - Contamination: {self.likely_contamination:,}",
            f"  - Novel sequences: {self.novel_sequences:,}",
        ]

        if self.top_sequences:
            lines.extend([
                "",
                "Top 10 Most Frequent Unmapped Sequences:",
            ])
            for seq, count in self.top_sequences[:10]:
                display_seq = seq[:40] + "..." if len(seq) > 40 else seq
                lines.append(f"  {count:>8,}x  {display_seq}")

        if self.common_motifs:
            lines.extend([
                "",
                "Common Motifs (DNA):",
            ])
            for motif, count in self.common_motifs[:10]:
                lines.append(f"  {count:>8,}x  {motif}")

        return "\n".join(lines)


class EnrichmentAnalyzer:
    """
    Analyze enriched sequences that don't map to the input library.

    Performs clustering, variant calling, and source attribution.
    """

    def __init__(
        self,
        library_sequences: Optional[dict[str, str]] = None,
        expected_length: Optional[int] = None,
        similarity_threshold: float = 0.9,
        min_cluster_size: int = 2,
    ):
        """
        Initialize the analyzer.

        Args:
            library_sequences: Dict mapping member_id -> sequence (for variant detection)
            expected_length: Expected length of variable region
            similarity_threshold: Minimum similarity (0-1) for clustering
            min_cluster_size: Minimum sequences to form a cluster
        """
        self.library_sequences = library_sequences or {}
        self.expected_length = expected_length
        self.similarity_threshold = similarity_threshold
        self.min_cluster_size = min_cluster_size

    def analyze(
        self,
        sequences: list[str],
        sample_name: str = "sample",
    ) -> UnmappedAnalysis:
        """
        Analyze a collection of unmapped sequences.

        Args:
            sequences: List of unmapped sequence strings
            sample_name: Name for this analysis

        Returns:
            UnmappedAnalysis with comprehensive results
        """
        result = UnmappedAnalysis(sample_name=sample_name)
        result.total_unmapped = len(sequences)

        # Count unique sequences
        seq_counts = Counter(sequences)
        result.unique_unmapped = len(seq_counts)

        # Top sequences
        result.top_sequences = seq_counts.most_common(100)

        # Basic statistics
        self._calculate_basic_stats(seq_counts, result)

        # Cluster similar sequences
        result.clusters = self._cluster_sequences(seq_counts)

        # Find library variants
        if self.library_sequences:
            result.library_variants = self._find_library_variants(seq_counts)

        # Motif analysis
        self._analyze_motifs(seq_counts, result)

        # Source attribution
        self._attribute_sources(seq_counts, result)

        return result

    def _calculate_basic_stats(self, seq_counts: Counter, result: UnmappedAnalysis):
        """Calculate basic sequence statistics."""
        for seq, count in seq_counts.items():
            length = len(seq)
            result.length_distribution[length] = result.length_distribution.get(length, 0) + count

            gc_count = sum(1 for b in seq.upper() if b in 'GC')
            gc_pct = int(100 * gc_count / length) if length > 0 else 0
            result.gc_distribution[gc_pct] = result.gc_distribution.get(gc_pct, 0) + count

    def _cluster_sequences(self, seq_counts: Counter) -> list[SequenceCluster]:
        """
        Cluster similar sequences using greedy clustering.

        Uses a simple greedy approach - more sophisticated algorithms
        (like CD-HIT) could be used for larger datasets.
        """
        # Sort by abundance (most common first)
        sorted_seqs = sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)

        clusters = []
        assigned = set()
        cluster_id = 0

        for seq, count in sorted_seqs:
            if seq in assigned:
                continue

            # Start new cluster with this sequence as representative
            cluster = SequenceCluster(
                id=cluster_id,
                representative=seq,
                members=[seq],
                count=count,
            )

            # Find similar sequences
            for other_seq, other_count in sorted_seqs:
                if other_seq in assigned or other_seq == seq:
                    continue

                if len(other_seq) != len(seq):
                    continue

                similarity = 1 - (hamming_distance(seq, other_seq) / len(seq))
                if similarity >= self.similarity_threshold:
                    cluster.members.append(other_seq)
                    cluster.count += other_count
                    assigned.add(other_seq)

            assigned.add(seq)

            if cluster.size >= self.min_cluster_size:
                # Calculate consensus
                cluster.consensus = self._calculate_consensus(cluster.members)
                clusters.append(cluster)
                cluster_id += 1

        return clusters

    def _calculate_consensus(self, sequences: list[str]) -> str:
        """Calculate consensus sequence from a list of sequences."""
        if not sequences:
            return ""

        length = len(sequences[0])
        consensus = []

        for pos in range(length):
            bases = [seq[pos].upper() for seq in sequences if pos < len(seq)]
            if bases:
                most_common = Counter(bases).most_common(1)[0][0]
                consensus.append(most_common)

        return ''.join(consensus)

    def _find_library_variants(self, seq_counts: Counter) -> list[VariantCall]:
        """
        Find sequences that are variants of library members.

        Identifies point mutations, small indels, and frameshifts.
        """
        variants = []

        for seq, count in seq_counts.items():
            best_match = None
            best_distance = float('inf')
            best_ref_id = None

            # Find closest library member
            for lib_id, lib_seq in self.library_sequences.items():
                if abs(len(seq) - len(lib_seq)) > 3:  # Allow small indels
                    continue

                if len(seq) == len(lib_seq):
                    dist = hamming_distance(seq, lib_seq)
                else:
                    # For different lengths, use edit distance approximation
                    dist = abs(len(seq) - len(lib_seq)) + self._min_mismatches(seq, lib_seq)

                if dist < best_distance:
                    best_distance = dist
                    best_match = lib_seq
                    best_ref_id = lib_id

            # If close match found, call variant
            if best_match and best_distance > 0 and best_distance <= 3:
                variant = self._call_variant(best_ref_id, best_match, seq, count)
                if variant:
                    variants.append(variant)

        # Sort by count
        variants.sort(key=lambda v: v.count, reverse=True)
        return variants

    def _min_mismatches(self, seq1: str, seq2: str) -> int:
        """Calculate minimum mismatches between sequences of different length."""
        shorter, longer = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)
        min_mm = len(longer)

        for offset in range(len(longer) - len(shorter) + 1):
            mm = sum(1 for i, b in enumerate(shorter) if longer[offset + i] != b)
            min_mm = min(min_mm, mm)

        return min_mm

    def _call_variant(
        self,
        ref_id: str,
        ref_seq: str,
        var_seq: str,
        count: int
    ) -> Optional[VariantCall]:
        """Call a specific variant between reference and variant sequence."""
        if len(ref_seq) == len(var_seq):
            # Find substitution(s)
            for i, (r, v) in enumerate(zip(ref_seq, var_seq)):
                if r != v:
                    return VariantCall(
                        reference_id=ref_id,
                        reference_seq=ref_seq,
                        variant_seq=var_seq,
                        variant_type='substitution',
                        position=i + 1,
                        ref_base=r,
                        alt_base=v,
                        count=count,
                    )
        elif len(var_seq) < len(ref_seq):
            return VariantCall(
                reference_id=ref_id,
                reference_seq=ref_seq,
                variant_seq=var_seq,
                variant_type='deletion',
                position=0,
                ref_base=f"{len(ref_seq) - len(var_seq)}bp",
                alt_base='-',
                count=count,
            )
        else:
            return VariantCall(
                reference_id=ref_id,
                reference_seq=ref_seq,
                variant_seq=var_seq,
                variant_type='insertion',
                position=0,
                ref_base='-',
                alt_base=f"{len(var_seq) - len(ref_seq)}bp",
                count=count,
            )

        return None

    def _analyze_motifs(self, seq_counts: Counter, result: UnmappedAnalysis):
        """Analyze common motifs in unmapped sequences."""
        # DNA k-mers
        kmer_counts: Counter = Counter()
        aa_counts: Counter = Counter()

        for seq, count in seq_counts.items():
            seq = seq.upper()

            # Count 6-mers
            for i in range(len(seq) - 5):
                kmer = seq[i:i+6]
                if 'N' not in kmer:
                    kmer_counts[kmer] += count

            # Translate and count amino acids
            try:
                peptide = translate_dna(seq)
                for aa in peptide:
                    if aa != 'X' and aa != '*':
                        aa_counts[aa] += count
            except Exception:
                pass

        result.common_motifs = kmer_counts.most_common(50)
        result.amino_acid_frequency = dict(aa_counts)

    def _attribute_sources(self, seq_counts: Counter, result: UnmappedAnalysis):
        """
        Attempt to attribute unmapped sequences to likely sources.
        """
        for seq, count in seq_counts.items():
            source = self._classify_sequence_source(seq)

            if source == 'mutation':
                result.likely_mutations += count
            elif source == 'chimera':
                result.likely_chimeras += count
            elif source == 'contamination':
                result.likely_contamination += count
            else:
                result.novel_sequences += count

    def _classify_sequence_source(self, seq: str) -> str:
        """
        Classify the likely source of an unmapped sequence.
        """
        seq = seq.upper()

        # Check if it's a close variant of library member
        if self.library_sequences:
            for lib_seq in self.library_sequences.values():
                if len(seq) == len(lib_seq):
                    dist = hamming_distance(seq, lib_seq)
                    if dist <= 2:
                        return 'mutation'

        # Check for chimeric sequences (partial matches to multiple library members)
        if self.library_sequences and len(seq) >= 20:
            half_len = len(seq) // 2
            first_half = seq[:half_len]
            second_half = seq[half_len:]

            first_match = False
            second_match = False

            for lib_seq in self.library_sequences.values():
                if first_half in lib_seq:
                    first_match = True
                if second_half in lib_seq:
                    second_match = True

            if first_match and second_match:
                return 'chimera'

        # Check length (unexpected length often indicates contamination)
        if self.expected_length:
            if abs(len(seq) - self.expected_length) > 6:
                return 'contamination'

        # Check for common contamination patterns
        contamination_patterns = [
            'AGATCGGAAGAG',  # Illumina adapter
            'CTGTCTCTTATA',  # Nextera adapter
            'AAAAAAAAAA',    # Poly-A
            'TTTTTTTTTT',    # Poly-T
            'GGGGGGGGGG',    # Poly-G (common Illumina artifact)
        ]

        for pattern in contamination_patterns:
            if pattern in seq:
                return 'contamination'

        return 'novel'


def analyze_enrichment_patterns(
    clusters: list[SequenceCluster],
    library_sequences: dict[str, str],
) -> dict:
    """
    Analyze patterns in enriched unmapped sequence clusters.

    Args:
        clusters: List of sequence clusters
        library_sequences: Dictionary of library sequences

    Returns:
        Dictionary with pattern analysis results
    """
    results = {
        'total_clusters': len(clusters),
        'mutation_hotspots': [],
        'enriched_motifs': [],
        'length_variants': [],
    }

    # Analyze mutation patterns
    position_mutations: dict[int, Counter] = defaultdict(Counter)

    for cluster in clusters:
        rep = cluster.representative
        for lib_seq in library_sequences.values():
            if len(rep) == len(lib_seq):
                for pos, (r, l) in enumerate(zip(rep, lib_seq)):
                    if r != l:
                        position_mutations[pos][f"{l}>{r}"] += cluster.count

    # Find mutation hotspots
    hotspots = []
    for pos, mutations in position_mutations.items():
        total = sum(mutations.values())
        if total > 10:
            hotspots.append({
                'position': pos + 1,
                'total_mutations': total,
                'top_mutations': mutations.most_common(3),
            })

    results['mutation_hotspots'] = sorted(hotspots, key=lambda x: x['total_mutations'], reverse=True)

    # Analyze enriched amino acid motifs
    peptide_counts: Counter = Counter()
    for cluster in clusters:
        try:
            peptide = translate_dna(cluster.representative)
            # Extract 3-mer peptide motifs
            for i in range(len(peptide) - 2):
                motif = peptide[i:i+3]
                if '*' not in motif and 'X' not in motif:
                    peptide_counts[motif] += cluster.count
        except Exception:
            pass

    results['enriched_motifs'] = peptide_counts.most_common(20)

    return results


def generate_enrichment_plots(analysis: UnmappedAnalysis, output_dir: Path | str):
    """
    Generate visualization plots for unmapped sequence analysis.

    Args:
        analysis: UnmappedAnalysis results
        output_dir: Directory to save plots
    """
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Length distribution
    if analysis.length_distribution:
        fig, ax = plt.subplots(figsize=(10, 6))
        lengths = sorted(analysis.length_distribution.keys())
        counts = [analysis.length_distribution[l] for l in lengths]
        ax.bar(lengths, counts, width=1, edgecolor='none')
        ax.set_xlabel('Sequence length (bp)')
        ax.set_ylabel('Count')
        ax.set_title(f'Unmapped Sequence Length Distribution - {analysis.sample_name}')
        plt.tight_layout()
        plt.savefig(output_dir / 'unmapped_length_distribution.png', dpi=150)
        plt.close()

    # 2. Source attribution pie chart
    source_data = {
        'Library mutations': analysis.likely_mutations,
        'Chimeric sequences': analysis.likely_chimeras,
        'Contamination': analysis.likely_contamination,
        'Novel sequences': analysis.novel_sequences,
    }
    source_data = {k: v for k, v in source_data.items() if v > 0}

    if source_data:
        fig, ax = plt.subplots(figsize=(8, 8))
        colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
        wedges, texts, autotexts = ax.pie(
            source_data.values(),
            labels=source_data.keys(),
            autopct='%1.1f%%',
            colors=colors[:len(source_data)],
            startangle=90,
        )
        ax.set_title(f'Unmapped Sequence Sources - {analysis.sample_name}')
        plt.tight_layout()
        plt.savefig(output_dir / 'unmapped_sources.png', dpi=150)
        plt.close()

    # 3. Top sequence frequency
    if analysis.top_sequences:
        fig, ax = plt.subplots(figsize=(12, 6))
        top_20 = analysis.top_sequences[:20]
        indices = range(len(top_20))
        counts = [c for _, c in top_20]

        ax.bar(indices, counts)
        ax.set_xlabel('Sequence rank')
        ax.set_ylabel('Count')
        ax.set_title(f'Top 20 Unmapped Sequences - {analysis.sample_name}')
        ax.set_xticks(indices)
        ax.set_xticklabels([str(i+1) for i in indices])
        plt.tight_layout()
        plt.savefig(output_dir / 'top_unmapped_sequences.png', dpi=150)
        plt.close()

    # 4. Cluster size distribution
    if analysis.clusters:
        fig, ax = plt.subplots(figsize=(10, 6))
        cluster_sizes = [c.count for c in analysis.clusters]

        # Bin large clusters
        bins = list(range(0, min(100, max(cluster_sizes) + 10), 5)) + [max(cluster_sizes) + 1]
        ax.hist(cluster_sizes, bins=bins, edgecolor='black')
        ax.set_xlabel('Cluster size (total reads)')
        ax.set_ylabel('Number of clusters')
        ax.set_title(f'Unmapped Sequence Cluster Sizes - {analysis.sample_name}')
        plt.tight_layout()
        plt.savefig(output_dir / 'cluster_sizes.png', dpi=150)
        plt.close()

    # 5. GC content distribution
    if analysis.gc_distribution:
        fig, ax = plt.subplots(figsize=(10, 6))
        gc_bins = sorted(analysis.gc_distribution.keys())
        counts = [analysis.gc_distribution[g] for g in gc_bins]
        ax.bar(gc_bins, counts, width=1, edgecolor='none')
        ax.axvline(x=50, color='r', linestyle='--', alpha=0.5, label='50% GC')
        ax.set_xlabel('GC content (%)')
        ax.set_ylabel('Count')
        ax.set_title(f'Unmapped Sequence GC Distribution - {analysis.sample_name}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(output_dir / 'unmapped_gc_distribution.png', dpi=150)
        plt.close()

    # 6. Amino acid frequency (if available)
    if analysis.amino_acid_frequency:
        fig, ax = plt.subplots(figsize=(12, 6))
        aa_order = 'ACDEFGHIKLMNPQRSTVWY'
        aa_counts = [analysis.amino_acid_frequency.get(aa, 0) for aa in aa_order]

        ax.bar(list(aa_order), aa_counts)
        ax.set_xlabel('Amino acid')
        ax.set_ylabel('Count')
        ax.set_title(f'Amino Acid Frequency in Unmapped Sequences - {analysis.sample_name}')
        plt.tight_layout()
        plt.savefig(output_dir / 'amino_acid_frequency.png', dpi=150)
        plt.close()

    print(f"Enrichment analysis plots saved to {output_dir}")
