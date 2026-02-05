"""
Library mapping module for phage display NGS analysis.

Maps sequencing reads to known library members and calculates enrichment.
Similar to RNA-seq read counting but adapted for phage display experiments.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Iterator
from collections import defaultdict
import json
import csv

import numpy as np

from .utils import (
    read_fastq,
    read_fasta,
    SequenceRecord,
    LibraryMember,
    extract_variable_region,
    reverse_complement,
    hamming_distance,
    translate_dna,
    calculate_enrichment_ratio,
)


@dataclass
class MappingStats:
    """Statistics from a single mapping run."""
    sample_name: str
    total_reads: int = 0
    mapped_reads: int = 0
    unmapped_reads: int = 0
    ambiguous_reads: int = 0  # Reads mapping to multiple library members equally well
    low_quality_reads: int = 0  # Reads filtered due to quality

    @property
    def mapping_rate(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return self.mapped_reads / self.total_reads

    def summary(self) -> str:
        return (
            f"Sample: {self.sample_name}\n"
            f"  Total reads: {self.total_reads:,}\n"
            f"  Mapped: {self.mapped_reads:,} ({self.mapping_rate:.1%})\n"
            f"  Unmapped: {self.unmapped_reads:,}\n"
            f"  Ambiguous: {self.ambiguous_reads:,}\n"
            f"  Low quality filtered: {self.low_quality_reads:,}"
        )


@dataclass
class MappingResult:
    """
    Complete results from library mapping.

    Contains counts per library member and enrichment calculations.
    """
    library_name: str
    samples: list[str] = field(default_factory=list)
    stats: dict[str, MappingStats] = field(default_factory=dict)

    # Count matrices: library_member_id -> {sample -> count}
    counts: dict[str, dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))

    # Normalized counts (CPM - counts per million)
    cpm: dict[str, dict[str, float]] = field(default_factory=dict)

    # Enrichment ratios (log2 fold change)
    enrichment: dict[str, dict[str, float]] = field(default_factory=dict)

    # Unmapped sequences for further analysis
    unmapped_sequences: dict[str, list[str]] = field(default_factory=lambda: defaultdict(list))

    def get_count_matrix(self) -> tuple[list[str], list[str], np.ndarray]:
        """
        Return count matrix as numpy array.

        Returns:
            Tuple of (library_member_ids, sample_names, count_matrix)
        """
        members = sorted(self.counts.keys())
        samples = self.samples

        matrix = np.zeros((len(members), len(samples)), dtype=int)
        for i, member in enumerate(members):
            for j, sample in enumerate(samples):
                matrix[i, j] = self.counts[member].get(sample, 0)

        return members, samples, matrix

    def calculate_cpm(self):
        """Calculate counts per million (CPM) normalization."""
        # Get total counts per sample
        sample_totals = defaultdict(int)
        for member_counts in self.counts.values():
            for sample, count in member_counts.items():
                sample_totals[sample] += count

        # Calculate CPM
        self.cpm = {}
        for member, member_counts in self.counts.items():
            self.cpm[member] = {}
            for sample in self.samples:
                count = member_counts.get(sample, 0)
                total = sample_totals.get(sample, 1)
                self.cpm[member][sample] = (count / total) * 1_000_000

    def calculate_enrichment(self, input_sample: str, selected_samples: Optional[list[str]] = None):
        """
        Calculate log2 enrichment ratios relative to input sample.

        Args:
            input_sample: Name of the input/naive library sample
            selected_samples: List of selected samples to compare (defaults to all non-input)
        """
        if selected_samples is None:
            selected_samples = [s for s in self.samples if s != input_sample]

        # Ensure CPM is calculated
        if not self.cpm:
            self.calculate_cpm()

        # Get totals for enrichment calculation
        sample_totals = {}
        for sample in self.samples:
            sample_totals[sample] = sum(
                self.counts[m].get(sample, 0) for m in self.counts
            )

        self.enrichment = {}
        for member in self.counts:
            self.enrichment[member] = {}
            count_input = self.counts[member].get(input_sample, 0)

            for sample in selected_samples:
                count_selected = self.counts[member].get(sample, 0)
                self.enrichment[member][sample] = calculate_enrichment_ratio(
                    count_selected,
                    count_input,
                    sample_totals.get(sample, 1),
                    sample_totals.get(input_sample, 1)
                )

    def get_top_enriched(
        self,
        sample: str,
        n: int = 100,
        min_count: int = 10
    ) -> list[tuple[str, float, int]]:
        """
        Get top enriched library members for a sample.

        Args:
            sample: Sample name
            n: Number of top members to return
            min_count: Minimum count threshold

        Returns:
            List of (member_id, enrichment_ratio, count) tuples
        """
        if not self.enrichment:
            return []

        results = []
        for member, enrichments in self.enrichment.items():
            if sample not in enrichments:
                continue
            count = self.counts[member].get(sample, 0)
            if count >= min_count:
                results.append((member, enrichments[sample], count))

        results.sort(key=lambda x: x[1], reverse=True)
        return results[:n]

    def save_counts(self, filepath: Path | str):
        """Save count matrix to TSV file."""
        filepath = Path(filepath)

        members, samples, matrix = self.get_count_matrix()

        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['library_member'] + samples)
            for i, member in enumerate(members):
                writer.writerow([member] + list(matrix[i, :]))

    def save_enrichment(self, filepath: Path | str):
        """Save enrichment results to TSV file."""
        filepath = Path(filepath)

        if not self.enrichment:
            return

        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            # Get all samples from enrichment
            all_samples = set()
            for member_enrichment in self.enrichment.values():
                all_samples.update(member_enrichment.keys())
            samples = sorted(all_samples)

            writer.writerow(['library_member'] + [f'log2FC_{s}' for s in samples])
            for member in sorted(self.enrichment.keys()):
                row = [member]
                for sample in samples:
                    row.append(f"{self.enrichment[member].get(sample, 0):.4f}")
                writer.writerow(row)

    def save_unmapped(self, filepath: Path | str, sample: Optional[str] = None):
        """Save unmapped sequences to FASTA file."""
        filepath = Path(filepath)

        samples_to_save = [sample] if sample else self.unmapped_sequences.keys()

        with open(filepath, 'w') as f:
            for s in samples_to_save:
                seqs = self.unmapped_sequences.get(s, [])
                for i, seq in enumerate(seqs):
                    f.write(f">{s}_unmapped_{i+1}\n{seq}\n")

    def summary(self) -> str:
        """Generate summary of mapping results."""
        lines = [
            f"Library Mapping Results: {self.library_name}",
            "=" * 50,
            f"Samples: {len(self.samples)}",
            f"Library members with reads: {len(self.counts)}",
        ]

        for sample in self.samples:
            if sample in self.stats:
                lines.append(f"\n{self.stats[sample].summary()}")

        return "\n".join(lines)


class LibraryMapper:
    """
    Map NGS reads to phage display library members.

    Supports exact matching and fuzzy matching with allowed mismatches.
    Handles variable regions by extracting them using flanking sequences.
    """

    def __init__(
        self,
        library_path: Optional[Path | str] = None,
        upstream_flank: Optional[str] = None,
        downstream_flank: Optional[str] = None,
        max_mismatches: int = 0,
        min_quality: int = 20,
        check_reverse_complement: bool = True,
    ):
        """
        Initialize the library mapper.

        Args:
            library_path: Path to FASTA file with library sequences
            upstream_flank: Conserved sequence upstream of variable region
            downstream_flank: Conserved sequence downstream of variable region
            max_mismatches: Maximum allowed mismatches for mapping
            min_quality: Minimum average quality score for a read
            check_reverse_complement: Also check reverse complement of reads
        """
        self.upstream_flank = upstream_flank
        self.downstream_flank = downstream_flank
        self.max_mismatches = max_mismatches
        self.min_quality = min_quality
        self.check_reverse_complement = check_reverse_complement

        # Library sequences indexed by sequence
        self.library: dict[str, LibraryMember] = {}
        self.library_by_id: dict[str, LibraryMember] = {}

        if library_path:
            self.load_library(library_path)

    def load_library(self, library_path: Path | str):
        """
        Load library sequences from FASTA file.

        Args:
            library_path: Path to FASTA file
        """
        self.library = {}
        self.library_by_id = {}

        for record in read_fasta(library_path):
            member = LibraryMember.from_fasta_record(record)
            seq = member.sequence.upper()

            # If flanking sequences specified, extract variable region
            if self.upstream_flank and self.downstream_flank:
                var_region = extract_variable_region(
                    seq,
                    self.upstream_flank,
                    self.downstream_flank,
                    allow_mismatches=0
                )
                if var_region:
                    seq = var_region

            self.library[seq] = member
            self.library_by_id[member.id] = member

        print(f"Loaded {len(self.library)} library members")

    def add_library_sequence(self, seq_id: str, sequence: str, peptide: Optional[str] = None):
        """Add a single sequence to the library."""
        member = LibraryMember(id=seq_id, sequence=sequence, peptide=peptide)
        self.library[sequence.upper()] = member
        self.library_by_id[seq_id] = member

    def _extract_query_sequence(self, read: SequenceRecord) -> Optional[str]:
        """Extract the query sequence from a read (variable region if flanks specified)."""
        seq = read.sequence.upper()

        if self.upstream_flank and self.downstream_flank:
            var_region = extract_variable_region(
                seq,
                self.upstream_flank,
                self.downstream_flank,
                allow_mismatches=1
            )
            if var_region:
                return var_region

            # Try reverse complement
            if self.check_reverse_complement:
                rc_seq = reverse_complement(seq)
                var_region = extract_variable_region(
                    rc_seq,
                    self.upstream_flank,
                    self.downstream_flank,
                    allow_mismatches=1
                )
                if var_region:
                    return var_region

            return None
        else:
            return seq

    def _find_match(self, query: str) -> tuple[Optional[str], int]:
        """
        Find matching library member for a query sequence.

        Returns:
            Tuple of (member_id or None, number of mismatches)
        """
        query = query.upper()

        # Exact match first
        if query in self.library:
            return self.library[query].id, 0

        # Fuzzy matching if allowed
        if self.max_mismatches > 0:
            best_match = None
            best_mismatches = self.max_mismatches + 1
            match_count = 0

            for lib_seq, member in self.library.items():
                if len(lib_seq) != len(query):
                    continue

                mismatches = hamming_distance(query, lib_seq)
                if mismatches < best_mismatches:
                    best_match = member.id
                    best_mismatches = mismatches
                    match_count = 1
                elif mismatches == best_mismatches:
                    match_count += 1

            if best_mismatches <= self.max_mismatches:
                if match_count == 1:
                    return best_match, best_mismatches
                else:
                    return None, -1  # Ambiguous match

        return None, -1

    def _passes_quality_filter(self, read: SequenceRecord) -> bool:
        """Check if read passes quality filter."""
        scores = read.quality_scores
        if scores is None:
            return True
        return np.mean(scores) >= self.min_quality

    def map_sample(
        self,
        fastq_path: Path | str,
        sample_name: Optional[str] = None,
        max_unmapped: int = 100000,
    ) -> tuple[dict[str, int], MappingStats, list[str]]:
        """
        Map reads from a single FASTQ file to library.

        Args:
            fastq_path: Path to FASTQ file
            sample_name: Name for this sample
            max_unmapped: Maximum unmapped sequences to store

        Returns:
            Tuple of (counts_dict, mapping_stats, unmapped_sequences)
        """
        fastq_path = Path(fastq_path)
        if sample_name is None:
            sample_name = fastq_path.stem.replace('.fastq', '')

        stats = MappingStats(sample_name=sample_name)
        counts: dict[str, int] = defaultdict(int)
        unmapped: list[str] = []

        for read in read_fastq(fastq_path):
            stats.total_reads += 1

            # Quality filter
            if not self._passes_quality_filter(read):
                stats.low_quality_reads += 1
                continue

            # Extract query sequence
            query = self._extract_query_sequence(read)
            if query is None:
                stats.unmapped_reads += 1
                if len(unmapped) < max_unmapped:
                    unmapped.append(read.sequence)
                continue

            # Find match
            member_id, mismatches = self._find_match(query)

            if member_id is not None:
                counts[member_id] += 1
                stats.mapped_reads += 1
            elif mismatches == -1 and self.max_mismatches > 0:
                stats.ambiguous_reads += 1
            else:
                stats.unmapped_reads += 1
                if len(unmapped) < max_unmapped:
                    unmapped.append(query)

            # Progress indicator
            if stats.total_reads % 100000 == 0:
                print(f"  Processed {stats.total_reads:,} reads...")

        print(f"Completed {sample_name}: {stats.mapped_reads:,}/{stats.total_reads:,} mapped ({stats.mapping_rate:.1%})")
        return dict(counts), stats, unmapped

    def map_samples(
        self,
        fastq_paths: dict[str, Path | str],
        library_name: str = "library",
    ) -> MappingResult:
        """
        Map multiple samples to the library.

        Args:
            fastq_paths: Dictionary mapping sample names to FASTQ paths
            library_name: Name for this mapping result

        Returns:
            MappingResult with counts and statistics
        """
        result = MappingResult(library_name=library_name)
        result.samples = list(fastq_paths.keys())

        for sample_name, fastq_path in fastq_paths.items():
            print(f"\nMapping {sample_name}...")
            counts, stats, unmapped = self.map_sample(fastq_path, sample_name)

            result.stats[sample_name] = stats
            result.unmapped_sequences[sample_name] = unmapped

            for member_id, count in counts.items():
                result.counts[member_id][sample_name] = count

        # Calculate CPM normalization
        result.calculate_cpm()

        return result


def compare_samples(
    result: MappingResult,
    input_sample: str,
    selected_sample: str,
    min_count: int = 10,
    min_fold_change: float = 2.0,
) -> list[dict]:
    """
    Compare two samples and identify enriched sequences.

    Args:
        result: MappingResult from mapping
        input_sample: Name of input/naive sample
        selected_sample: Name of selected/enriched sample
        min_count: Minimum count in selected sample
        min_fold_change: Minimum log2 fold change

    Returns:
        List of dictionaries with enrichment data
    """
    result.calculate_enrichment(input_sample, [selected_sample])

    enriched = []
    for member_id, enrichments in result.enrichment.items():
        if selected_sample not in enrichments:
            continue

        log2fc = enrichments[selected_sample]
        count_selected = result.counts[member_id].get(selected_sample, 0)
        count_input = result.counts[member_id].get(input_sample, 0)
        cpm_selected = result.cpm.get(member_id, {}).get(selected_sample, 0)

        if count_selected >= min_count and log2fc >= min_fold_change:
            enriched.append({
                'member_id': member_id,
                'log2_fold_change': log2fc,
                'count_input': count_input,
                'count_selected': count_selected,
                'cpm_selected': cpm_selected,
            })

    # Sort by fold change
    enriched.sort(key=lambda x: x['log2_fold_change'], reverse=True)
    return enriched


def generate_mapping_plots(result: MappingResult, output_dir: Path | str):
    """
    Generate visualization plots for mapping results.

    Args:
        result: MappingResult to visualize
        output_dir: Directory to save plots
    """
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Mapping rate bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    samples = list(result.stats.keys())
    mapping_rates = [result.stats[s].mapping_rate * 100 for s in samples]
    colors = ['green' if r > 50 else 'orange' if r > 20 else 'red' for r in mapping_rates]

    bars = ax.bar(samples, mapping_rates, color=colors)
    ax.set_ylabel('Mapping rate (%)')
    ax.set_xlabel('Sample')
    ax.set_title('Read Mapping Rates')
    ax.set_ylim(0, 100)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / 'mapping_rates.png', dpi=150)
    plt.close()

    # 2. Read count stacked bar
    fig, ax = plt.subplots(figsize=(10, 6))
    mapped = [result.stats[s].mapped_reads for s in samples]
    unmapped = [result.stats[s].unmapped_reads for s in samples]
    ambiguous = [result.stats[s].ambiguous_reads for s in samples]
    low_qual = [result.stats[s].low_quality_reads for s in samples]

    x = range(len(samples))
    width = 0.6

    ax.bar(x, mapped, width, label='Mapped', color='green')
    ax.bar(x, unmapped, width, bottom=mapped, label='Unmapped', color='red')
    ax.bar(x, ambiguous, width, bottom=[m+u for m,u in zip(mapped, unmapped)],
           label='Ambiguous', color='orange')
    ax.bar(x, low_qual, width, bottom=[m+u+a for m,u,a in zip(mapped, unmapped, ambiguous)],
           label='Low quality', color='gray')

    ax.set_ylabel('Number of reads')
    ax.set_xlabel('Sample')
    ax.set_title('Read Classification')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / 'read_classification.png', dpi=150)
    plt.close()

    # 3. Library member coverage heatmap (top members)
    if len(result.samples) >= 2 and result.cpm:
        # Get top 50 most abundant members
        member_totals = {m: sum(c.values()) for m, c in result.counts.items()}
        top_members = sorted(member_totals, key=member_totals.get, reverse=True)[:50]

        if top_members:
            fig, ax = plt.subplots(figsize=(12, 10))

            matrix = np.zeros((len(top_members), len(result.samples)))
            for i, member in enumerate(top_members):
                for j, sample in enumerate(result.samples):
                    # Use log10(CPM + 1) for visualization
                    cpm_val = result.cpm.get(member, {}).get(sample, 0)
                    matrix[i, j] = np.log10(cpm_val + 1)

            im = ax.imshow(matrix, aspect='auto', cmap='viridis')
            ax.set_xticks(range(len(result.samples)))
            ax.set_xticklabels(result.samples, rotation=45, ha='right')
            ax.set_ylabel('Library member')
            ax.set_title('Library Member Abundance (log10 CPM)')
            plt.colorbar(im, ax=ax, label='log10(CPM + 1)')
            plt.tight_layout()
            plt.savefig(output_dir / 'abundance_heatmap.png', dpi=150)
            plt.close()

    # 4. Enrichment volcano plot (if enrichment calculated)
    if result.enrichment and len(result.enrichment) > 0:
        for sample in list(result.enrichment.values())[0].keys():
            fig, ax = plt.subplots(figsize=(10, 8))

            log2fc = []
            log_counts = []
            for member, enrichments in result.enrichment.items():
                if sample in enrichments:
                    fc = enrichments[sample]
                    count = result.counts[member].get(sample, 0)
                    if count > 0:
                        log2fc.append(fc)
                        log_counts.append(np.log10(count + 1))

            if log2fc:
                colors = ['red' if fc > 2 else 'blue' if fc < -2 else 'gray'
                          for fc in log2fc]
                ax.scatter(log2fc, log_counts, c=colors, alpha=0.5, s=20)
                ax.axvline(x=2, color='gray', linestyle='--', alpha=0.5)
                ax.axvline(x=-2, color='gray', linestyle='--', alpha=0.5)
                ax.set_xlabel('Log2 Fold Change')
                ax.set_ylabel('Log10(Count + 1)')
                ax.set_title(f'Enrichment Plot - {sample}')
                plt.tight_layout()
                plt.savefig(output_dir / f'enrichment_{sample}.png', dpi=150)
                plt.close()

    print(f"Mapping plots saved to {output_dir}")
