"""
Read quality assessment for phage display NGS data.

Provides FastQC-like quality metrics including:
- Per-base quality scores
- Read length distribution
- GC content analysis
- Quality score distribution
- Sequence duplication levels
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from collections import Counter
import json

import numpy as np

from .utils import read_fastq, SequenceRecord


@dataclass
class QualityReport:
    """
    Container for read quality assessment results.

    Similar to a FastQC report but tailored for phage display data.
    """
    sample_name: str
    total_reads: int = 0
    total_bases: int = 0

    # Read length statistics
    min_length: int = 0
    max_length: int = 0
    mean_length: float = 0.0
    length_distribution: dict[int, int] = field(default_factory=dict)

    # Quality score statistics
    mean_quality: float = 0.0
    per_position_quality: list[dict] = field(default_factory=list)  # {pos, mean, median, q1, q3, min, max}
    quality_score_distribution: dict[int, int] = field(default_factory=dict)

    # GC content
    mean_gc_content: float = 0.0
    gc_distribution: dict[int, int] = field(default_factory=dict)  # GC% bin -> count

    # Sequence duplication
    unique_sequences: int = 0
    duplication_levels: dict[int, int] = field(default_factory=dict)  # duplication level -> count
    most_common_sequences: list[tuple[str, int]] = field(default_factory=list)

    # N content
    reads_with_n: int = 0
    per_position_n_content: list[float] = field(default_factory=list)

    # Phage display specific
    variable_region_lengths: dict[int, int] = field(default_factory=dict)
    reads_with_valid_flanks: int = 0

    # Quality flags
    warnings: list[str] = field(default_factory=list)
    failures: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert report to dictionary for JSON serialization."""
        return {
            "sample_name": self.sample_name,
            "total_reads": self.total_reads,
            "total_bases": self.total_bases,
            "read_length": {
                "min": self.min_length,
                "max": self.max_length,
                "mean": self.mean_length,
                "distribution": self.length_distribution,
            },
            "quality_scores": {
                "mean": self.mean_quality,
                "per_position": self.per_position_quality,
                "distribution": self.quality_score_distribution,
            },
            "gc_content": {
                "mean": self.mean_gc_content,
                "distribution": self.gc_distribution,
            },
            "sequence_duplication": {
                "unique_sequences": self.unique_sequences,
                "duplication_levels": self.duplication_levels,
                "most_common": self.most_common_sequences[:10],
            },
            "n_content": {
                "reads_with_n": self.reads_with_n,
                "per_position": self.per_position_n_content,
            },
            "phage_display": {
                "variable_region_lengths": self.variable_region_lengths,
                "reads_with_valid_flanks": self.reads_with_valid_flanks,
            },
            "warnings": self.warnings,
            "failures": self.failures,
        }

    def save_json(self, filepath: Path | str):
        """Save report to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def summary(self) -> str:
        """Generate a text summary of the quality report."""
        lines = [
            f"Quality Report: {self.sample_name}",
            "=" * 50,
            f"Total reads: {self.total_reads:,}",
            f"Total bases: {self.total_bases:,}",
            f"Read length: {self.min_length}-{self.max_length} bp (mean: {self.mean_length:.1f})",
            f"Mean quality score: {self.mean_quality:.1f}",
            f"Mean GC content: {self.mean_gc_content:.1%}",
            f"Unique sequences: {self.unique_sequences:,} ({100*self.unique_sequences/max(1,self.total_reads):.1f}%)",
            f"Reads with N bases: {self.reads_with_n:,}",
        ]

        if self.reads_with_valid_flanks > 0:
            lines.append(f"Reads with valid flanks: {self.reads_with_valid_flanks:,}")

        if self.warnings:
            lines.append("\nWarnings:")
            for w in self.warnings:
                lines.append(f"  - {w}")

        if self.failures:
            lines.append("\nFailures:")
            for f in self.failures:
                lines.append(f"  - {f}")

        return "\n".join(lines)


class QualityAssessor:
    """
    Assess quality of NGS reads from phage display experiments.

    Usage:
        assessor = QualityAssessor()
        report = assessor.assess("sample.fastq.gz")
        print(report.summary())
    """

    def __init__(
        self,
        upstream_flank: Optional[str] = None,
        downstream_flank: Optional[str] = None,
        sample_size: Optional[int] = None,
        seed: int = 42
    ):
        """
        Initialize the quality assessor.

        Args:
            upstream_flank: Conserved sequence upstream of variable region (for phage display)
            downstream_flank: Conserved sequence downstream of variable region
            sample_size: If set, randomly sample this many reads for faster analysis
            seed: Random seed for reproducible sampling
        """
        self.upstream_flank = upstream_flank
        self.downstream_flank = downstream_flank
        self.sample_size = sample_size
        self.seed = seed

    def assess(
        self,
        fastq_path: Path | str,
        sample_name: Optional[str] = None
    ) -> QualityReport:
        """
        Perform quality assessment on a FASTQ file.

        Args:
            fastq_path: Path to FASTQ file
            sample_name: Name for this sample (defaults to filename)

        Returns:
            QualityReport with comprehensive quality metrics
        """
        fastq_path = Path(fastq_path)
        if sample_name is None:
            sample_name = fastq_path.stem.replace('.fastq', '')

        # Collect all reads (or sample if specified)
        reads = list(self._load_reads(fastq_path))

        if not reads:
            report = QualityReport(sample_name=sample_name)
            report.failures.append("No reads found in file")
            return report

        # Calculate all metrics
        report = QualityReport(sample_name=sample_name)
        report.total_reads = len(reads)

        self._calculate_length_stats(reads, report)
        self._calculate_quality_stats(reads, report)
        self._calculate_gc_stats(reads, report)
        self._calculate_duplication_stats(reads, report)
        self._calculate_n_content(reads, report)

        if self.upstream_flank and self.downstream_flank:
            self._calculate_phage_display_stats(reads, report)

        self._evaluate_quality_flags(report)

        return report

    def _load_reads(self, fastq_path: Path) -> list[SequenceRecord]:
        """Load reads, optionally sampling."""
        all_reads = list(read_fastq(fastq_path))

        if self.sample_size and len(all_reads) > self.sample_size:
            rng = np.random.default_rng(self.seed)
            indices = rng.choice(len(all_reads), size=self.sample_size, replace=False)
            return [all_reads[i] for i in indices]

        return all_reads

    def _calculate_length_stats(self, reads: list[SequenceRecord], report: QualityReport):
        """Calculate read length statistics."""
        lengths = [read.length for read in reads]
        report.total_bases = sum(lengths)
        report.min_length = min(lengths)
        report.max_length = max(lengths)
        report.mean_length = np.mean(lengths)
        report.length_distribution = dict(Counter(lengths))

    def _calculate_quality_stats(self, reads: list[SequenceRecord], report: QualityReport):
        """Calculate quality score statistics."""
        # All quality scores
        all_scores = []
        # Per-position quality scores
        max_len = report.max_length
        position_scores = [[] for _ in range(max_len)]

        for read in reads:
            scores = read.quality_scores
            if scores is None:
                continue
            all_scores.extend(scores)
            for pos, score in enumerate(scores):
                position_scores[pos].append(score)

        if all_scores:
            report.mean_quality = np.mean(all_scores)
            report.quality_score_distribution = dict(Counter(all_scores))

        # Per-position statistics
        for pos, scores in enumerate(position_scores):
            if not scores:
                continue
            scores_arr = np.array(scores)
            report.per_position_quality.append({
                "position": pos + 1,
                "mean": float(np.mean(scores_arr)),
                "median": float(np.median(scores_arr)),
                "q1": float(np.percentile(scores_arr, 25)),
                "q3": float(np.percentile(scores_arr, 75)),
                "min": int(np.min(scores_arr)),
                "max": int(np.max(scores_arr)),
                "count": len(scores),
            })

    def _calculate_gc_stats(self, reads: list[SequenceRecord], report: QualityReport):
        """Calculate GC content statistics."""
        gc_contents = [read.gc_content for read in reads]
        report.mean_gc_content = np.mean(gc_contents)

        # Bin GC content into 1% bins
        gc_bins = Counter(int(gc * 100) for gc in gc_contents)
        report.gc_distribution = dict(gc_bins)

    def _calculate_duplication_stats(self, reads: list[SequenceRecord], report: QualityReport):
        """Calculate sequence duplication statistics."""
        sequence_counts = Counter(read.sequence for read in reads)
        report.unique_sequences = len(sequence_counts)

        # Duplication levels
        dup_levels = Counter(count for count in sequence_counts.values())
        report.duplication_levels = dict(dup_levels)

        # Most common sequences
        report.most_common_sequences = sequence_counts.most_common(100)

    def _calculate_n_content(self, reads: list[SequenceRecord], report: QualityReport):
        """Calculate N base content."""
        max_len = report.max_length
        position_n_counts = [0] * max_len
        position_totals = [0] * max_len

        for read in reads:
            seq = read.sequence.upper()
            if 'N' in seq:
                report.reads_with_n += 1
            for pos, base in enumerate(seq):
                position_totals[pos] += 1
                if base == 'N':
                    position_n_counts[pos] += 1

        report.per_position_n_content = [
            n / total if total > 0 else 0
            for n, total in zip(position_n_counts, position_totals)
        ]

    def _calculate_phage_display_stats(self, reads: list[SequenceRecord], report: QualityReport):
        """Calculate phage display specific statistics."""
        from .utils import extract_variable_region

        var_lengths = []
        for read in reads:
            var_region = extract_variable_region(
                read.sequence,
                self.upstream_flank,
                self.downstream_flank
            )
            if var_region is not None:
                report.reads_with_valid_flanks += 1
                var_lengths.append(len(var_region))

        report.variable_region_lengths = dict(Counter(var_lengths))

    def _evaluate_quality_flags(self, report: QualityReport):
        """Evaluate quality and set warning/failure flags."""
        # Check mean quality score
        if report.mean_quality < 20:
            report.failures.append(f"Low mean quality score: {report.mean_quality:.1f}")
        elif report.mean_quality < 28:
            report.warnings.append(f"Moderate mean quality score: {report.mean_quality:.1f}")

        # Check duplication levels
        dup_rate = 1 - (report.unique_sequences / max(1, report.total_reads))
        if dup_rate > 0.9:
            report.warnings.append(f"Very high duplication rate: {dup_rate:.1%}")
        elif dup_rate > 0.7:
            report.warnings.append(f"High duplication rate: {dup_rate:.1%}")

        # Check N content
        if report.reads_with_n > 0.1 * report.total_reads:
            report.warnings.append(f"High proportion of reads with N bases: {report.reads_with_n/report.total_reads:.1%}")

        # Check per-position quality dropoff
        if report.per_position_quality:
            last_positions = report.per_position_quality[-10:]
            if last_positions:
                end_quality = np.mean([p["mean"] for p in last_positions])
                if end_quality < 20:
                    report.warnings.append(f"Low quality at read ends: mean Q{end_quality:.1f}")

        # Check GC content (typically 40-60% for balanced libraries)
        if report.mean_gc_content < 0.3 or report.mean_gc_content > 0.7:
            report.warnings.append(f"Unusual GC content: {report.mean_gc_content:.1%}")

        # Phage display specific checks
        if self.upstream_flank and self.downstream_flank:
            flank_rate = report.reads_with_valid_flanks / max(1, report.total_reads)
            if flank_rate < 0.5:
                report.warnings.append(f"Low proportion of reads with valid flanking sequences: {flank_rate:.1%}")
            if flank_rate < 0.1:
                report.failures.append(f"Very few reads with valid flanking sequences: {flank_rate:.1%}")


def generate_quality_plots(report: QualityReport, output_dir: Path | str):
    """
    Generate quality visualization plots using matplotlib.

    Args:
        report: QualityReport to visualize
        output_dir: Directory to save plots
    """
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Per-position quality boxplot
    if report.per_position_quality:
        fig, ax = plt.subplots(figsize=(12, 6))
        positions = [p["position"] for p in report.per_position_quality]
        means = [p["mean"] for p in report.per_position_quality]
        q1s = [p["q1"] for p in report.per_position_quality]
        q3s = [p["q3"] for p in report.per_position_quality]

        ax.fill_between(positions, q1s, q3s, alpha=0.3, label='IQR')
        ax.plot(positions, means, 'b-', linewidth=2, label='Mean')
        ax.axhline(y=28, color='g', linestyle='--', alpha=0.5, label='Good (Q28)')
        ax.axhline(y=20, color='r', linestyle='--', alpha=0.5, label='Poor (Q20)')
        ax.set_xlabel('Position in read (bp)')
        ax.set_ylabel('Quality score (Phred)')
        ax.set_title(f'Per-position Quality Scores - {report.sample_name}')
        ax.legend()
        ax.set_ylim(0, 42)
        plt.tight_layout()
        plt.savefig(output_dir / 'per_position_quality.png', dpi=150)
        plt.close()

    # 2. Read length distribution
    if report.length_distribution:
        fig, ax = plt.subplots(figsize=(10, 6))
        lengths = sorted(report.length_distribution.keys())
        counts = [report.length_distribution[l] for l in lengths]
        ax.bar(lengths, counts, width=1, edgecolor='none')
        ax.set_xlabel('Read length (bp)')
        ax.set_ylabel('Count')
        ax.set_title(f'Read Length Distribution - {report.sample_name}')
        plt.tight_layout()
        plt.savefig(output_dir / 'length_distribution.png', dpi=150)
        plt.close()

    # 3. GC content distribution
    if report.gc_distribution:
        fig, ax = plt.subplots(figsize=(10, 6))
        gc_bins = sorted(report.gc_distribution.keys())
        counts = [report.gc_distribution[b] for b in gc_bins]
        ax.bar(gc_bins, counts, width=1, edgecolor='none')
        ax.axvline(x=50, color='r', linestyle='--', alpha=0.5, label='50% GC')
        ax.set_xlabel('GC content (%)')
        ax.set_ylabel('Count')
        ax.set_title(f'GC Content Distribution - {report.sample_name}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(output_dir / 'gc_distribution.png', dpi=150)
        plt.close()

    # 4. Duplication levels
    if report.duplication_levels:
        fig, ax = plt.subplots(figsize=(10, 6))
        # Bin high duplication levels
        dup_binned = {}
        for level, count in report.duplication_levels.items():
            if level <= 10:
                dup_binned[str(level)] = count
            elif level <= 50:
                key = ">10"
                dup_binned[key] = dup_binned.get(key, 0) + count
            elif level <= 100:
                key = ">50"
                dup_binned[key] = dup_binned.get(key, 0) + count
            else:
                key = ">100"
                dup_binned[key] = dup_binned.get(key, 0) + count

        labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '>10', '>50', '>100']
        values = [dup_binned.get(l, 0) for l in labels]

        ax.bar(range(len(labels)), values, tick_label=labels)
        ax.set_xlabel('Duplication level')
        ax.set_ylabel('Number of sequences')
        ax.set_title(f'Sequence Duplication Levels - {report.sample_name}')
        plt.tight_layout()
        plt.savefig(output_dir / 'duplication_levels.png', dpi=150)
        plt.close()

    # 5. Quality score distribution
    if report.quality_score_distribution:
        fig, ax = plt.subplots(figsize=(10, 6))
        scores = sorted(report.quality_score_distribution.keys())
        counts = [report.quality_score_distribution[s] for s in scores]
        ax.bar(scores, counts, width=1, edgecolor='none')
        ax.axvline(x=28, color='g', linestyle='--', alpha=0.5, label='Good (Q28)')
        ax.axvline(x=20, color='r', linestyle='--', alpha=0.5, label='Poor (Q20)')
        ax.set_xlabel('Quality score (Phred)')
        ax.set_ylabel('Count')
        ax.set_title(f'Quality Score Distribution - {report.sample_name}')
        ax.legend()
        plt.tight_layout()
        plt.savefig(output_dir / 'quality_distribution.png', dpi=150)
        plt.close()

    print(f"Quality plots saved to {output_dir}")
