#!/usr/bin/env python3
"""
Command-line interface for phage display NGS analysis pipeline.

Usage:
    python -m phage_ngs_analysis.cli quality input.fastq.gz -o qc_report/
    python -m phage_ngs_analysis.cli map --library library.fasta --samples sample1.fq,sample2.fq
    python -m phage_ngs_analysis.cli analyze-unmapped unmapped.fasta -o unmapped_analysis/
    python -m phage_ngs_analysis.cli pipeline --config config.yaml
"""

import argparse
import sys
from pathlib import Path
import json


def cmd_quality(args):
    """Run quality assessment on FASTQ files."""
    from .quality_assessment import QualityAssessor, generate_quality_plots

    assessor = QualityAssessor(
        upstream_flank=args.upstream_flank,
        downstream_flank=args.downstream_flank,
        sample_size=args.sample_size,
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    for fastq_path in args.input:
        fastq_path = Path(fastq_path)
        print(f"\nAssessing quality: {fastq_path.name}")

        report = assessor.assess(fastq_path)
        print(report.summary())

        # Save report
        report_name = fastq_path.stem.replace('.fastq', '')
        report.save_json(output_dir / f"{report_name}_quality.json")

        # Generate plots if requested
        if args.plots:
            plot_dir = output_dir / f"{report_name}_plots"
            generate_quality_plots(report, plot_dir)


def cmd_map(args):
    """Map reads to library and calculate enrichment."""
    from .library_mapping import LibraryMapper, generate_mapping_plots

    # Initialize mapper
    mapper = LibraryMapper(
        library_path=args.library,
        upstream_flank=args.upstream_flank,
        downstream_flank=args.downstream_flank,
        max_mismatches=args.max_mismatches,
        min_quality=args.min_quality,
        check_reverse_complement=not args.no_reverse_complement,
    )

    # Parse sample files
    sample_paths = {}
    for sample_spec in args.samples:
        if ':' in sample_spec:
            name, path = sample_spec.split(':', 1)
        else:
            path = sample_spec
            name = Path(path).stem.replace('.fastq', '').replace('.fq', '')
        sample_paths[name] = path

    # Run mapping
    print(f"\nMapping {len(sample_paths)} samples to library...")
    result = mapper.map_samples(sample_paths, library_name=args.library_name or "library")

    # Calculate enrichment if input sample specified
    if args.input_sample:
        print(f"\nCalculating enrichment relative to {args.input_sample}...")
        result.calculate_enrichment(args.input_sample)

    # Save results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    result.save_counts(output_dir / "counts.tsv")
    print(f"Saved counts to {output_dir / 'counts.tsv'}")

    if result.enrichment:
        result.save_enrichment(output_dir / "enrichment.tsv")
        print(f"Saved enrichment to {output_dir / 'enrichment.tsv'}")

    # Save unmapped sequences
    for sample in result.samples:
        if result.unmapped_sequences.get(sample):
            result.save_unmapped(output_dir / f"{sample}_unmapped.fasta", sample)

    # Generate plots
    if args.plots:
        generate_mapping_plots(result, output_dir / "plots")

    # Print summary
    print("\n" + result.summary())

    # Print top enriched if calculated
    if result.enrichment and args.input_sample:
        for sample in result.samples:
            if sample == args.input_sample:
                continue
            print(f"\nTop 10 enriched in {sample}:")
            for member_id, log2fc, count in result.get_top_enriched(sample, n=10):
                print(f"  {member_id}: log2FC={log2fc:.2f}, count={count}")


def cmd_analyze_unmapped(args):
    """Analyze unmapped sequences."""
    from .enrichment_analysis import EnrichmentAnalyzer, generate_enrichment_plots
    from .utils import read_fasta, read_fastq

    # Load library if provided
    library_sequences = {}
    if args.library:
        from .utils import read_fasta
        for record in read_fasta(args.library):
            library_sequences[record.id] = record.sequence

    analyzer = EnrichmentAnalyzer(
        library_sequences=library_sequences,
        expected_length=args.expected_length,
        similarity_threshold=args.similarity_threshold,
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    for input_path in args.input:
        input_path = Path(input_path)
        print(f"\nAnalyzing unmapped sequences: {input_path.name}")

        # Load sequences
        sequences = []
        if input_path.suffix in ['.fasta', '.fa', '.fna'] or '.fasta' in str(input_path):
            for record in read_fasta(input_path):
                sequences.append(record.sequence)
        else:
            for record in read_fastq(input_path):
                sequences.append(record.sequence)

        sample_name = input_path.stem.replace('_unmapped', '')
        analysis = analyzer.analyze(sequences, sample_name)

        print(analysis.summary())

        # Save results
        analysis.save_json(output_dir / f"{sample_name}_unmapped_analysis.json")

        # Generate plots
        if args.plots:
            generate_enrichment_plots(analysis, output_dir / f"{sample_name}_plots")


def cmd_pipeline(args):
    """Run complete analysis pipeline from config file."""
    import yaml

    config_path = Path(args.config)
    with open(config_path) as f:
        config = yaml.safe_load(f)

    output_dir = Path(config.get('output_dir', 'phage_ngs_output'))
    output_dir.mkdir(parents=True, exist_ok=True)

    # Quality assessment
    if config.get('quality_assessment', {}).get('enabled', True):
        print("\n" + "="*60)
        print("STEP 1: Quality Assessment")
        print("="*60)

        from .quality_assessment import QualityAssessor, generate_quality_plots

        qc_config = config.get('quality_assessment', {})
        assessor = QualityAssessor(
            upstream_flank=config.get('flanking_sequences', {}).get('upstream'),
            downstream_flank=config.get('flanking_sequences', {}).get('downstream'),
            sample_size=qc_config.get('sample_size'),
        )

        qc_dir = output_dir / "quality"
        qc_dir.mkdir(exist_ok=True)

        for sample in config.get('samples', []):
            fastq_path = sample.get('fastq')
            if fastq_path:
                report = assessor.assess(fastq_path, sample.get('name'))
                print(report.summary())
                report.save_json(qc_dir / f"{sample.get('name', 'sample')}_quality.json")
                if qc_config.get('generate_plots', True):
                    generate_quality_plots(report, qc_dir / f"{sample.get('name', 'sample')}_plots")

    # Library mapping
    if config.get('library_mapping', {}).get('enabled', True):
        print("\n" + "="*60)
        print("STEP 2: Library Mapping")
        print("="*60)

        from .library_mapping import LibraryMapper, generate_mapping_plots

        map_config = config.get('library_mapping', {})
        mapper = LibraryMapper(
            library_path=config.get('library'),
            upstream_flank=config.get('flanking_sequences', {}).get('upstream'),
            downstream_flank=config.get('flanking_sequences', {}).get('downstream'),
            max_mismatches=map_config.get('max_mismatches', 0),
            min_quality=map_config.get('min_quality', 20),
        )

        sample_paths = {}
        for sample in config.get('samples', []):
            sample_paths[sample.get('name')] = sample.get('fastq')

        result = mapper.map_samples(sample_paths)

        # Calculate enrichment
        input_sample = config.get('input_sample')
        if input_sample:
            result.calculate_enrichment(input_sample)

        # Save results
        map_dir = output_dir / "mapping"
        map_dir.mkdir(exist_ok=True)

        result.save_counts(map_dir / "counts.tsv")
        if result.enrichment:
            result.save_enrichment(map_dir / "enrichment.tsv")

        for sample in result.samples:
            if result.unmapped_sequences.get(sample):
                result.save_unmapped(map_dir / f"{sample}_unmapped.fasta", sample)

        if map_config.get('generate_plots', True):
            generate_mapping_plots(result, map_dir / "plots")

        print(result.summary())

    # Unmapped sequence analysis
    if config.get('unmapped_analysis', {}).get('enabled', True):
        print("\n" + "="*60)
        print("STEP 3: Unmapped Sequence Analysis")
        print("="*60)

        from .enrichment_analysis import EnrichmentAnalyzer, generate_enrichment_plots
        from .utils import read_fasta

        unmapped_config = config.get('unmapped_analysis', {})

        # Load library for variant detection
        library_sequences = {}
        if config.get('library'):
            for record in read_fasta(config.get('library')):
                library_sequences[record.id] = record.sequence

        analyzer = EnrichmentAnalyzer(
            library_sequences=library_sequences,
            expected_length=unmapped_config.get('expected_length'),
            similarity_threshold=unmapped_config.get('similarity_threshold', 0.9),
        )

        unmapped_dir = output_dir / "unmapped_analysis"
        unmapped_dir.mkdir(exist_ok=True)

        # Analyze unmapped sequences from mapping step
        map_dir = output_dir / "mapping"
        for unmapped_file in map_dir.glob("*_unmapped.fasta"):
            sequences = [r.sequence for r in read_fasta(unmapped_file)]
            if sequences:
                sample_name = unmapped_file.stem.replace('_unmapped', '')
                analysis = analyzer.analyze(sequences, sample_name)
                print(analysis.summary())
                analysis.save_json(unmapped_dir / f"{sample_name}_analysis.json")
                if unmapped_config.get('generate_plots', True):
                    generate_enrichment_plots(analysis, unmapped_dir / f"{sample_name}_plots")

    print("\n" + "="*60)
    print(f"Pipeline complete! Results saved to: {output_dir}")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="Phage Display NGS Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quality assessment
  python -m phage_ngs_analysis.cli quality sample.fastq.gz -o qc_output/

  # Map reads to library
  python -m phage_ngs_analysis.cli map --library lib.fasta \\
      --samples input:input.fq,round1:r1.fq,round2:r2.fq \\
      --input-sample input -o mapping_output/

  # Analyze unmapped sequences
  python -m phage_ngs_analysis.cli analyze-unmapped unmapped.fasta \\
      --library lib.fasta -o unmapped_output/

  # Run complete pipeline
  python -m phage_ngs_analysis.cli pipeline --config analysis_config.yaml
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Quality assessment command
    qc_parser = subparsers.add_parser('quality', help='Assess read quality')
    qc_parser.add_argument('input', nargs='+', help='Input FASTQ file(s)')
    qc_parser.add_argument('-o', '--output', default='qc_output', help='Output directory')
    qc_parser.add_argument('--upstream-flank', help='Upstream flanking sequence')
    qc_parser.add_argument('--downstream-flank', help='Downstream flanking sequence')
    qc_parser.add_argument('--sample-size', type=int, help='Subsample reads for faster analysis')
    qc_parser.add_argument('--plots', action='store_true', help='Generate visualization plots')

    # Mapping command
    map_parser = subparsers.add_parser('map', help='Map reads to library')
    map_parser.add_argument('--library', required=True, help='Library FASTA file')
    map_parser.add_argument('--samples', required=True, help='Sample files (name:path,name:path,...)')
    map_parser.add_argument('-o', '--output', default='mapping_output', help='Output directory')
    map_parser.add_argument('--library-name', help='Name for this library')
    map_parser.add_argument('--upstream-flank', help='Upstream flanking sequence')
    map_parser.add_argument('--downstream-flank', help='Downstream flanking sequence')
    map_parser.add_argument('--max-mismatches', type=int, default=0, help='Max mismatches for mapping')
    map_parser.add_argument('--min-quality', type=int, default=20, help='Min average quality score')
    map_parser.add_argument('--input-sample', help='Input sample name for enrichment calculation')
    map_parser.add_argument('--no-reverse-complement', action='store_true', help='Disable reverse complement checking')
    map_parser.add_argument('--plots', action='store_true', help='Generate visualization plots')

    # Unmapped analysis command
    unmapped_parser = subparsers.add_parser('analyze-unmapped', help='Analyze unmapped sequences')
    unmapped_parser.add_argument('input', nargs='+', help='Input FASTA/FASTQ file(s) with unmapped sequences')
    unmapped_parser.add_argument('-o', '--output', default='unmapped_output', help='Output directory')
    unmapped_parser.add_argument('--library', help='Library FASTA for variant detection')
    unmapped_parser.add_argument('--expected-length', type=int, help='Expected sequence length')
    unmapped_parser.add_argument('--similarity-threshold', type=float, default=0.9, help='Clustering similarity threshold')
    unmapped_parser.add_argument('--plots', action='store_true', help='Generate visualization plots')

    # Pipeline command
    pipeline_parser = subparsers.add_parser('pipeline', help='Run complete analysis pipeline')
    pipeline_parser.add_argument('--config', required=True, help='Pipeline configuration YAML file')

    args = parser.parse_args()

    if args.command == 'quality':
        cmd_quality(args)
    elif args.command == 'map':
        cmd_map(args)
    elif args.command == 'analyze-unmapped':
        cmd_analyze_unmapped(args)
    elif args.command == 'pipeline':
        cmd_pipeline(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
