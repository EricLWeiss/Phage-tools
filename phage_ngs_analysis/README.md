# Phage Display NGS Analysis Pipeline

A Python toolkit for analyzing Next-Generation Sequencing (NGS) data from phage display experiments. This pipeline helps you understand what sequences are enriched during biopanning selection rounds.

## What This Pipeline Does

1. **Assesses Read Quality** - Checks if your sequencing data is good enough for analysis (similar to FastQC)
2. **Maps Reads to Your Library** - Counts how many times each library member appears in each sample
3. **Calculates Enrichment** - Identifies which sequences became more abundant after selection
4. **Analyzes Unmapped Sequences** - Investigates sequences that don't match your input library (mutations, artifacts, contamination)

## Installation

### Step 1: Install Python

If you don't have Python installed:
- **Mac**: Download from [python.org](https://www.python.org/downloads/) or install via Homebrew: `brew install python`
- **Windows**: Download from [python.org](https://www.python.org/downloads/) and run the installer (check "Add Python to PATH")
- **Linux**: Usually pre-installed. If not: `sudo apt install python3 python3-pip`

### Step 2: Download This Package

Download or clone this repository to your computer:
```bash
git clone https://github.com/EricLWeiss/Phage-tools.git
cd Phage-tools
```

### Step 3: Install Dependencies

Open a terminal/command prompt in the Phage-tools folder and run:
```bash
pip install -r requirements.txt
```

## Preparing Your Data

### What You Need

1. **FASTQ files** - Your sequencing data (can be gzip compressed: `.fastq.gz`)
2. **Library FASTA file** - A file containing all the sequences in your input library

### FASTQ Files

Your sequencing facility provides these. They contain the actual sequence reads. Example format:
```
@READ_001
GCTAGCATGCATGCATGCATGCGGATCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
```

### Library FASTA File

Create a text file listing all your library member sequences:
```
>member_001
ATGCATGCATGCATGC
>member_002
GCTAGCTAGCTAGCTA
>member_003
TATATATATATATATA
```

### Flanking Sequences (Important!)

If your reads contain constant regions flanking the variable/randomized region, identify these sequences. For example, if your construct is:

```
5'--[constant]--[VARIABLE REGION]--[constant]--3'
      GCTAGC     NNNNNNNNNNNNNN     GGATCC
```

You'll need:
- **Upstream flank**: `GCTAGC` (sequence just before the variable region)
- **Downstream flank**: `GGATCC` (sequence just after the variable region)

## Running the Analysis

### Option A: Step-by-Step Analysis

#### 1. Quality Assessment

Check if your sequencing data is good quality:

```bash
python -m phage_ngs_analysis.cli quality \
    your_sample.fastq.gz \
    -o quality_results/ \
    --plots
```

**What to look for in results:**
- Mean quality score should be >28 (good) or >20 (acceptable)
- Check for quality drop-off at read ends
- High duplication is normal for phage display (you're selecting for specific sequences!)

#### 2. Map Reads to Library

Count how many reads match each library member:

```bash
python -m phage_ngs_analysis.cli map \
    --library your_library.fasta \
    --samples naive:naive.fastq.gz,round1:round1.fastq.gz,round2:round2.fastq.gz \
    --input-sample naive \
    --upstream-flank GCTAGC \
    --downstream-flank GGATCC \
    --max-mismatches 1 \
    -o mapping_results/ \
    --plots
```

**Parameters explained:**
- `--library`: Your library FASTA file
- `--samples`: Your FASTQ files in format `name:path,name:path,...`
- `--input-sample`: Which sample is your naive/input library (for enrichment calculation)
- `--upstream-flank` / `--downstream-flank`: Constant sequences flanking your variable region
- `--max-mismatches`: Allow fuzzy matching (1-2 mismatches accounts for sequencing errors)

#### 3. Analyze Unmapped Sequences

Understand sequences that don't match your library:

```bash
python -m phage_ngs_analysis.cli analyze-unmapped \
    mapping_results/round2_unmapped.fasta \
    --library your_library.fasta \
    -o unmapped_results/ \
    --plots
```

### Option B: Run Everything at Once

Create a configuration file (copy and edit `phage_ngs_analysis/example_config.yaml`):

```yaml
output_dir: my_experiment_results

library: path/to/my_library.fasta

flanking_sequences:
  upstream: GCTAGC
  downstream: GGATCC

input_sample: naive

samples:
  - name: naive
    fastq: data/naive_library.fastq.gz
  - name: round1
    fastq: data/selection_round1.fastq.gz
  - name: round2
    fastq: data/selection_round2.fastq.gz
  - name: round3
    fastq: data/selection_round3.fastq.gz

quality_assessment:
  enabled: true
  generate_plots: true

library_mapping:
  enabled: true
  max_mismatches: 1
  generate_plots: true

unmapped_analysis:
  enabled: true
  generate_plots: true
```

Then run:
```bash
python -m phage_ngs_analysis.cli pipeline --config my_config.yaml
```

## Understanding Your Results

### Output Files

```
my_experiment_results/
├── quality/
│   ├── naive_quality.json          # Quality metrics
│   └── naive_plots/                # Quality visualization plots
│       ├── per_position_quality.png
│       ├── length_distribution.png
│       └── gc_distribution.png
│
├── mapping/
│   ├── counts.tsv                  # Read counts per library member
│   ├── enrichment.tsv              # Log2 fold-change values
│   ├── round1_unmapped.fasta       # Sequences that didn't map
│   └── plots/
│       ├── mapping_rates.png
│       ├── abundance_heatmap.png
│       └── enrichment_round2.png
│
└── unmapped_analysis/
    ├── round2_analysis.json        # Analysis of unmapped sequences
    └── round2_plots/
        └── unmapped_sources.png    # Pie chart of likely sources
```

### Key Output Files Explained

#### `counts.tsv` - Raw Read Counts
A table showing how many reads matched each library member in each sample:

| library_member | naive | round1 | round2 | round3 |
|----------------|-------|--------|--------|--------|
| member_001     | 150   | 890    | 4521   | 12340  |
| member_002     | 142   | 95     | 12     | 0      |
| member_003     | 148   | 2100   | 8900   | 25000  |

#### `enrichment.tsv` - Fold-Change Values
Log2 fold-change relative to your input sample. **Positive values = enriched, Negative = depleted**:

| library_member | log2FC_round1 | log2FC_round2 | log2FC_round3 |
|----------------|---------------|---------------|---------------|
| member_001     | 2.5           | 4.8           | 6.2           |
| member_002     | -0.6          | -3.2          | -inf          |
| member_003     | 3.8           | 5.9           | 7.4           |

**Interpreting fold-change:**
- `log2FC = 2` means 4x more abundant (2² = 4)
- `log2FC = 3` means 8x more abundant (2³ = 8)
- `log2FC = -2` means 4x less abundant
- Large positive values = strong enrichment (good binders!)

### What Good Results Look Like

1. **Quality Assessment**
   - Mean quality >28 (green) or >20 (acceptable)
   - Most reads at expected length
   - GC content appropriate for your library design

2. **Mapping**
   - Mapping rate >50% (depends on library complexity and flanking sequence accuracy)
   - Clear enrichment pattern across selection rounds
   - Top enriched sequences show progressively higher counts

3. **Unmapped Analysis**
   - Most unmapped sequences are "library mutations" (1-2 base errors)
   - Low contamination percentage
   - Novel sequences may indicate interesting variants or recombination

### Troubleshooting

| Problem | Possible Cause | Solution |
|---------|----------------|----------|
| Low mapping rate (<20%) | Wrong flanking sequences | Check your construct design, try without flanking sequences |
| No enrichment observed | Selection didn't work | Check your biopanning protocol |
| High "contamination" in unmapped | Adapter sequences present | Trim adapters before analysis |
| All reads same quality score | Quality scores not in Phred+33 format | Check FASTQ format with sequencing facility |

## Getting Help

- Check the example configuration file: `phage_ngs_analysis/example_config.yaml`
- For bugs or feature requests: [GitHub Issues](https://github.com/EricLWeiss/Phage-tools/issues)

## Citation

If you use this tool in your research, please cite this repository.
