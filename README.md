# HiFiSim

**HiFiSim** is a Python-based command-line tool for simulating synthetic PacBio HiFi reads.  
It introduces sequencing errors (insertions, deletions, substitutions), supports different read length distributions, and generates FASTA or FASTQ files with optional quality scores. This tool is useful for testing bioinformatics pipelines where real data may be too sensitive, expensive, or incomplete.

---

## Requirements

- Python 3.7+
- Standard libraries only: `argparse`, `os`, `random`, `sys`

Project files:
- `HiFiSimCLI.py` – Command-line interface
- `FileGlance.py` – Glance at FASTA or FASTQ files
- `ReadAnalysis.py` – Analyzes statistics of given reads
- `SimReads.py` – Generates synthetic reads
- `GenerateQscores.py` – Models quality score profiles

---

## Quick Start

Clone the repo:

```bash
git clone https://github.com/your-username/hifireadsim.git
cd hifireadsim
```

Run the CLI tool:

```bash
python HiFiSimCLI.py --help
```

---

## Usage Example

```bash
python HiFiSimCLI.py \
  --input reference_genome.fa \
  --outputname test_reads \
  --outputformat fastq \
  --numreads 1000 \
  --errorrate 0.02 \
  --insertf 0.2 \
  --deletef 0.2 \
  --subf 0.6 \
  --readlengthmode lognormal \
  --lengthmean 3000 \
  --lengthsd 500 \
  --qscoremode normal \
  --qscoreparams 25 5
```

---

## CLI Arguments

| Argument             | Required | Description |
|----------------------|----------|-------------|
| `--input`            | ✅        | Reference genome file or folder |
| `--outputname`       | ✅        | Prefix for output file |
| `--outputformat`     | ✅        | Output format: `fasta` or `fastq` |
| `--numreads`         | ❌        | Number of reads to simulate (auto-selected if not provided) |
| `--errorrate`        | ❌        | Overall error rate (0.0–1.0) |
| `--insertf`          | ❌        | Fraction of insertions |
| `--deletef`          | ❌        | Fraction of deletions |
| `--subf`             | ❌        | Fraction of substitutions |
| `--readlengthmode`   | ❌        | Choose `lognormal` or `empirical` |
| `--lengthmean`       | ❌        | Mean read length (for lognormal) |
| `--lengthsd`         | ❌        | Standard deviation of read length (for lognormal) |
| `--qscoremode`       | ❌        | Quality score model: `normal` or `weighted` |
| `--qscoreparams`     | ❌        | Parameters depending on mode:<br>- `normal`: mean sd<br>- `weighted`: score1 weight1 score2 weight2 ... |

---

## Output

- `test_reads.fasta` or `test_reads.fastq`: Synthetic HiFi reads.
- Each read contains random variation and modeled sequencing errors.

---

## Applications

- Testing genome assembly pipelines
- Evaluating contamination handling
- Metagenomic simulations
- Benchmarking variant callers (SNPs/indels)
- Resequencing with reference genomes

---

## Known Limitations

- Single-end reads only (no paired-end support)
- No structural variation simulation
- Output is stochastic unless seeded
- No built-in visualization or coverage control (yet)

---

## Planned Features

- Paired-end read simulation
- Improved modeling of real error distributions
- Sequence-specific quality decay
- Built-in visualization and summary reports

---

## Author

Created by **Cianna Grama**  
GitHub: [@cianna-grama](https://github.com/cianna-grama)

---

## License

MIT License – free for personal or academic use.

