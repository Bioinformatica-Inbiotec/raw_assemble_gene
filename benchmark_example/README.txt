# Benchmark Example: `assemble_gene` vs SPAdes

This folder provides a reproducible example to benchmark the performance of the `assemble_gene` tool against the SPAdes assembler, using simulated reads for a single gene of interest.

The benchmark evaluates:
- Number of contigs generated
- Alignment length and identity to the reference gene
- Execution time (wall-clock)
- Maximum memory usage

---

## Requirements

Before running this benchmark, ensure you have the following installed:

- Python 3.7+
- [SPAdes](https://cab.spbu.ru/software/spades/) (≥ v3.15)
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- `/usr/bin/time` (Linux system utility)
- Python packages:
  - `biopython`
  - `openpyxl`

You can install the required Python packages with:

```bash
pip install biopython openpyxl

python3 benchmark.py \
  --reads1 example/hspb1_simulated1.fq \
  --reads2 example/hspb1_simulated2.fq \
  --query example/hspb1.fasta \
  --extension 100
Output File: benchmark_results.xlsx
This spreadsheet includes the following columns for each tool:

Assembler

Contigs (number of sequences in output)

Alignment Length (vs. reference gene)

Identity (%) (BLAST match percentage)

Elapsed Time (s) (wall-clock execution time)

Max Memory (KB) (peak memory used)

Notes
The file assemble_gene.py must be located in the same folder as benchmark.py or be available in your system PATH.

The output file from assemble_gene must end with putative_gene_with_utrs.fasta, as it is automatically detected by the benchmark script.

Contact
For questions or suggestions, please open an issue on the GitHub repository or contact the authors directly.

---

