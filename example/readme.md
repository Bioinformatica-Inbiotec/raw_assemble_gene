# Example: Assembling the hspb1 gene with simulated Illumina reads

This folder contains a complete example to test the `assemble_gene` tool using simulated paired-end Illumina reads.

## Files included

- `hspb1_with_utrs.fasta`: input gene sequence including 200 bp upstream and downstream (putative UTRs).
- `hspb1_simulated1.fq`: simulated R1 reads (150 bp, 50× coverage).
- `hspb1_simulated2.fq`: simulated R2 reads (paired-end).
- `assemble_gene.py`: the assembly script (minimal version).
- `expected_output/hspb1.putative_gene_with_utrs.fasta`: the expected final output after running the tool.

## How to run this example

Make sure you have all dependencies installed (e.g., `spades.py`, `minimap2`, `samtools`, `bedtools`, and Python packages like `biopython`).

From inside the `example/` directory, run:

```bash
python3 assemble_gene.py \
  --query hspb1_with_utrs.fasta \
  --reads1 hspb1_simulated1.fq \
  --reads2 hspb1_simulated2.fq \
  --extension 200 \
  --output_prefix hspb1_example
