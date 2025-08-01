# Example: Gene assembly from simulated reads

This folder contains a reproducible test case for the `assemble_gene` tool.

## Files

- `hspb1_with_utrs.fasta`: Target gene including 200 bp upstream and downstream extensions.
- `hspb1_simulated1.fq`, `hspb1_simulated2.fq`: Paired-end Illumina reads (150 bp), 50× coverage.
- `expected_output/putative_gene_with_utrs.fasta`: Expected result of the gene assembly.

## How to run

```bash
assemble_gene \
    --r1 hspb1_simulated1.fq \
    --r2 hspb1_simulated2.fq \
    --reference hspb1_with_utrs.fasta \
    --output output_dir \
    --extension
