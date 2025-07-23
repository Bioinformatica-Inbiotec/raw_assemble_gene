# raw_assemble_gene
Assemble a target gene from raw sequencing reads without full genome assembly

**Assemble_Gene** is a bioinformatics tool to recover a gene of interest from unassembled Illumina reads.  
It aligns reads to a reference gene, filters them, assembles the gene region, and identifies the best contig.

## Installation
We recommend using conda with Bioconda (coming soon). For now, manually install dependencies:

```bash
conda install -c bioconda bwa samtools blast spades biopython
```

## Usage

```bash
python assemble_gene.py --reads example/reads.fastq --reference example/ref_gene.fasta --output putative_gene.fasta
```

## Output
A FASTA file with the best assembled contig corresponding to the gene of interest.

## Example
Included is a test dataset in `example/` with simulated Illumina reads and a reference gene.
