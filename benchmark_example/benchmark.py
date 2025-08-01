#!/usr/bin/env python3
import os
import subprocess
import argparse
from datetime import datetime
from openpyxl import Workbook
from Bio import SeqIO
from Bio.Blast import NCBIXML
import glob

def run_and_measure(cmd, log_path):
    full_cmd = f"/usr/bin/time -v -o {log_path} {cmd}"
    print(f"Running: {full_cmd}")
    subprocess.run(full_cmd, shell=True, check=True)

def count_contigs(fasta_file):
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

def run_blast(query, subject, out_xml):
    blast_cmd = f"blastn -query {query} -subject {subject} -outfmt 5 -out {out_xml}"
    subprocess.run(blast_cmd, shell=True, check=True)

def parse_blast(xml_file):
    with open(xml_file) as f:
        blast_record = next(NCBIXML.parse(f), None)
    if blast_record and blast_record.alignments:
        hsp = blast_record.alignments[0].hsps[0]
        return {
            "identity": hsp.identities / hsp.align_length * 100,
            "length": hsp.align_length
        }
    return {"identity": 0, "length": 0}

def write_excel(results, out_path):
    wb = Workbook()
    ws = wb.active
    ws.append(["Assembler", "Contigs", "Alignment Length", "Identity (%)", "Elapsed Time (s)", "Max Memory (KB)"])
    for r in results:
        ws.append([r["name"], r["contigs"], r["length"], r["identity"], r["time"], r["memory"]])
    wb.save(out_path)

def parse_time_log(time_log):
    elapsed, memory = None, None
    with open(time_log) as f:
        for line in f:
            if "Elapsed (wall clock) time" in line:
                parts = line.strip().split(":")
                try:
                    if len(parts) == 3:
                        elapsed = int(parts[1]) * 60 + float(parts[2])
                    else:
                        elapsed = float(parts[-1])
                except:
                    elapsed = None
            elif "Maximum resident set size" in line:
                memory = int(line.strip().split(":")[-1])
    return elapsed, memory

def find_assemble_gene_output():
    matches = glob.glob("*.putative_gene_with_utrs.fasta")
    if not matches:
        raise FileNotFoundError("No se encontró el archivo *.putative_gene_with_utrs.fasta generado por assemble_gene.")
    return matches[0]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads1", required=True)
    parser.add_argument("--reads2", required=False)
    parser.add_argument("--query", required=True, help="FASTA file with gene of interest")
    parser.add_argument("--extension", required=True, type=int)
    args = parser.parse_args()

    os.makedirs("benchmark_outputs", exist_ok=True)
    results = []

    # assemble_gene
    assemble_cmd = (
        f"python3 assemble_gene.py "
        f"--reads1 {args.reads1} "
        f"{'--reads2 ' + args.reads2 if args.reads2 else ''} "
        f"--query {args.query} "
        f"--output_prefix benchmark_outputs/assemble_gene "
        f"--extension {args.extension}"
    )
    run_and_measure(assemble_cmd, "benchmark_outputs/assemble_gene_time.log")
    assembled_file = find_assemble_gene_output()
    contigs = count_contigs(assembled_file)
    run_blast(args.query, assembled_file, "benchmark_outputs/assemble_gene_blast.xml")
    blast_result = parse_blast("benchmark_outputs/assemble_gene_blast.xml")
    time_sec, memory_kb = parse_time_log("benchmark_outputs/assemble_gene_time.log")
    results.append({
        "name": "assemble_gene",
        "contigs": contigs,
        "length": blast_result["length"],
        "identity": round(blast_result["identity"], 2),
        "time": round(time_sec, 2),
        "memory": memory_kb
    })

    # SPAdes
    os.makedirs("benchmark_outputs/spades_out", exist_ok=True)
    spades_cmd = (
        f"spades.py "
        f"-1 {args.reads1} "
        f"{'-2 ' + args.reads2 if args.reads2 else ''} "
        f"-o benchmark_outputs/spades_out "
        f"--only-assembler"
    )
    run_and_measure(spades_cmd, "benchmark_outputs/spades_time.log")
    contigs = count_contigs("benchmark_outputs/spades_out/contigs.fasta")
    run_blast(args.query, "benchmark_outputs/spades_out/contigs.fasta", "benchmark_outputs/spades_blast.xml")
    blast_result = parse_blast("benchmark_outputs/spades_blast.xml")
    time_sec, memory_kb = parse_time_log("benchmark_outputs/spades_time.log")
    results.append({
        "name": "SPAdes",
        "contigs": contigs,
        "length": blast_result["length"],
        "identity": round(blast_result["identity"], 2),
        "time": round(time_sec, 2),
        "memory": memory_kb
    })

    write_excel(results, "benchmark_outputs/benchmark_results.xlsx")
    print("Benchmark completo. Resultados guardados en 'benchmark_outputs/benchmark_results.xlsx'")

if __name__ == "__main__":
    main()




# Comando principal:
    # python3 benchmark.py \
    # --reads1 example/hspb1_simulated1.fq \
    # --reads2 example/hspb1_simulated2.fq \
    # --query example/hspb1.fasta \
    # --extension 100
    # --output_excel benchmark_results.xlsx
