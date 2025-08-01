#!/usr/bin/python3
import os
import sys
import argparse
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Align

def run_command(command):
    print(f"Running: {command}")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        sys.exit(f"Command failed: {command}")

def get_aligned_reads(query, reads1, reads2, minimap2_args, output_prefix):
    sam_file = f"{output_prefix}.sam"
    command = ["minimap2", "-ax", "sr", *minimap2_args.split(), query, reads1, reads2]
    with open(sam_file, "w") as out:
        subprocess.run(command, stdout=out, check=True)

    fastq_file = f"{output_prefix}_reads.fastq"
    command = ["samtools", "fastq", "-f", "0x2", sam_file]
    with open(fastq_file, "w") as out:
        subprocess.run(command, stdout=out, check=True)

    if not os.path.exists(fastq_file) or os.path.getsize(fastq_file) == 0:
        raise ValueError(f"El archivo {fastq_file} no se generó o está vacío.")

    return fastq_file

def run_spades(aligned_reads, outdir):
    print(f"Archivo de entrada para SPAdes: {aligned_reads}")
    command = f"spades.py -s {aligned_reads} -o {outdir}"
    run_command(command)

    contigs_path = os.path.join(outdir, "contigs.fasta")
    if not os.path.exists(contigs_path):
        raise FileNotFoundError("ERROR: contigs.fasta no fue generado por SPAdes.")
    return contigs_path

def find_best_contig(query_seq, contigs_path):
    from Bio import Align
    from Bio import SeqIO

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    min_score = len(query_seq) * 0.75  # Umbral más flexible

    best_score = 0
    best_contig = None

    for contig in SeqIO.parse(contigs_path, "fasta"):
        score = aligner.score(contig.seq, query_seq.seq)
        if score > best_score:
            best_score = score
            best_contig = contig

    if best_score < min_score:
        print(f"[WARNING] {query_seq.id}: No se encontró un contig con similitud completa. Se utilizará el mejor fragmento parcial (score={best_score:.1f}, requerido>{min_score:.1f}).")
        if best_contig is None:
            raise ValueError(f"{query_seq.id}: No se encontró ningún contig con similitud significativa.")

    return best_contig


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
import subprocess
import os

def extend_contig_with_utrs(best_contig, query_seq, extension_length, spades_dir, reads1, reads2):
    from Bio import Align
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import os
    import subprocess

    contig_seq = str(best_contig.seq)
    gene_seq = str(query_seq.seq)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    alignment = aligner.align(contig_seq, gene_seq)[0]

    if alignment.score == 0 or len(alignment.aligned) == 0:
        raise ValueError("El gen no se alinea con el contig. Posiblemente el ensamblado falló.")

    contig_aligned_regions, _ = alignment.aligned
    start = contig_aligned_regions[0][0]
    end = contig_aligned_regions[-1][1]

    # Calcular las regiones extendidas deseadas
    extended_start = max(0, start - extension_length)
    extended_end = min(len(contig_seq), end + extension_length)
    extended_seq = contig_seq[extended_start:extended_end]

    # Coordenadas relativas dentro de la secuencia extendida
    gene_start = start - extended_start
    gene_end = gene_start + (end - start)

    # Crear secuencia anotada con UTRs en minúscula
    annotated_seq = (
        extended_seq[:gene_start].lower() +
        extended_seq[gene_start:gene_end].upper() +
        extended_seq[gene_end:].lower()
    )

    # Validar si se pudo extender, y si no, avisar
    if gene_start == 0:
        print(f"[WARNING] {query_seq.id}: No se pudo extender el 5'-UTR (no hay suficiente secuencia antes del gen).")
    if gene_end == len(annotated_seq):
        print(f"[WARNING] {query_seq.id}: No se pudo extender el 3'-UTR (no hay suficiente secuencia después del gen).")

    return SeqRecord(Seq(annotated_seq), id=query_seq.id, description="putative_gene_with_utrs")

def analyze_coverage(contigs_fasta, reads1, reads2, outdir):
    print("Analizando cobertura del ensamblado...")

    if not contigs_fasta or not os.path.exists(contigs_fasta):
        raise FileNotFoundError(f"No se encontró el archivo de contigs: {contigs_fasta}")
    if not reads1 or not os.path.exists(reads1):
        raise FileNotFoundError(f"No se encontró el archivo reads1: {reads1}")
    if not reads2 or not os.path.exists(reads2):
        raise FileNotFoundError(f"No se encontró el archivo reads2: {reads2}")

    os.makedirs(outdir, exist_ok=True)
    bam_file = os.path.join(outdir, "mapped.bam")
    sorted_bam = os.path.join(outdir, "mapped.sorted.bam")
    cov_file = os.path.join(outdir, "coverage_summary.tsv")
    plot_file = os.path.join(outdir, "coverage_plot.png")
    stats_file = os.path.join(outdir, "assembly_stats.tsv")

    subprocess.run(["bwa", "index", contigs_fasta], check=True)
    p1 = subprocess.Popen(["bwa", "mem", contigs_fasta, reads1, reads2], stdout=subprocess.PIPE)
    with open(bam_file, "wb") as bam_out:
        p2 = subprocess.Popen(["samtools", "view", "-bS", "-"], stdin=p1.stdout, stdout=bam_out)
        p1.stdout.close()
        p2.communicate()

    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    depth_out = subprocess.check_output(["samtools", "depth", sorted_bam])
    depth_lines = depth_out.decode().strip().split("\n")

    data = {}
    for line in depth_lines:
        contig, pos, depth = line.split("\t")
        data.setdefault(contig, []).append(int(depth))

    summary = []
    for contig, depths in data.items():
        summary.append({
            "contig": contig,
            "mean_coverage": sum(depths) / len(depths),
            "max_coverage": max(depths),
            "min_coverage": min(depths)
        })
        plt.plot(range(len(depths)), depths, label=contig)

    pd.DataFrame(summary).to_csv(cov_file, sep="\t", index=False)

    with open(stats_file, "a") as f:
        f.write("\n# Cobertura por contig\n")
        f.write("contig\tmean_coverage\tmax_coverage\tmin_coverage\n")
        for row in summary:
            f.write(f"{row['contig']}\t{row['mean_coverage']:.2f}\t{row['max_coverage']}\t{row['min_coverage']}\n")

    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.title("Contig Coverage")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()

def assembly_stats(contigs_fasta, outdir):
    print("Calculando estadísticas del ensamblado...")
    lengths = [len(rec.seq) for rec in SeqIO.parse(contigs_fasta, "fasta")]
    if not lengths:
        raise ValueError("No se encontraron contigs en el ensamblado.")

    total_contigs = len(lengths)
    total_length = sum(lengths)
    max_length = max(lengths)
    lengths.sort(reverse=True)
    acc = 0
    half = total_length / 2
    for l in lengths:
        acc += l
        if acc >= half:
            n50 = l
            break

    stats_file = os.path.join(outdir, "assembly_stats.tsv")
    with open(stats_file, "w") as f:
        f.write("total_contigs\ttotal_length\tmax_contig_length\tN50\n")
        f.write(f"{total_contigs}\t{total_length}\t{max_length}\t{n50}\n")

    print(f"Estadísticas escritas en {stats_file}")

def main():
    parser = argparse.ArgumentParser(description="Asamblea y extensión de gen desde lecturas")
    parser.add_argument("--query", required=True, help="Archivo FASTA con genes objetivo")
    parser.add_argument("--reads1", required=True, help="Lecturas forward")
    parser.add_argument("--reads2", required=True, help="Lecturas reverse")
    parser.add_argument("--original_reads1", help="Lecturas forward originales")
    parser.add_argument("--original_reads2", help="Lecturas reverse originales")
    parser.add_argument("--minimap2_args", default="", help="Argumentos extra para minimap2")
    parser.add_argument("--extension", type=int, default=0, help="Tamaño de extensión UTR")
    parser.add_argument("--output_prefix", default="genes", help="Prefijo general para reporte final")
    args = parser.parse_args()

    original_reads1 = args.original_reads1 or args.reads1
    original_reads2 = args.original_reads2 or args.reads2

    status_lines = []

    for query_seq in SeqIO.parse(args.query, "fasta"):
        gene_id = query_seq.id
        try:
            aligned_reads = get_aligned_reads(
                args.query, args.reads1, args.reads2,
                args.minimap2_args, gene_id
            )

            spades_dir = os.path.join(f"{gene_id}.output", "spades_output")
            os.makedirs(spades_dir, exist_ok=True)

            contigs_path = run_spades(aligned_reads, spades_dir)
            best_contig = find_best_contig(query_seq, contigs_path)
            annotated = extend_contig_with_utrs(
                best_contig, query_seq, args.extension,
                spades_dir, original_reads1, original_reads2
            )

            output_fasta = f"{gene_id}.putative_gene_with_utrs.fasta"
            SeqIO.write(annotated, output_fasta, "fasta")
            print(f"[OK] {gene_id} guardado en {output_fasta}")
            status_lines.append(f"{gene_id}\tOK")

            analyze_coverage(contigs_path, original_reads1, original_reads2, spades_dir)
            assembly_stats(contigs_path, spades_dir)

        except Exception as e:
            print(f"[ERROR] {gene_id}: {e}")
            status_lines.append(f"{gene_id}\tERROR\t{e}")

    # Nombre del reporte final general
    with open(f"{args.output_prefix}.status_report.txt", "w") as f:
        f.write("GeneID\tStatus\tDetails\n")
        for line in status_lines:
            f.write(line + "\n")
    print(f"Informe final guardado en: {args.output_prefix}.status_report.txt")

if __name__ == "__main__":
    main()
