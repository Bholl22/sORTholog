from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Read the seed file and create Index database to better parse
parser = SeqIO.parse(snakemake.input.fasta_seed, "fasta")

with open(snakemake.output.fasta, "wt") as w_file:
    for protein in parser:
        for protein_id in snakemake.params.list_proteinId:
            if protein_id in protein.id:
                SeqIO.write(protein, w_file, "fasta")
