from Bio import SeqIO
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Get all the proteins of the database inside a dict-like structure
all_index_fasta = SeqIO.index(snakemake.input.protein_fasta, "fasta")

# List of possible truncatenated ids
list_cutname = []

# Read line by line the protein table hit without loading it in memory
with open(snakemake.input.list_all_prot, "rt") as r_file:
    with open(snakemake.output.fasta, "wt") as w_file:
        for line in r_file:
            if not line.startswith("protein_id"):
                tmp_line = line.split()
                proteins_of_interest = tmp_line[0]

                if proteins_of_interest in all_index_fasta:
                    SeqIO.write(all_index_fasta[proteins_of_interest], w_file, "fasta")
                else:
                    list_cutname.append(proteins_of_interest)

        # Now parse the file to get the trucatenate fasta sequences
        parser = SeqIO.parse(snakemake.input.protein_fasta, "fasta")

        for prot in parser:
            for cutname in list_cutname:
                if prot.id.startswith(cutname):
                    SeqIO.write(prot, w_file, "fasta")
