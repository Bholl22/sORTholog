from Bio import SeqIO
import sys, os

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

index_fasta = SeqIO.index(snakemake.input.fasta_protein, "fasta")

open_output = [open(output, "wt") for output in snakemake.output]
list_seed = [os.path.basename(fasta).split(".")[0] for fasta in snakemake.output]

with open(snakemake.input.PAtab, "rt") as r_file:
    header = r_file.readline().split()
    index_protein = header.index("protein_id")
    index_genome = header.index("genome_id")
    index_seed = header.index("seed")
    index_PA = header.index("PA")

    for line in r_file:
        split_line = line.rstrip().split("\t")

        if split_line[index_PA] != "0":
            protein_ids = split_line[index_protein].split(",")
            index_file = list_seed.index(split_line[index_seed])

            for protein_id in protein_ids:
                tmp_protein = f"{protein_id}--{split_line[index_genome]}"
                protein = index_fasta[tmp_protein]
                protein.id = protein_id
                protein.name = ""
                protein.description = protein.description.split(f"{protein.id} ")[-1]

                SeqIO.write(protein, open_output[index_file], "fasta")

for open_file in open_output:
    open_file.close()
