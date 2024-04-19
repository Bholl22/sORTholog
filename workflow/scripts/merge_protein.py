import pandas as pd
import sys
from Bio import SeqIO

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

num_files = len(snakemake.input.annotations)

snakemake.output.taxid_db

annotations_id_dict = {}
annotations_name_dict = {}

# In annotations[0] if could be an empty string or the table provided by the user
if snakemake.input.annotations[0]:
    annotations_df = pd.read_table(snakemake.input.annotations[0])
    annotations_id_dict = annotations_df.set_index("protein_id").genome_id.to_dict()

    if "genome_name" in annotations_df.columns:
        annotations_name_dict = annotations_df.set_index(
            "protein_id"
        ).genome_name.to_dict()

# First there is always all_databases[0] because there is at least the perso database
parser_fasta = SeqIO.parse(snakemake.input.all_databases[0], "fasta")

tmp_df = []

with open(snakemake.output.taxid_db, "wt") as fasta_file:
    for prot in parser_fasta:
        length = len(prot.seq)

        # Security to remove protein id from description
        description = prot.description.replace(f"{prot.id} ", "")
        description = prot.description.replace(f"{prot.id}", "")

        if annotations_name_dict:
            genome_name = annotations_name_dict[prot.id]
        else:
            # if no genome name create a fake one
            genome_name = ""

        if annotations_id_dict:
            genome_id = annotations_id_dict[prot.id]
            protein_id = f"{prot.id}--{genome_id}"
        else:
            id_split = prot.id.split("--")
            genome_id = id_split[1]
            protein_id = prot.id

        # Adding the information inside the datafame to be
        tmp_df.append(
            {
                "protein_id": protein_id,
                "genome_id": genome_id,
                "length": length,
                "protein_description": description,
                "genome_name": genome_name,
            }
        )

        prot.id = protein_id
        prot.description = description
        prot.name = ""

        SeqIO.write(prot, fasta_file, "fasta")

    merge_df = pd.DataFrame(tmp_df)

    if num_files == 2:
        prot_taxid = snakemake.input.annotations[1]
        taxid_df = pd.read_table(prot_taxid)

        merge_df = pd.concat([taxid_df, merge_df])

        parser = SeqIO.parse(snakemake.input.all_databases[1], "fasta")

        for prot in parser:
            SeqIO.write(prot, fasta_file, "fasta")

merge_df.to_csv(snakemake.output.tsv, index=False, sep="\t")
