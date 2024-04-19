import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Opening blastout
blast_names = [
    "qacc",
    "qlen",
    "qseq",
    "qstart",
    "qend",
    "sacc",
    "slen",
    "sseq",
    "sstart",
    "send",
    "length",
    "pident",
    "evalue",
    "bitscore",
    "qcovs",
    "qcovhsp",
    "ssciname",
    "sblastname",
    "stitle",
]

all_sacc = []

# Reading the file line by line
with open(snakemake.input.psiblast, "rt") as r_file:
    for line in r_file:
        # Cleaning blastout
        if (
            "Search has CONVERGED!" not in line
            and not line.startswith("#")
            and not line.startswith("\n")
        ):
            line_split = line.split()
            sacc = line_split[5]

            # Getting the list of protein matches unique
            if sacc not in all_sacc:
                all_sacc.append(sacc)

with open(snakemake.output.list_all_prot, "w") as w_file:
    w_file.write("protein_id\n")

    for sacc in all_sacc:
        w_file.write(f"{sacc}\n")
