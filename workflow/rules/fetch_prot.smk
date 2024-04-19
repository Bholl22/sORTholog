# Rules to fetch proteins fasta files from NCBI

##########################################################################
##########################################################################


rule fetch_proteins_database:
    input:
        taxid,
    output:
        fasta_final=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "taxid_all_together.fasta"
        ),
        assembly_output=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "summary_assembly_taxid.tsv"
        ),
        new_taxid_file=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "taxid_checked.txt"
        ),
        output_table_protein=os.path.join(
            OUTPUT_FOLDER, "databases", "all_taxid", "protein_table.tsv"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER, "logs", "fetch_proteins", "fetch_proteins_database.log"
        ),
    params:
        output_database_folder=lambda x, output: os.path.dirname(
            output.output_table_protein
        ),
        section=section,
        assembly_levels=assembly_levels,
        refseq_categories=refseq_categories,
        update_db=config["update_db"],
    resources:
        cpus=5,
        time=3000,
        mem_mb=20000,
    threads: 5
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/fetch_proteome.py"


##########################################################################
##########################################################################


rule fetch_fasta_from_seed:
    input:
        seed_file=start_seed_file,
    output:
        fasta_seed=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.fasta"),
        new_seed_file=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "new_seeds.tsv"),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fetch_proteins",
            "fetch_fasta_from_seed.log",
        ),
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/fetch_seed.py"


##########################################################################
##########################################################################
