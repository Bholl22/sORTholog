# Rule to make fasta from table

##########################################################################
##########################################################################


rule make_fasta:
    input:
        protein_fasta=merge_db,
        list_all_prot=tsv_prot,
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "reduce_taxid",
            f"all_proteins_reduced.fasta",
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "format_table", "make_fasta.log"),
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/hit2fasta.py"


##########################################################################
##########################################################################


rule make_seed_psiblast:
    input:
        fasta_seed=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.fasta"),
    output:
        fasta=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.psiblast.fasta"),
    params:
        list_proteinId=PSIBLAST,
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "fetch_proteins", "make_seed_psiblast.log"),
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/reduce_seeds.py"


##########################################################################
##########################################################################


rule extract_protein:
    input:
        PAtab=os.path.join(OUTPUT_FOLDER, "results", "patab_melt.tsv"),
        fasta_protein=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
        ),
    output:
        expand(
            os.path.join(OUTPUT_FOLDER, "results", "fasta", "{seed}.fasta"),
            seed=seed_table.seed,
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fetch_proteins",
            "extract_protein.log",
        ),
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/extract_protein.py"


##########################################################################
##########################################################################


rule merge_fasta:
    input:
        taxid_fasta=speedup,
        seed_fasta=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.fasta"),
    output:
        fasta_for_blast=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fetch_proteins",
            "extract_protein.log",
        ),
    shell:
        """cat '{input.seed_fasta}' '{input.taxid_fasta}' > '{output.fasta_for_blast}' 2> {log}"""


##########################################################################
##########################################################################


rule merge_databases:
    input:
        all_databases=list_starting_database,
        annotations=annotationTable,
    output:
        taxid_db=temp(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "merge_databases",
                "databases_all_together.fasta",
            )
        ),
        tsv=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_databases", "protein_table.merged.tsv"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fetch_proteins",
            "merge_databases.log",
        ),
    conda:
        "../envs/biopython_ete3.yaml"
    script:
        "../scripts/merge_protein.py"


##########################################################################
##########################################################################
