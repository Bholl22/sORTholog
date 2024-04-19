# Module containing all the ncbi-blast related rules

##########################################################################
##########################################################################


rule psiblast:
    input:
        seed=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "seeds.psiblast.fasta"),
        taxid_db=merge_db,
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "psiblast",
            f"psiblast--eval_{e_val_psiblast:.0e}_raw.out",
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "psiblast", "psiblast.log"),
    params:
        e_value=e_val_psiblast,
        iteration=iteration_psiblast,
    resources:
        cpus=5,
        mem_mb=20000,
        time=7200,
    conda:
        "../envs/blast.yaml"
    envmodules:
        "ncbi_blast/2.10.1",
    threads: 5
    shell:
        """
        makeblastdb -dbtype prot -in '{input.taxid_db}' &> {log}

        psiblast -query '{input.seed}' -db '{input.taxid_db}' -evalue {params.e_value} \
                 -outfmt '7 qacc qlen qseq qstart qend sacc slen sseq sstart send length pident evalue bitscore qcovs' \
                 -num_threads {threads} -num_iterations {params.iteration} -out '{output}' \
                 -num_alignments 250000 2>> {log}

        rm {input.taxid_db}.*
        """


##########################################################################
##########################################################################


# rule blast:
#     input:
#         fasta_for_blast=os.path.join(
#             OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
#         ),
#     output:
#         blast_out=os.path.join(
#             OUTPUT_FOLDER, "processing_files", "blast", "blastp--blast_evalue_1e-2.out"
#         ),
#     log:
#         os.path.join(OUTPUT_FOLDER, "logs", "blast", "blast.log"),
#     resources:
#         cpus=5,
#         mem_mb=20000,
#         time=43200,
#     conda:
#         "../envs/blast.yaml"
#     envmodules:
#         "ncbi_blast/2.10.1",
#     threads: 5
#     shell:
#         """
#         makeblastdb -dbtype prot -in '{input}' &> {log}

#         blastp -query '{input}' -db '{input}' -evalue 0.01 \
#                -outfmt 6 -out '{output.blast_out}' -num_threads {threads} -num_alignments 25000 2>> {log}

#         rm '{input}'.*
#         """

##########################################################################
##########################################################################


rule blast:
    input:
        fasta_for_blast=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER, "processing_files", "blast", "blastp--blast_evalue_1e-2.out"
        ),
    params:
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processing_files", "blast", "tmp"
        ),        
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "blast", "blast.log"),
    resources:
        cpus=5,
        mem_mb=20000,
        time=43200,
    conda:
        "../envs/blast.yaml"
    envmodules:
        "ncbi_blast/2.10.1",
    threads: 5
    script:
        "../scripts/blastp_wrapper.py"


##########################################################################
##########################################################################
