# Rules to deal with or create HMM profiles

##########################################################################
##########################################################################


rule hmmsearch:
    input:
        hmm=HMM,
        taxid_db=merge_db,
    output:
        domtblout=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmsearch",
            f"hmmsearch--eval_{e_val_HMM:.0e}.out",
        ),
        seed_hmm=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "hmm_db",
            "merged_HMM.hmm",
        ),
    params:
        hmm_type=hmm_type,
        e_value=e_val_HMM,
    resources:
        cpus=5,
        mem_mb=10000,
        time=7200,
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "hmmer", "hmmsearch.log"),
    threads: 5
    conda:
        "../envs/hmmer.yaml"
    envmodules:
        "hmmer/3.3.2",
    shell:
        """
        cat {input.hmm} > '{output.seed_hmm}' 2> '{log}'

        hmmsearch --domtblout {output.domtblout} {params.hmm_type} {params.e_value} --cpu {threads} \
            {output.seed_hmm} {input.taxid_db} &>> '{log}'
        """


##########################################################################
##########################################################################
