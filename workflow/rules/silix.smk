# Rule for silix

##########################################################################
##########################################################################


rule silix:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
        ),
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "blast",
            "split_blast_out",
            "filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.out",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "silix",
            "fnodes_files",
            "filtered_blast--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.fnodes",
        ),
    params:
        minimum_overlap=silix_dict[cov_min],
        minimum_percId=silix_dict[pid_min],
        minimum_length=length_min,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "silix",
            "{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.silix.log",
        ),
    conda:
        "../envs/silix.yaml"
    envmodules:
        "silix/1.2.11",
    shell:
        """
        if [ -s {input.blast_out} ]
        then   
            sh -c 'silix "{input.fasta}" "{input.blast_out}" -f "{wildcards.seed}"\
               -i "{wildcards.pid}" -r "{wildcards.coverage}" -q "{params.minimum_overlap}"\
               -s "{params.minimum_percId}" -l "{params.minimum_length}" > "{output}" 2> {log}'
        else
            touch '{output}'
        fi
        """


##########################################################################
##########################################################################
