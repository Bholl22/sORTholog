# Rule to make fasta from table

##########################################################################
##########################################################################

rule diamond_database:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER, "databases", "merge_fasta", "all_protein_with_seeds.fasta"
        ),
    output:
        db=os.path.join(
            OUTPUT_FOLDER, "databases", "diamond", "all_protein_with_seeds.dmnd"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "diamond",
            "diamond_database.log",
        ),
    conda:
        "../envs/diamond.yaml"
    shell:
        """diamond makedb --in {input.fasta} -d {output.db} 2> {log}"""


##########################################################################
##########################################################################

rule diamond_cluster:
    input:
        db=os.path.join(OUTPUT_FOLDER,"databases","diamond","all_protein_with_seeds.dmnd"),
    output:
        cluster=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processing_files",
                "diamond",
                "all_protein_with_seeds_{gene_constrains}.cluster",
            ),
            gene_constrains=gene_constrains,
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "diamond",
            "diamond_cluster.log",
        ),
    params:
        dir = os.path.join(OUTPUT_FOLDER,"processing_files","diamond"),
    resources:
        cpus=5,
        mem_mb=25000,
        time=7200,
    conda:
        "../envs/diamond.yaml"
    threads: 5
    shell:
        """
        mkdir -p {params.dir}
        for OUTPUT in {output.cluster}
        do
            EVALUE=$(basename $OUTPUT .cluster | awk -F 'evalue_' '{{print ($2)}}' | awk -F '_' '{{print ($1)}}')
            COVERAGE=$(basename $OUTPUT .cluster | awk -F 'cov_' '{{print ($2)}}' | awk -F '_' '{{print ($1)}}')
            PIDENT=$(basename $OUTPUT .cluster | awk -F 'pid_' '{{print ($2)}}' | awk -F '_' '{{print ($1)}}')
        
            diamond deepclust -d {input.db} -o $OUTPUT -M 4G -p {threads} --header --cluster-steps ultra-sensitive \
                --approx-id $PIDENT --member-cover $COVERAGE -e $EVALUE \
                --masking 0 --soft-masking 0 --no-block-size-limit --ext full \
                2>> {log}
                
            echo $OUTPUT
        done
         """


##########################################################################
##########################################################################

rule diamond_realign:
    input:
        cluster=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "diamond",
            "all_protein_with_seeds_{seed}_evalue_{eval}_cov_{cov}_pid_{pid}.cluster"
        ),
        db=os.path.join(
            OUTPUT_FOLDER,"databases","diamond","all_protein_with_seeds.dmnd"
        ),
    output:
        table=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "diamond",
            "all_protein_with_seeds_{seed}_evalue_{eval}_cov_{cov}_pid_{pid}.cluster.tsv"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "diamond",
            "diamond_realign_{seed}_evalue_{eval}_cov_{cov}_pid_{pid}.log",
        ),
    conda:
        "../envs/diamond.yaml"
    threads: 5
    shell:
        """diamond realign --clusters {input.cluster} -d {input.db} -o {output.table} \
         -f 6 -M 40G -p {threads} --header \
         2> {log}"""


##########################################################################
##########################################################################
