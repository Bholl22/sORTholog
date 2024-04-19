# Rule for plugin threshold

##########################################################################
##########################################################################


rule blast2threshold_table:
    input:
        seed_file=os.path.join(OUTPUT_FOLDER, "databases", "seeds", "new_seeds.tsv"),
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "diamond",
            "all_protein_with_seeds_{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.cluster.tsv"
        ),
        diamond=os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "diamond",
            "all_protein_with_seeds_{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.cluster.flushed",
        ),
        protein_file=proteinTable,
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "analysis_thresholds",
            "tables",
            "table_hits_family--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.tsv",
        ),
    params:
        option_cov=cov_min,
        option_pid=pid_min,
        minimum_length=length_min,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "analysis_thresholds",
            "blast2threshold_table--{seed}_evalue_{eval}_cov_{coverage}_pid_{pid}.log",
        ),
    resources:
        mem_mb=10000,
        time=720,
    conda:
        "../envs/pandas_plots.yaml"
    script:
        "../scripts/blast2threshold_table.py"


##########################################################################
##########################################################################


rule report_threshold:
    input:
        table_hits = expand(
            os.path.join(
                OUTPUT_FOLDER,
                "analysis_thresholds",
                "tables",
                "table_hits_family--{gene_constrains}.tsv",
            ),
            gene_constrains=gene_constrains,
        ),
        seeds = os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "seeds",
            "new_seeds.tsv"
        ),
    output:
        figure = os.path.join(
            OUTPUT_FOLDER,
            "analysis_thresholds",
            "report_figure_thresholds.html",
        ),
        seed_detection = os.path.join(
            OUTPUT_FOLDER,
            "analysis_thresholds",
            "threshold_detection_seeds.tsv",
        ) if config['threshold_detection']['active'] else [],

    params:
        css=workflow.source_path("../report/threshold_report.css"),
        round_value=round_value,
        min_lines=min_lines,
        threshold_detection=config['threshold_detection']['active'],
        var_on=config['threshold_detection']['dimentions'],
    log:
        os.path.join(
            OUTPUT_FOLDER, "logs", "analysis_thresholds", "analysis_thresholds_fig.log"
        ),
    conda:
        "../envs/plotly.yaml"
    script:
        "../scripts/blast2threshold_fig.py"


##########################################################################
##########################################################################
