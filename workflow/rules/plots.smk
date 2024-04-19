# Rules for plotting the results

##########################################################################
##########################################################################


rule plots:
    input:
        final_table=os.path.join(OUTPUT_FOLDER, "results", "patab_melt.tsv"),
    output:
        png=report(
            os.path.join(OUTPUT_FOLDER, "results", "plots", "gene_PA.png"),
            caption="../report/PA_plot.rst",
            category="Plots",
        ),
        pdf=report(
            os.path.join(OUTPUT_FOLDER, "results", "plots", "gene_PA.pdf"),
            caption="../report/PA_plot.rst",
            category="Plots",
        ),
    resources:
        mem_mb=10000,
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "plots", "plots.log"),
    conda:
        "../envs/pandas_plots.yaml"
    script:
        "../scripts/plot-PA.py"


##########################################################################
##########################################################################


rule user_plots:
    input:
        final_table=config["PAtab_table"],
    output:
        png=os.path.join(OUTPUT_FOLDER, "results", "plots", "user_gene_PA.png"),
        pdf=os.path.join(OUTPUT_FOLDER, "results", "plots", "user_gene_PA.pdf"),
    params:
        color=config["default_values_plot"]["color"],
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "plots", "user_plots.log"),
    conda:
        "../envs/pandas_plots.yaml"
    script:
        "../scripts/quick_plotPA.py"


##########################################################################
##########################################################################
