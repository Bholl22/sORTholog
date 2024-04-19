##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
import pandas as pd
import numpy as np
from snakemake.utils import validate
from snakemake.workflow import workflow

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################


def get_final_output():
    """
    Generate final output name
    """
    final_output = multiext(
        os.path.join(OUTPUT_FOLDER, "results", "plots", "gene_PA"), ".png", ".pdf"
    )

    if config['threshold_analysis']:
        final_output.append(
            os.path.join(OUTPUT_FOLDER, "analysis_thresholds", "report_figure_thresholds.html")
        )

    return final_output


##########################################################################


def infer_gene_constrains(seed_df):
    """
    Infer gene_constrains from default config value or table
    """

    list_constrains = []

    for index, row in seed_df.iterrows():
        if "evalue" in seed_df.columns and not pd.isna(row.evalue):
            tmp_evalue = row.evalue
        else:
            tmp_evalue = config["default_blast_options"]["e_val"]
            seed_df.at[index, "evalue"] = tmp_evalue

        if "coverage" in seed_df.columns and not pd.isna(row.coverage):
            tmp_coverage = row.coverage
        else:
            tmp_coverage = config["default_blast_options"]["cov"]
            seed_df.at[index, "coverage"] = tmp_coverage

        if "pident" in seed_df.columns and not pd.isna(row.pident):
            tmp_pident = row.pident
        else:
            tmp_pident = config["default_blast_options"]["pid"]
            seed_df.at[index, "pident"] = tmp_pident

        tmp_text = (
            f"{row.seed}_evalue_{tmp_evalue:.0e}_cov_{tmp_coverage}_pid_{tmp_pident}"
        )

        list_constrains.append(tmp_text)

    return list_constrains, seed_df


##########################################################################


def check_color_seed(seed_df):
    """
    Infer color if color is not set by the user in the seed's file
    """

    if "color" not in seed_df.columns:
        seed_df["color"] = config["default_values_plot"]["color"]
    else:
        seed_df.fillna(
            value={"color": config["default_values_plot"]["color"]}, inplace=True
        )

    return seed_df


##########################################################################


def compare_seed_table(seed_df, new_seed_file, start_seed_file, seed_dtypes):
    """
    Compare the seed and new seed if exists to update the new_seed
    Restart the pipeline from start if:
        - New seed file not found
        - Seed file and new seed file don't have the same number of seeds
        - Protein id does not match
    Else:
        - Update new seed file
    """

    columns2change = ["seed", "evalue", "pident", "coverage", "color"]

    if os.path.isfile(new_seed_file):
        new_seed_df = pd.read_table(new_seed_file, dtype=seed_dtypes)
        start_seed_df = pd.read_table(start_seed_file, dtype=seed_dtypes)

        # Because bug it might happen that the empty slot are NA
        start_seed_df.fillna("", inplace=True)
        seed_df.fillna("", inplace=True)
        new_seed_df.fillna("", inplace=True)

        seed_df = seed_df.astype(seed_dtypes)

        msg_rerun = "If you want to rerun sORTolog with a new seed file, please change the project name or remove the old project folder"

        # If seed is added
        if seed_df.shape[0] != start_seed_df.shape[0]:
            # seed_df.to_csv(start_seed_file, sep="\t", index=False)
            sys.exit(
                f"WARNING:: Your seed file is different from your last run\nWhy? Because a seed was removed or added\n{msg_rerun}"
            )
        # If protein name change
        elif not seed_df.protein_id.equals(start_seed_df.protein_id):
            # seed_df.to_csv(start_seed_file, sep="\t", index=False)
            sys.exit(
                f"WARNING:: Your seed file is different from your last run\nWhy? Because not the same protein ids\n{msg_rerun}"
            )
        # If hmm name change
        elif not seed_df.hmm.equals(start_seed_df.hmm):
            # seed_df.to_csv(start_seed_file, sep="\t", index=False)
            sys.exit(
                f"WARNING:: Your seed file is different from your last run\nWhy? Because HMM column changed (added or deleted hmm)\n{msg_rerun}"
            )
        # If something else change
        elif not seed_df[columns2change].equals(new_seed_df[columns2change]):
            # Update new seed with information of seed
            new_seed_df.update(seed_df[columns2change])
            new_seed_df.to_csv(new_seed_file, sep="\t", index=False)
    else:
        seed_df.to_csv(start_seed_file, sep="\t", index=False)

    return


##########################################################################


def get_list_hmm(seed_df):
    """
    Gather the list of HMM files from a folder and make sure they are in the seed table
    Update the seed table with the proper hmm profile file name
    """

    if "hmm" not in seed_df.columns:
        seed_df["hmm"] = ""
    else:
        seed_df.fillna(value={"hmm": ""}, inplace=True)

    # Check the index where hmm are
    index_hmm = seed_df.hmm != ""

    list_psiblast = seed_df[~index_hmm].protein_id.tolist()
    list_hmm = seed_df[index_hmm].hmm.tolist()

    return list_hmm, list_psiblast, seed_df


##########################################################################


def create_folder(mypath):
    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


##########################################################################


def check_annotation(mypath):
    """
    Check the presence of annotation file or if at least the name are in
    right format in the fasta file
    :param mypath: path to perso_database if exists
    :type: string
    :return: The path of the annotation file if exists and validated or empty string
             if fasta formated
    """

    with open(mypath) as r_file:
        first_header = r_file.readline()

        if "perso_annotation" in config and os.path.isfile(config["perso_annotation"]):
            perso_annotation = pd.read_table(config["perso_annotation"], dtype="string")
            validate(perso_annotation, schema="../schemas/annotations.schema.yaml")

            return config["perso_annotation"]
        elif "--" in first_header:
            return ""
        else:
            sys.exit(
                "ERROR: Please provided an annotation file for your database or \
                     format the header of the fasta file as: sequence_name--genome_id \
                     description"
            )


##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

# path to seeds sheet (TSV format, columns: seed, protein_id, ...)
seed_file = config["seed"]

# Validation of the seed file
seed_dtypes = {
    "seed": "string",
    "protein_id": "string",
    "hmm": "string",
    "evalue": np.float64,
    "pident": np.float64,
    "coverage": np.float64,
    "color": "string",
}

seed_table = pd.read_table(seed_file, dtype=seed_dtypes)

# Definition of the requirements for each seed
gene_constrains, seed_table = infer_gene_constrains(seed_table)

# Check color of the seeds
seed_table = check_color_seed(seed_table)

HMM, PSIBLAST, seed_table = get_list_hmm(seed_table)

validate(seed_table, schema="../schemas/seeds.schema.yaml")

##########################################################################

if "taxid" in config and os.path.isfile(config["taxid"]):
    # path to taxonomic id to search seeds in (TSV format, columns: TaxId, NCBIGroups)
    taxid = config["taxid"]

    # Validation of the taxid file
    taxid_dtypes = {
        "TaxId": "Int64",
        "NCBIGroups": "string",
    }

    taxid_table = pd.read_table(taxid, dtype=taxid_dtypes)

    validate(taxid_table, schema="../schemas/taxid.schema.yaml")
else:
    # Create empty DataFrame if file doesn't exists
    taxid_table = pd.DataFrame()

##########################################################################
##########################################################################
##
##                        Core configuration
##
##########################################################################
##########################################################################

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_basedir_short__"] = os.path.basename(workflow.basedir)
config["__workflow_workdir__"] = os.getcwd()
#Original
if workflow.config_settings.config_args:
    tmp_config_arg = '" '.join(workflow.config_settings.config_args).replace("=", '="')
    config["__config_args__"] = f' -C {tmp_config_arg}"'
else:
    config["__config_args__"] = ""

#Modified version

#ChatGPT Version
# if workflow.config_args:
#     tmp_config_arg = '" "'.join(workflow.config_args).replace("=", '="')
#     config["__config_args__"] = f' -C "{tmp_config_arg}"'
# else:
#     config["__config_args__"] = ""


with open(os.path.join(workflow.basedir, "../config/VERSION"), "rt") as version:
    url = "https://github.com/vdclab/sORTholog/releases/tag"
    config["__workflow_version__"] = version.readline()
    config["__workflow_version_link__"] = f"{url}/{config['__workflow_version__']}"


##########################################################################
##########################################################################
##
##                           Options
##
##########################################################################
##########################################################################

# Name your project
project_name = config["project_name"]

# Result folder
OUTPUT_FOLDER = os.path.join(config["output_folder"], project_name)
# Adding to config for report
config["__output_folder__"] = os.path.abspath(OUTPUT_FOLDER)

# Psiblast default e-value thershold
e_val_psiblast = config["default_psiblast_options"]["psiblast_e_val"]

# Psiblast default iteration thershold
iteration_psiblast = config["default_psiblast_options"]["iteration"]

# Silix option coverage
cov_min = config["silix_options"]["cov_min"]

# Silix option percentage identity
pid_min = config["silix_options"]["pid_min"]

# Silix option minimum length
length_min = config["silix_options"]["length_min"]

# Dictionary translation option for silix
silix_dict = {
    "mean": "0",
    "subject": "1",
    "query": "-1",
    "shortest": "2",
    "longest": "-2",
    "HSP": "3",
}

# Option for ncbi_genome_download
section = config["ndg_options"]["section"]

# Values for assembly_levels :
assembly_levels = config["ndg_options"]["assembly_levels"]

# Values for refseq_categories :
refseq_categories = config["ndg_options"]["refseq_categories"]

# Name of the file with all the taxids
starting_database = os.path.join(
    OUTPUT_FOLDER, "databases", "all_taxid", "taxid_all_together.fasta"
)

merge_db = os.path.join(
    OUTPUT_FOLDER,
    "databases",
    "merge_databases",
    "databases_all_together.fasta",
)

# Get the output protein table name to activate the needed rules
protein_table_taxid = os.path.join(
    OUTPUT_FOLDER, "databases", "all_taxid", "protein_table.tsv"
)

protein_table_merge = os.path.join(
    OUTPUT_FOLDER, "databases", "merge_databases", "protein_table.merged.tsv"
)


proteinTable = protein_table_merge

# Check if there is a database specified in the config file
if (
    "perso_database" in config
    and os.path.isfile(config["perso_database"])
    and "taxid" in config
    and not taxid_table.empty
):
    list_starting_database = [config["perso_database"], starting_database]
    annotationTable = [check_annotation(config["perso_database"]), protein_table_taxid]
elif "taxid" in config and not taxid_table.empty:
    list_starting_database = starting_database
    merge_db = starting_database
    proteinTable = protein_table_taxid
    annotationTable = protein_table_taxid
elif "perso_database" in config and os.path.isfile(config["perso_database"]):
    list_starting_database = [config["perso_database"]]
    annotationTable = [check_annotation(config["perso_database"])]
else:
    sys.exit("ERROR: Missing input file, no perso_database nor taxid table found")

# Compare seed_table and new_seed_table (if exists) to update e_val, cov, pident
new_seed_file = os.path.join(OUTPUT_FOLDER, "databases", "seeds", "new_seeds.tsv")

# Create a file as input of the first rule that change on if seeds.tsv change in required value
create_folder(os.path.join(OUTPUT_FOLDER, "databases", "seeds"))
start_seed_file = os.path.join(OUTPUT_FOLDER, "databases", "seeds", "start_seeds.tsv")

compare_seed_table(seed_table, new_seed_file, start_seed_file, seed_dtypes)

# HMM profile folder
hmm_folder = config["hmm_profiles"]

# Check HMMs exist
HMM = [os.path.join(hmm_folder, hmm) for hmm in HMM]

for hmm_file in HMM:
    if not os.path.isfile(hmm_file):
        sys.exit(f"ERROR:: The provided hmm file does not exists: {hmm_file}")

# HMM default e-value threshold
e_val_HMM = config["default_hmmsearch_options"]["e_val"]

# HMM type of filtering
hmm_type = "-E" if config["default_hmmsearch_options"]["focus"] == "full" else "--domE"

# Seepup option that create a reduce dataset using a psiblast step with the seed
if config["speedup"]:
    speedup = os.path.join(
        OUTPUT_FOLDER,
        "databases",
        "reduce_taxid",
        f"all_proteins_reduced.fasta",
    )

    # Once we know that there is a speedup just need to check if HMM or not
    if HMM and PSIBLAST:
        tsv_prot = os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "reduce_taxid",
            "list_all_proteins.tsv",
        )
    elif PSIBLAST:
        tsv_prot = os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "psiblast",
            f"list_all_proteins_psiblast--eval_{e_val_psiblast:.0e}.tsv",
        )
    else:
        tsv_prot = os.path.join(
            OUTPUT_FOLDER,
            "processing_files",
            "hmmsearch",
            f"list_all_proteins_hmmsearch--eval_{e_val_HMM:.0e}.tsv",
        )
else:
    speedup = merge_db
    tsv_prot = ""

if "default_threshold" in config:
    round_value=config["default_threshold"]["round_value"]
    min_lines=config["default_threshold"]["min_lines"]
else :
    round_value = False
    min_lines = False
