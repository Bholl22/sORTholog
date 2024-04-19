import pandas as pd
import numpy as np
import sys, os

##########################################################################


def transform_proteinid(df):
    """
    get back the protein name without genomeid
    """

    for index, row in df.iterrows():
        protein_id = row.protein_id.split("--")[0]
        df.at[index, "protein_id"] = protein_id

    return df


##########################################################################


def proteins2csv(sub_df):
    """
    Test the object to not be null and concatenate the protein_id to one
    """

    if not pd.isna(sub_df).all():
        return ",".join(sub_df.tolist())
    else:
        return ""


##########################################################################


def find_neighbors(df_patab, protein_tab, nb_check):
    '''
    investigate if two matches in a single genomes are encoded close to each other
    :param patab: dataframe of presence abscence table
    :param protein_tab: protein list of all genomes, ordered by sequences in the genomes
    :param nb_check: number of neighbors to investigates
    :return: dataframe of presence abscence table updated for neighbors
    '''

    for assembly_id in df_patab.assembly_id.unique():
        patab_filtered = df_patab[df_patab['assembly_id'] == assembly_id]
        protein_tab_filtered = protein_tab[protein_tab['genome_id'] == assembly_id].reset_index()
        protein_tab_filtered['protein_id'] = protein_tab_filtered.protein_id.apply(lambda x: x.split('--')[0])
        all_prot = [*set(patab_filtered.protein_id.to_list())]
        if np.nan in all_prot:
            all_prot.remove(np.nan)
        elif "" in all_prot:
            all_prot.remove("")

        for protein in all_prot:
            index_patab = patab_filtered[patab_filtered["protein_id"] == protein].index[0]
            index_prot_tab = protein_tab_filtered[protein_tab_filtered["protein_id"] == protein].index[0]
            tot_number = protein_tab_filtered.shape[0] - 1

            if index_prot_tab <= nb_check:
                index2gather = [*range(0, index_prot_tab + nb_check + 1, 1)]
                index2gather.extend(range(tot_number - 10 + index_prot_tab, tot_number + 1, 1))

            elif index_prot_tab >= tot_number - nb_check:
                index2gather = [*range(index_prot_tab - nb_check, tot_number + 1, 1)]
                index2gather.extend(range(0, tot_number - index_prot_tab + 1, 1))

            else:
                index2gather = [*range(index_prot_tab - nb_check, index_prot_tab + nb_check + 1, 1)]

            index2gather.remove(index_prot_tab)
            neighbor_list = []

            for index in index2gather:
                neighbor_list.append(protein_tab_filtered.at[index, 'protein_id'])

            for other_prot in all_prot:

                if other_prot == protein:
                    next

                elif other_prot in neighbor_list:
                    df_patab.at[index_patab, 'neighbors'] = True

    return df_patab


##########################################################################

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Seed preparing
seed_table = pd.read_table(snakemake.input.seed_file)
seed_list = seed_table.seed.to_list()

seed_color_dict = seed_table.set_index("seed").color.to_dict()

# list of all proteins
all_proteins = pd.read_table(snakemake.input.protein_table, dtype="string")
all_proteins.fillna("", inplace=True)

# fnodes opening
fam_id_table = pd.DataFrame()

for diamond_file in snakemake.input.diamonds:
    tmp_df = pd.read_table(diamond_file, sep='\t')
    fam_id_table = pd.concat([fam_id_table, tmp_df])

# add the protein information from the protein table and genome table
fam_id_table = fam_id_table.merge(all_proteins, on="protein_id")

# Retrieved protein_id
fam_id_table = transform_proteinid(fam_id_table)

# Table with number
patab = pd.crosstab(index=fam_id_table["genome_id"], columns=fam_id_table["seed"])

# To add missing seed if not find
seed_missing = [seed for seed in seed_list if seed not in patab.columns]
patab.loc[:, seed_missing] = 0

patab = patab[seed_list].reset_index()
# Add the genome name to the table in case needed
patab = patab.merge(
    fam_id_table[["genome_id", "genome_name"]].drop_duplicates(), on="genome_id"
)

# Tmp tab for order genomes
patab_tmp = patab.copy()

# Find max value to change
print(patab_tmp[seed_list])
max_value = patab_tmp[seed_list].values.max()

# Change value > 2 by 1
patab_tmp = patab_tmp.replace(to_replace=range(2, max_value + 1), value=1)

# Create a temporary columns to know how many seed they have
patab_tmp["total"] = patab_tmp[seed_list].sum(axis=1)

# Sort the columns by seed_list and genomes name
ascending_bool = [False] + [False] * len(seed_list) + [True]
patab_tmp = patab_tmp.sort_values(
    by=["total"] + seed_list + ["genome_name"], ascending=ascending_bool
)

patab = patab.loc[patab_tmp.index, :]

patab = patab.melt(
    id_vars=["genome_id", "genome_name"], var_name="seed", value_name="PA"
)

# Order the table by genome id to be more readable after
patab = patab.set_index("genome_id").loc[patab.genome_id.unique(), :].reset_index()

# Put color instead of number
for index, row in patab.iterrows():
    # Use the fact that 0 == False in python to test if it's 1 or 0
    if row.PA:
        patab.at[index, "color"] = seed_color_dict[row.seed]
    else:
        patab.at[index, "color"] = "#FFFFFF"  # White color

patab = patab.merge(
    fam_id_table[["genome_id", "seed", "protein_id"]],
    on=["genome_id", "seed"],
    how="left",
)

# Change genome_id to assembly_id
patab = patab.rename(columns={"genome_id": "assembly_id"})

# neighbor analysis
patab['neighbors'] = False
if snakemake.params.neighbor_analysis:
    patab = find_neighbors(df_patab=patab,
                           protein_tab=all_proteins,
                           nb_check=snakemake.params.neighbor_number,
                           )

# save the table
patab.to_csv(snakemake.output.final_table, sep="\t", index=False)

# pivot format with the name of the protein in cells
patab_table = patab.pivot_table(
    index="assembly_id",
    columns="seed",
    values="protein_id",
    aggfunc=proteins2csv,
    sort=False,
    dropna=False,
).reset_index()

fam_id_table = fam_id_table.rename(columns={"genome_id": "assembly_id"})

patab_table = patab_table.merge(
    fam_id_table[["assembly_id", "genome_name"]].drop_duplicates(), on="assembly_id"
)

patab_table[["assembly_id", "genome_name"] + seed_list].to_csv(
    snakemake.output.final_table_2, sep="\t", index=False
)
