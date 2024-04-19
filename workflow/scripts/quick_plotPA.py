import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import FancyBboxPatch
import colorsys
import sys

# It seems there is a bug if another backend is used
import matplotlib

matplotlib.use("Agg")

##########################################################################


def adjust_lightness(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


##########################################################################


# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# Plot parameters
plt.rcParams["text.color"] = "#131516"
plt.rcParams["svg.fonttype"] = "none"  # Editable SVG text
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.weight"] = "light"

# Name in PAtab genome_id, seed, PA, color, genome_name
patab = pd.read_table(
    snakemake.input.final_table, na_filter=False, index_col=0, sep="\t"
)
patab = patab.replace(".+", "1", regex=True).replace("", "0")

# Dict position genomes and gene
num_genome = patab.shape[0]
list_genome = patab.index.tolist()
dict_pos_genome = {
    list_genome[index]: num_genome - index - 1 for index in range(num_genome)
}

# here fist genome on top
num_seed = patab.shape[1]
list_seed = patab.columns.tolist()
dict_pos_seed = {list_seed[index]: index for index in range(num_seed)}


# Melt the table
patab = patab.reset_index()
old_index_col = patab.columns[0]
patab = patab.rename(columns={old_index_col: "genome_id"})
patab = patab.melt(id_vars="genome_id", var_name="seed", value_name="PA")
print(patab[patab.genome_id == "2 genomes"])
patab["color"] = patab.PA.map({"1": snakemake.params.color, "0": "white"})

# Try to have the magic figure size
leftmargin = 0.5  # inches
rightmargin = 0.3  # inches
topmargin = bottommargin = 0.1  # inches
categorysize = 0.25  # inches

figwidth = leftmargin + rightmargin + (num_seed + 1) * categorysize
figheight = topmargin + bottommargin + (num_genome + 1) * categorysize

# To update the size of the square
size_rec = 0.8

# figsize = (width, height) plosBio ((7.5, 8.75))
fig, ax = plt.subplots(1, 1, figsize=(figwidth, figheight))
fig.subplots_adjust(
    left=leftmargin / figwidth,
    right=1 - rightmargin / figwidth,
    top=1 - topmargin / figheight,
    bottom=bottommargin / figheight,
)

label_format = {"fontweight": "bold"}

for _, row in patab.iterrows():
    # Change the border's shade to a darker color infer from the background color
    if snakemake.config["default_values_plot"]["colored_border"]:
        edge_color = (
            "#2F3D44" if row.color == "#FFFFFF" else adjust_lightness(row.color)
        )
    else:
        edge_color = "#131516"

    # Change the border's shape to a round version
    if snakemake.config["default_values_plot"]["round_border"]:
        boxstyle = "round,pad=-0.0040,rounding_size=2"
    else:
        boxstyle = "round,pad=-0.04"

    # Change the border's shape to a round version
    if snakemake.config["default_values_plot"]["round_border"]:
        boxstyle = "round,pad=-0.0040,rounding_size=2"
    else:
        boxstyle = "round,pad=-0.04"

    ax.add_artist(
        FancyBboxPatch(
            xy=(
                dict_pos_seed[row.seed] - size_rec / 2,
                dict_pos_genome[row.genome_id] - size_rec / 2,
            ),
            facecolor=row.color,
            boxstyle=boxstyle,
            mutation_scale=0.2,
            width=size_rec,
            height=size_rec,
            edgecolor=edge_color,
            lw=1,
        )
    )

plt.yticks(range(num_genome), list_genome[::-1], **label_format)
plt.xticks(range(num_seed), list_seed, **label_format)

ax.tick_params(axis="both", which="both", length=0)  # No tick markers
ax.set_ylabel("")  # No ylabel
ax.xaxis.tick_top()  # xticklabels on top
ax.xaxis.set_label_position("top")
plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, ha="center")  # Rotate x labels

for pos in ["top", "bottom", "left", "right"]:
    ax.spines[pos].set_visible(False)  # Remove border

plt.xlim(-0.5, num_seed - 0.5)
plt.ylim(-0.5, num_genome - 0.5)

for plot_name in snakemake.output:
    plt.savefig(plot_name, bbox_inches="tight", dpi=300)
