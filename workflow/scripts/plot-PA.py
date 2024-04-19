import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import FancyBboxPatch
import colorsys
import sys
import numpy as np

# It seems there is a bug if another backend is used
import matplotlib

matplotlib.use("Agg")

##########################################################################


def contrasting_text_color(hex_str):
    """
    Input a string without hash sign of RGB hex digits to compute
    complementary contrasting color such as for fonts
    """

    (r, g, b) = (hex_str[1:3], hex_str[3:5], hex_str[5:])

    luminosity = (
        1 - (int(r, 16) * 0.299 + int(g, 16) * 0.587 + int(b, 16) * 0.114) / 255
    )

    return "#131516" if luminosity < 0.5 else "white"


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


def invert_color(hex_color):
    """
    Invert the color represented by the given hexadecimal color code.

    :param hex_color: Hexadecimal color code (e.g. "#FFFFFF" for white).
    :return: Inverted color represented as a hexadecimal color code.
    """
    # Convert the hex color code to RGB
    r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)

    # Invert the color by subtracting the RGB values from 255
    inverted_r, inverted_g, inverted_b = 255 - r, 255 - g, 255 - b

    # Convert the inverted RGB values to hex
    inverted_hex_color = "#{:02x}{:02x}{:02x}".format(inverted_r, inverted_g, inverted_b)

    return inverted_hex_color


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

# Name in PAtab assembly_id, seed, PA, color, genome_name
patab_dtypes = {
    "assembly_id": "object",
    "seed": "object",
    "PA": np.int8,
    "color": "string",
    "genome_name": "string",
    "neighbors": "bool",
}

patab = pd.read_table(snakemake.input.final_table, dtype=patab_dtypes)
patab.fillna("", inplace=True)

# Dict position genomes and gene
list_genome = patab.assembly_id.unique().tolist()

### Amelioration purposes : The genome name could be followed by the id in parenthesis
### or we could ask the user for id or name in the figure

# As some genome name could be the same I need to create a dict to convert the name
genomeId_2_genomeName = patab.set_index("assembly_id").genome_name.to_dict()
list_genome_name = [
    f"{genomeId_2_genomeName[genome]} ({genome})" for genome in list_genome
]

num_genome = len(list_genome)
# here fist genome on top
dict_pos_genome = {
    list_genome[index]: num_genome - index - 1 for index in range(num_genome)
}

list_seed = patab.seed.unique().tolist()
num_seed = len(list_seed)
dict_pos_seed = {list_seed[index]: index for index in range(num_seed)}

# Try to have the magic figure size
leftmargin = 0.5  # inches
rightmargin = 0.3  # inches
topmargin = 0.5
bottommargin = 0.1  # inches
categorysize = 0.25  # inches

figwidth = leftmargin + rightmargin + (num_seed + 1) * categorysize
figheight = topmargin + bottommargin + (num_genome + 1) * categorysize

# Allows to controle the size of the box marker
size_rec = 0.8

# calculating dpi max pixel == 2**16
dpi = 300 if 300 * figheight < 2**16 else 2**16 // (figheight + 1)

# figsize = (width, height) plosBio ((7.5, 8.75))
fig, ax = plt.subplots(1, 1, figsize=(figwidth, figheight), dpi=dpi)
fig.subplots_adjust(
    left=leftmargin / figwidth,
    right=1 - rightmargin / figwidth,
    top=1 - topmargin / figheight,
    bottom=bottommargin / figheight,
)

label_format = {"fontweight": "bold"}
neighbor_line = snakemake.config["neighborhoods"]["line"]

for _, row in patab.iterrows():
    edge_color = "#131516"
    line_shape = '-'

    if row.neighbors:
        line_shape = neighbor_line

    if snakemake.config["default_values_plot"]["colored_border"]:
        # Change the border's shade to a darker color infer from the background color
        edge_color = (
            "#2F3D44" if row.color == "#FFFFFF" else adjust_lightness(row.color)
        )


    # Change the border's shape to a round version
    if snakemake.config["default_values_plot"]["round_border"]:
        boxstyle = "round,pad=-0.0040,rounding_size=2"
    else:
        boxstyle = "round,pad=-0.04"

    ax.add_artist(
        FancyBboxPatch(
            xy=(
                dict_pos_seed[row.seed] - size_rec / 2,
                dict_pos_genome[row.assembly_id] - size_rec / 2,
            ),
            facecolor=row.color,
            boxstyle=boxstyle,
            mutation_scale=0.2,
            width=size_rec,
            height=size_rec,
            edgecolor=edge_color,
            linestyle=line_shape,
            lw=1,
        )
    )

    if row.PA > 1:
        ax.text(
            x=dict_pos_seed[row.seed],
            y=dict_pos_genome[row.assembly_id],
            s=str(row.PA),
            color=contrasting_text_color(row.color),
            ha="center",
            va="center",
            fontweight="heavy",
        )

plt.yticks(range(num_genome), list_genome_name[::-1], **label_format)
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
    plt.savefig(plot_name, bbox_inches="tight")