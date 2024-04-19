import pandas as pd
import sys
import numpy as np
from numba import jit
import numpy.lib.recfunctions as rfn

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################
# Try to load it from file but not good now because can't access to the file
# if deployed worlflow by snakedeploy
##########################################################################


def iterrator_on_blast_hsp(blast_out):
    """Iterate over blast output. It considers that the output
    is in outfmt 6 and that all the hsp should be one after the
    other

    Args:
        blast_out (string): Path to a blast tabulated output

    Yields:
        pandas.DataFrame: A pandas Dataframe of the two proteins and their hsps
    """

    current_pair = ""
    list_line = []

    with open(blast_out) as r_file:

        for line in r_file:
            split_line = line.rstrip().split("\t")

            if (split_line[0], split_line[1]) == current_pair:
                list_line.append(split_line)
            elif current_pair == "":
                current_pair = (split_line[0], split_line[1])
                list_line = [split_line]
            else:
                yield list_line
                current_pair = (split_line[0], split_line[1])
                list_line = [split_line]

        yield list_line


##########################################################################


@jit(nopython=True)
def get_orientation(sstart, send):
    """Returns the orientation between two points.

    Args:
        sstart (numpy.array): Start of the sequence
        send (numpy.array): End of the sequence

    Returns:
        numpy.array: 1 or -1 depending on the orientation
    """

    return (send - sstart) // np.abs((send - sstart))


##########################################################################


@jit(nopython=True)
def update_values(old_values, sens):
    """Update sequence start/end with depending of the orientation.

    Args:
        old_values (numpy.array): Position in the sequence
        sens (numpy.array): 1 or -1 depending on the orientation

    Returns:
        np.array: Position in the sequence
    """

    return old_values * sens


##########################################################################


@jit(nopython=True)
def length_aln_on_sequence(start, end):
    """Return the length of the sequence.

    Args:
        start (numpy.array): Start of the sequence
        end (numpy.array): End of the sequence

    Returns:
        np.array: Length of the protein
    """
    return end - start + 1


##########################################################################


@jit(nopython=True)
def calculateIdentities(percIdentity, length):
    """Return the length of the sequence.

    Args:
        percIdentity (numpy.array): Percentage of identity from blast
        length (numpy.array): length of the alignment

    Returns:
        np.array: number of identical amino acids in the alignment
    """
    return np.floor(percIdentity * length / 100 + 0.5)


##########################################################################


@jit(nopython=True)
def calculate_fraction(delta, lgHSP, bitscore, pid, pos):
    """Calculate the new score, identities and positives of the hsp2.

    Args:
        delta (int): Difference of size between old and new hsp2 length
        lgHSP (int): size of the hsp before trimming
        bitscore (float): Score of the alignment
        pid (int): Number of identical in the alignment
        pos (int): Number of positives in the alignment

    Returns:
        float, int, int, int: The updated values of the score, identities, positives and length
    """

    # Calculate new score, id and positive
    # Calculation: initial_value * (franction of length that has been preserved)
    fraction = 1 - (delta / lgHSP)

    new_score = np.floor(bitscore * fraction)
    new_id = np.floor(pid * fraction)
    new_pos = np.floor(pos * fraction)

    # Calculate new length
    new_length = lgHSP - delta

    # Set expect value to -1 : this value should not be used after
    # having changed HSPs boundaries
    # new_evalue = -1

    return new_score, new_id, new_pos, new_length


##########################################################################


def remove_overlap_query(hsp1, hsp2):
    """Remove overlap with the given hsp2 in the query.

    Args:
        hsp1 (np.array): Blast values of the given hsp (best)
        hsp2 (np.array): Blast values of the given hsp (questioning)

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    # Calculate where is the overlap and remove the overlapping part: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    if hsp2[6] < hsp1[6]:
        new_qstart = hsp2[6]
        new_qend = hsp1[6] - 1
        delta = hsp2[7] - new_qend
        new_sstart = hsp2[8]
        new_send = hsp2[9] - delta

    elif hsp2[7] > hsp1[7]:
        new_qstart = hsp1[7] + 1
        new_qend = hsp2[7]
        delta = new_qstart - hsp2[6]
        new_sstart = hsp2[8] + delta
        new_send = hsp2[9]

    # lgHSP: 17, bitscore: 11, id: 12, pos:13
    new_score, new_id, new_pos, new_length = calculate_fraction(
        delta=delta, lgHSP=hsp2[17], bitscore=hsp2[11], pid=hsp2[12], pos=hsp2[13]
    )

    return {
        11: new_score,
        17: new_length,
        7: new_qend,
        6: new_qstart,
        9: new_send,
        8: new_sstart,
        13: new_pos,
        12: new_id,
    }


##########################################################################


def remove_overlap_subject(hsp1, hsp2):
    """Remove overlap with the given hsp2 in the subject.

    Args:
        hsp1 (pandas.Series): Blast values of the given hsp (best)
        hsp2 (pandas.Series): Blast values of the given hsp (questioning)

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    # Calculate where is the overlap and remove the overlapping part: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    if hsp2[8] < hsp1[8]:
        new_sstart = hsp2[8]
        new_send = hsp1[8] - 1
        delta = hsp2[9] - new_send
        new_qstart = hsp2[6]
        new_qend = hsp2[7] - delta

    elif hsp2[9] > hsp1[9]:
        new_sstart = hsp1[9] + 1
        new_send = hsp2[9]
        delta = new_sstart - hsp2[8]
        new_qstart = hsp2[6] + delta
        new_qend = hsp2[7]

    # lgHSP: 17, bitscore: 11, id: 12, pos:13
    new_score, new_id, new_pos, new_length = calculate_fraction(
        delta=delta, lgHSP=hsp2[17], bitscore=hsp2[11], pid=hsp2[12], pos=hsp2[13]
    )

    return {
        11: new_score,
        17: new_length,
        7: new_qend,
        6: new_qstart,
        9: new_send,
        8: new_sstart,
        13: new_pos,
        12: new_id,
    }


##########################################################################


def checkHSPS(hsp1, hsp2, HSPMIN=100):
    """compare two HSPS in the blast output

    Args:
        hsp1 (np.array): Blast values of the given hsp (best)
        hsp2 (np.array): Blast values of the given hsp (questioning)
        HSPMIN (int, optional): Minumum length of the HSP. Defaults to 100.

    Returns:
        dict: Dictionnary with the updated values for the hsp2
    """
    dict_update = {18: 0}

    # Check if the hsp2 is in a different orientation than hsp1: 'sens': 14
    if hsp1[14] != hsp2[14]:
        # print(f'orientation wrong: {hsp1[14]} != {hsp2[14]}')
        return dict_update

    # Check is hsp2 inside hsp1 for the query sequence: 'qstart': 6, 'qend': 7,
    if hsp1[6] >= hsp2[6] and hsp1[7] <= hsp2[7]:
        # print(f'hsp2 inside hsp1 for query: {hsp1[6]} >= {hsp2[6]} and {hsp1[7]} <= {hsp2[7]}')
        return dict_update

    # Check is hsp1 inside hsp2 for the query sequence: 'qstart': 6, 'qend': 7,
    elif hsp1[6] <= hsp2[6] and hsp1[7] >= hsp2[7]:
        # print(f'hsp1 inside hsp2 for query: {hsp1[6]} <= {hsp2[6]} and {hsp1[7]} >= {hsp2[7]}')
        return dict_update

    # Check is hsp1 inside hsp2 for the subject sequence: 'sstart': 8, 'send': 9,
    elif hsp1[8] >= hsp2[8] and hsp1[9] <= hsp2[9]:
        # print(f'hsp2 inside hsp1 for subject: {hsp1[8]} >= {hsp2[8]} and {hsp1[9]} <= {hsp2[9]}')
        return dict_update

    # Check is hsp2 inside hsp1 for the subject sequence: 'sstart': 8, 'send': 9,
    elif hsp1[8] <= hsp2[8] and hsp1[9] >= hsp2[9]:
        # print(f'hsp1 inside hsp2 for subject: {hsp1[8]} <= {hsp2[8]} and {hsp1[9]} >= {hsp2[9]}')
        return dict_update

    # reject HSPs that are in different orientation: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    # Query:    ---- A ---- B -----   A = HSP1
    # Sbjct:    ---- B ---- A -----   B = HSP2

    if (hsp1[7] - hsp2[7]) * (hsp1[9] - hsp2[9]) < 0:
        # print(f'HSPs are in different orientation 1: ({hsp1[7]} - {hsp2[7]}) * ({hsp1[9]} - {hsp2[9]}) ===> {(hsp1[7] - hsp2[7]) * (hsp1[9] - hsp2[9])} < 0')
        return dict_update
    elif (hsp1[6] - hsp2[6]) * (hsp1[8] - hsp2[8]) < 0:
        # print(f'HSPs are in different orientation 2: ({hsp1[6]} - {hsp2[6]}) * ({hsp1[8]} - {hsp2[8]}) ===> {(hsp1[6] - hsp2[6]) * (hsp1[8] - hsp2[8])} < 0')
        return dict_update

    overlap_q = (hsp2[6] - hsp1[7]) * (hsp2[7] - hsp1[6])
    overlap_s = (hsp2[8] - hsp1[9]) * (hsp2[9] - hsp1[8])

    # Accept non-overlapping HSPs in correct orientation
    if overlap_q > 0 and overlap_s > 0:
        # print(f'No overlap in query and subject: {overlap_q} > 0 and {overlap_s} > 0')
        dict_update[18] = 1
        return dict_update

    # Test if the query is overlaping
    elif overlap_q < 0:
        # print(f'Overlap in query: {overlap_q} > 0')
        dict_update = remove_overlap_query(hsp1=hsp1, hsp2=hsp2)

        # update the hsp2 array with the new values
        for index_key in dict_update:
            hsp2[index_key] = dict_update[index_key]

        overlap_s = (hsp2[8] - hsp1[9]) * (hsp2[9] - hsp1[8])

    # Test if the subject is overlaping after possible update of an overlaping query
    if overlap_s < 0:
        # print(f'Overlap in subject: {overlap_s} > 0')
        dict_update = remove_overlap_subject(hsp1=hsp1, hsp2=hsp2)

        # update the hsp2 array with the new values
        for index_key in dict_update:
            hsp2[index_key] = dict_update[index_key]

    # Filter out HSPs that are too short
    if hsp2[17] < HSPMIN:
        # print(f'HSP too short: {hsp2[17]} < {HSPMIN}')
        dict_update[18] = 0
        return dict_update

    # Set status to 1 for consistent HSPs
    dict_update[18] = 1
    return dict_update


##########################################################################


def prepare_df_hsps(list_hsps, blast_dtypes, HSPMIN=100):
    """Prepare a dataframe containing HSPs and update values.

    Args:
        list_hsps (list of strings): List of the lines containing the HSPs for the same pair of proteins
        blast_dtypes (dict): Type of the different columns in the blast outfmt 6 output
        HSPMIN (int, optional): [description]. Defaults to 100.

    Returns:
        pandas.DataFrame: Reduce DataFrame with only the significant hsps
    """

    dtypes = np.dtype(blast_dtypes)

    df_hsps = np.array(list_hsps)
    df_hsps = rfn.unstructured_to_structured(df_hsps, dtype=dtypes)

    # Prepare the dtypes of the columns to add
    dtype2add = [
        ("id", np.int32),
        ("pos", np.int32),
        ("sens", np.int8),
        ("lg_aln_subject", np.int32),
        ("lg_aln_query", np.int32),
        ("lgHSP", np.int32),
        ("stat", np.int8),
    ]

    dtype2add = np.dtype(dtype2add)

    array2add = np.empty(df_hsps.shape, dtype=dtype2add)

    # Order HSPS by decreasing score
    # ind = np.argsort(df_hsps[:,11].astype(np.float64))[::-1]
    # ind = np.transpose([ind] * 18)
    # df_hsps = np.take_along_axis(df_hsps, ind, axis=0)
    ind = np.argsort(df_hsps, order="bitscore")[::-1]
    df_hsps = df_hsps[ind]

    # Calculate a more precise percentage of identity for each independent HSP: 'pident': 2, 'length': 3
    array2add["id"] = calculateIdentities(
        percIdentity=df_hsps["pident"], length=df_hsps["length"]
    )
    array2add["pos"] = array2add["id"]

    # Get the sens of the HSP on the subject: 'sstart': 8, 'send': 9
    array2add["sens"] = get_orientation(sstart=df_hsps["sstart"], send=df_hsps["send"])

    # Update the values of the position in case the sens is changed for the subject (as done in silix)
    df_hsps["sstart"] = update_values(
        old_values=df_hsps["sstart"], sens=array2add["sens"]
    )
    df_hsps["send"] = update_values(old_values=df_hsps["send"], sens=array2add["sens"])

    # Calculate the length of the alignment on the subject and query: 'qstart': 6, 'qend': 7, 'sstart': 8, 'send': 9,
    array2add["lg_aln_subject"] = length_aln_on_sequence(
        start=df_hsps["sstart"], end=df_hsps["send"]
    )
    array2add["lg_aln_query"] = length_aln_on_sequence(
        start=df_hsps["qstart"], end=df_hsps["qend"]
    )

    # Calculate length of the HSP: lg_aln_subject: 15, lg_aln_query:16
    array2add["lgHSP"] = np.max(
        [array2add["lg_aln_subject"], array2add["lg_aln_query"]], axis=0
    )

    # Put Stats value to know which HSPS to conserve
    array2add["stat"] = np.array([1] + [-1] * (df_hsps.shape[0] - 1))

    # Concat df_hsps and array2add
    df_hsps = rfn.merge_arrays([df_hsps, array2add], flatten=True, usemask=False)

    for hsp1_index in range(df_hsps.shape[0] - 1):
        for hsp2_index in range(hsp1_index + 1, df_hsps.shape[0]):
            if df_hsps["stat"][hsp1_index] and df_hsps["stat"][hsp2_index]:
                values2update = checkHSPS(
                    hsp1=df_hsps[hsp1_index], hsp2=df_hsps[hsp2_index], HSPMIN=HSPMIN
                )

                # print(hsp1_index, hsp2_index, values2update)
                for value in values2update:
                    df_hsps[hsp2_index][value] = values2update[value]

    return df_hsps, ind


##########################################################################
##########################################################################

# Seeds
seed_table = pd.read_table(snakemake.input.seed_file).set_index("protein_id")

# Protein lenght dict
# protein_dict = {}

# with open(snakemake.input.protein_file, 'rt') as r_file :
#     r_file.readline()
#     for line in r_file:
#         r_line = line.rstrip()
#         split_line = r_line.split('\t')

#         protein_dict[split_line[0]] = int(split_line[-1])

# protein_dict.update(seed_table.length
#                               .to_dict())

# Blast_out and preparation
blast_names = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

# Get the types of te columns for multiple HSPs dataframe
blast_dtypes = [
    ("qseqid", "S100"),
    ("sseqid", "S100"),
    ("pident", np.float64),
    ("length", np.int32),
    ("mismatch", np.int32),
    ("gapopen", np.int32),
    ("qstart", np.int32),
    ("qend", np.int32),
    ("sstart", np.int32),
    ("send", np.int32),
    ("evalue", np.float64),
    ("bitscore", np.float64),
]

# Calculating the max and min value of interest
max_eval = seed_table.evalue.max()
# min_pident = seed_table.pident.min()
# min_coverage = seed_table.coverage.min()

# Open all the output
list_open_output = [open(tsv_output, "wt") for tsv_output in snakemake.output]
num_output = len(list_open_output)

# Read the blast hsp by hsp
for sub_blast in iterrator_on_blast_hsp(blast_out=snakemake.input.blast_out):
    # Variable to get the lines
    line = ""

    # Get the number of hsps
    num_HSPs = len(sub_blast)

    if num_HSPs == 1:
        evalue_blast = float(sub_blast[0][10])
        line = "\t".join(sub_blast[0]) + "\n"
    else:
        df_hsps, reorder = prepare_df_hsps(
            list_hsps=sub_blast,
            blast_dtypes=blast_dtypes,
            HSPMIN=snakemake.params.minimum_length,
        )

        evalue_blast = df_hsps["evalue"][0]

        # Test which index I kept from initial
        true_false = df_hsps["stat"] == 1

        for index in range(df_hsps.shape[0]):
            if true_false[index]:
                # here index of the return df_hsps are potentially reorder so get the old index
                line += "\t".join(sub_blast[reorder[index]]) + "\n"

    if evalue_blast <= max_eval:
        # if evalue_blast <= max_eval and pident_blast >= min_pident :
        # coverage_blast = length / protein_dict[qseqid]

        # Calculating the coverage on the query
        # if coverage_blast >= min_coverage :
        for index in range(num_output):
            file_out_path = snakemake.output[index]

            constrains = file_out_path.split("--")[-1].split(".out")[0]

            constrains_split = constrains.split("_")

            # As i have the control about the end of the file name it is better to go this way in case "_" in the gene name
            evalue = float(constrains_split[-5])
            # coverage = float(constrains_split[-3])
            # pident = float(constrains_split[-1])

            # if evalue_blast <= evalue and pident_blast >= pident and coverage_blast >= coverage:
            if evalue_blast <= evalue:
                list_open_output[index].write(line)

for file_open in list_open_output:
    file_open.close()
