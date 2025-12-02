#############
# FUNCTIONS #
#############
def compute_probe_matrix(genotypage, sequencage, list_sample: list) -> list:
    """
    Compute concordance matrix for a probe, and return it in a list of list format,
    as well as a counter of heterozygote sequencing.
    """
    # All matrix emplacement, line by line.
    # Genotypage 00, 01, 11, NA, Total
    matrix = [
        [0, 0, 0, 0, 0],  # Sequençage 00
        [0, 0, 0, 0, 0],  # 01
        [0, 0, 0, 0, 0],  # 11
        [0, 0, 0, 0, 0],  # NA
        [0, 0, 0, 0, 0],  # Total
    ]
    for sample in list_sample:
        sample_genotypage = dict(genotypage.samples[sample])["GT"]
        sample_sequencage = dict(sequencage.samples[sample])["GT"]
        # If sequencage == NA, second to last col
        if None not in sample_sequencage:
            seq_a = sample_sequencage[0]
            seq_b = sample_sequencage[1]
            coord_line = seq_a + seq_b
        else:
            coord_line = -2
        # If genotypage == NA, col -2
        if sample_genotypage[0] is not None:
            coord_col = sample_genotypage[0] + sample_genotypage[1]
        else:
            coord_col = -2
        # When correct place as been determined, increment it.
        matrix[coord_line][coord_col] += 1
    # Matrix completed, now we just have to compute totals
    matrix = compute_matrix_totals(matrix)

    return matrix


def compute_matrix_totals(matrix: list) -> list:
    """
    Take a square matrix and return total for line and column on the last line/column
    """
    lim_line = len(matrix) - 1
    lim_col = len(matrix[0]) - 1
    # Line/Col matrix
    for i in range(0, lim_line):
        matrix[i][lim_col] = sum(matrix[i])
    for i in range(0, lim_col):
        matrix[lim_line][i] = sum([line[i] for line in matrix])
    # Grand total
    matrix[lim_line][lim_col] = sum(matrix[lim_line][:-1])

    return matrix
