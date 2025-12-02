def estimate_kappa_na(matrix: list) -> dict:
    """
    Take a concordance matrix and return a tuple containing p0, pc and kappa estimated
    on the whole matrix.
    """
    lim_line = len(matrix) - 1
    lim_col = len(matrix[0]) - 1
    if matrix[lim_line][lim_col] != 0:
        po = (
            sum([matrix[i][i] for i in range(0, min(lim_line, lim_col))])
            / matrix[lim_line][lim_col]
        )
        pe = sum(
            [
                matrix[lim_line][i] * matrix[i][lim_col]
                for i in range(0, min(lim_line, lim_col))
            ]
        ) / (matrix[lim_line][lim_col] ** 2)
        if 1 - pe != 0:
            kappa = (po - pe) / (1 - pe)
        else:
            kappa = float("nan")
    else:
        po = float("nan")
        pe = float("nan")
        kappa = float("nan")

    return {"Po": po, "Pe": pe, "kappa": kappa}


def estimate_kappa(matrix: list) -> dict:
    """
    Take a matrix, recalculate totals ignoring NA informations (before last one), and
    return a dictionnary containing Pa, PC and kappa estimated on this partial matrix.
    """
    lim_line = len(matrix) - 1
    lim_col = len(matrix[0]) - 1
    for i in range(0, lim_line - 1):
        matrix[i][lim_col] = sum(matrix[i][:-2])
    for i in range(0, lim_col - 1):
        matrix[lim_line][i] = sum([line[i] for line in matrix[:-2]])
    # Grand total:
    matrix[lim_line][lim_col] = sum(matrix[lim_line][:-2])
    if matrix[lim_line][lim_col] != 0:
        po = (
            sum([matrix[i][i] for i in range(0, min(lim_line, lim_col) - 1)])
            / matrix[lim_line][lim_col]
        )
        pe = sum(
            [
                matrix[lim_line][i] * matrix[i][lim_col]
                for i in range(0, min(lim_line, lim_col) - 1)
            ]
        ) / (matrix[lim_line][lim_col] ** 2)
        if 1 - pe != 0:
            kappa = (po - pe) / (1 - pe)
        else:
            kappa = float("nan")
    else:
        po = float("nan")
        pe = float("nan")
        kappa = float("nan")

    return {"Po": po, "Pe": pe, "kappa": kappa}
