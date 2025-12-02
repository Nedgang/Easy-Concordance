##########
# IMPORT #
##########
import polars as pl


#############
# FUNCTIONS #
#############
def moving_window(dataframe, window_size):
    """
    Take a dataframe, and return a version of it with kappa computed on moving frame.
    """
    if window_size % 2 == 0:
        raise ValueError(
            "Window size must be an uneven number to be centered on one SNP"
        )
    else:
        lim_window = (window_size - 1) // 2
    dataframe = dataframe.with_row_index().with_columns(
        window_kappa=0, window_kappa_na=0
    )
    for index in range(lim_window, len(dataframe) - lim_window):
        dataframe = dataframe.with_columns(
            window_kappa=pl.when(pl.col("index") == index)
            .then(
                dataframe[index - lim_window : index + lim_window]["kappa"]
                .drop_nans()
                .mean()
            )
            .otherwise(pl.col("window_kappa")),
            window_kappa_na=pl.when(pl.col("index") == index)
            .then(
                dataframe[index - lim_window : index + lim_window]["kappa_na"]
                .drop_nans()
                .mean()
            )
            .otherwise(pl.col("window_kappa_na")),
        )

    return dataframe[5:-5]
