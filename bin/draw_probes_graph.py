##########
# IMPORT #
##########
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns

from bin.utils.graph_utils import moving_window


########
# MAIN #
########
def draw_probes_graph(concordance_file_path: str, output_dir: str):
    """
    Draw distribution and manhattan plot for kappa and kappa_na in output_dir
    """
    dataframe = pl.read_csv(concordance_file_path, separator="\t")
    dataframe = dataframe.with_columns(pl.col("kappa").cast(pl.Float64).alias("kappa"))
    dataframe = dataframe.with_columns(
        pl.col("kappa_na").cast(pl.Float64).alias("kappa_na")
    )
    displot = sns.displot(dataframe["kappa"], bins=50)
    plt.savefig(output_dir + "/probes_distribution_kappa.png")
    plt.close()

    displot = sns.displot(dataframe["kappa_na"], bins=50)
    plt.savefig(output_dir + "/probes_distribution_kappa_na.png")
    plt.close()

    list_chrom = list(dataframe["chrom"].unique(maintain_order=True))
    new_dt = dataframe.filter(pl.col("chrom") == list_chrom[0])
    if len(list_chrom) > 1:
        for target_chr in list_chrom[1:]:
            cumulative_pos = new_dt["pos"].max()
            tmp_dt = dataframe.filter(pl.col("chrom") == target_chr).with_columns(
                (pl.col("pos") + cumulative_pos).alias("pos")
            )
            new_dt = pl.concat([new_dt, tmp_dt], how="vertical")

    window_dt = moving_window(new_dt, 11)
    scatterplot = sns.relplot(
        data=window_dt,
        x="pos",
        y="window_kappa",
        hue="chrom",
        palette="Set1",
        aspect=4,
        linewidth=0,
        s=4,
        legend=None,
    )
    scatterplot.ax.set_xlabel("Chromosome")
    scatterplot.ax.set_xticks([])
    plt.savefig(output_dir + "/probes_kappa_repartition_window_11snp.png")
    plt.close()

    scatterplot = sns.relplot(
        data=window_dt,
        x="pos",
        y="window_kappa_na",
        hue="chrom",
        palette="Set1",
        aspect=4,
        linewidth=0,
        s=4,
        legend=None,
    )
    scatterplot.ax.set_xlabel("Chromosome")
    scatterplot.ax.set_xticks([])
    plt.savefig(output_dir + "/probes_kappa_na_repartition_window_11snp.png")
