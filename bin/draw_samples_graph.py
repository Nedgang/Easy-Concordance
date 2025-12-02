##########
# IMPORT #
##########
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns


########
# MAIN #
########
def draw_samples_graph(concordance_file_path: str, output_dir: str):
    """
    Draw distribution for kappa and kappa_na in output_dir, and link between kappa/kappa_na
    and missing seq/genotyping, and seq/genotype heterozygosis.
    """
    dataframe = pl.read_csv(concordance_file_path, separator="\t")
    displot = sns.displot(dataframe["kappa"], bins=50)
    plt.savefig(output_dir + "/samples_distribution_kappa.png")
    plt.close()
    displot = sns.displot(dataframe["kappa_na"], bins=50)
    plt.savefig(output_dir + "/samples_distribution_kappa_na.png")
    plt.close()

    dataframe = dataframe.with_columns(
        (
            pl.col("T1.")
            / (pl.col("TM") - (pl.col("T.na") + pl.col("Tna.") - pl.col("N..")))
        ).alias("seq_heterozygosis"),
        (
            pl.col("T.1")
            / (pl.col("TM") - (pl.col("T.na") + pl.col("Tna.") - pl.col("N..")))
        ).alias("gen_heterozygosis"),
        (pl.col("Tna.") / pl.col("TM")).alias("%_missing_seq"),
        (pl.col("T.na") / pl.col("TM")).alias("%_missing_gen"),
    )

    scatterplot = sns.relplot(data=dataframe, x="seq_heterozygosis", y="kappa")
    plt.savefig(output_dir + "/samples_seq_heterozygosis_vs_kappa.png")
    plt.close()

    scatterplot = sns.relplot(data=dataframe, x="seq_heterozygosis", y="kappa_na")
    plt.savefig(output_dir + "/samples_seq_heterozygosis_vs_kappa_na.png")
    plt.close()

    scatterplot = sns.relplot(data=dataframe, x="gen_heterozygosis", y="kappa")
    plt.savefig(output_dir + "/samples_gen_heterozygosis_vs_kappa.png")
    plt.close()

    scatterplot = sns.relplot(data=dataframe, x="gen_heterozygosis", y="kappa_na")
    plt.savefig(output_dir + "/samples_gen_heterozygosis_vs_kappa_na.png")
    plt.close()

    scatterplot = sns.relplot(data=dataframe, x="%_missing_seq", y="kappa")
    plt.savefig(output_dir + "/samples_%_missing_seq_vs_kappa.png")
    plt.close()

    scatterplot = sns.relplot(data=dataframe, x="%_missing_seq", y="kappa_na")
    plt.savefig(output_dir + "/samples_%_missing_seq_vs_kappa_na.png")
    plt.close()

    scatterplot = sns.relplot(data=dataframe, x="%_missing_gen", y="kappa")
    plt.savefig(output_dir + "/samples_%_missing_gen_vs_kappa.png")
    plt.close()

    scatterplot = sns.relplot(data=dataframe, x="%_missing_gen", y="kappa_na")
    plt.savefig(output_dir + "/samples_%_missing_gen_vs_kappa_na.png")
    plt.close()
