# Easy_Concordance 🧬 💻
[![DOI](https://zenodo.org/badge/1108490133.svg)](https://doi.org/10.5281/zenodo.17865414)

## What is it?
Python tool to estimate concordance between genotyping and sequencing using Cohen's Kappa.
Two kappa are computed, one which take in account missing value (kappa_na), the other only when both genotyping and sequencing are known (kappa).

## Requirements
Python version 3.11+

Python packages:
- matplotlib
- polars
- pysam
- seaborn

## Install
```bash
  pip install matplotlib polars pysam seaborn
  git clone git@github.com:Nedgang/Easy-Concordance.git
```

## Input files format
- Genotyping and sequencing must be in VCF format (compressed bcf and vcf.gz also work).
- Each probe must be found in both files (see pre-processing section for more details).
- Both files must be sorted in the same way (see pre-processing section for more details).
- Samples will be read from genotyping file. Those must be present into sequencing file. You can have other samples in the sequencing file, they won't be evaluated. It is however strongly advised to remove them to reduce computation time.
- Samples does not need to be on the same order in the two files. They will be found via their ID read in the genotyping file.

## Usage
### Main actions
#### Concordance per probes
Estimate concordance for each probes, and product distribution and repartition on the genome figures.
```bash
  python easy_concordance.py probes -g GENOTYPING_FILE -s SEQUENCING_FILE [-i list_interest_samples] [-l limit_chr] [-e exclude_chr] [-t]
```
Options:
- -e, --exclude_chr: Chromosome to exclude from computation. Possible to pass multiple, separated by comma. Ex: chr12,chr13
- -i, --list_interest_samples: Text file containing the list of samples to use for concordance computation. Format: a sample per line.
- -l, --limit_chr: Chromosome at which to stop compute concordance. This chromosome, and the ones after it, won't be taken in account in the calculations. Faster than exclude_chr.
- -o, --output_dir: Path to directory to fill with results (default= working directory).
- -t, --test: Test run on 1000 probes.

#### Concordance per samples
Estimate concordance for each sample, and product distribution, comparison with heterozygosis and with % of missing genotyping/sequencing figures.
```bash
  python easy_concordance.py samples -g GENOTYPING_FILE -s SEQUENCING_FILE [-i list_interest_samples] [-l limit_chr] [-e exclude_chr] [-t]
```
Options:
- -e, --exclude_chr: Chromosome to exclude from computation. Possible to pass multiple, separated by comma. Ex: chr12,chr13
- -i, --list_interest_samples: Text file containing the list of samples to use for concordance computation. Format: a sample per line.
- -l, --limit_chr: Chromosome at which to stop compute concordance. This chromosome, and the ones after it, won't be taken in account in the calculations. Faster than exclude_chr.
- -o, --output_dir: Path to directory to fill with results (default= working directory).
- -t, --test: Test run on 1000 probes.

### Sub-actions
Those actions allow, in case of need, to redraw the figures without having to compute concordance once again.
#### Draw probes
Draw distribution and repartition figures for probes concordance.
```bash
  python easy_concordance.py draw_probes -c CONCORDANCE_DIR [-o output_dir]
```
Options:
- -o, --output_dir: Path to directory to fill with graphs (default= working directory).

#### Draw samples
Draw distribution, comparison with heterozygosis and with % of missing genotyping/sequencing figures for samples concordance.
```bash
  python easy_concordance.py draw_samples -c CONCORDANCE_DIR [-o output_dir]
```
Options:
- -o, --output_dir: Path to directory to fill with graphs (default= working directory).

### Examples
- Concordance per samples
```bash
  python easy_concordance/easy_concordance.py samples \
    -g genotyping_norm_sorted_filtered.bcf \
    -s sequencing_norm_sorted_filtered.bcf \
    -o concordance
```

- Concordance per probes
```bash
  python easy_concordance/easy_concordance.py probes \
    -g genotyping_norm_sorted_filtered.bcf \
    -s sequencing_norm_sorted_filtered.bcf \
    -o concordance
```

- Concordance per samples on autosomes
```bash
  python easy_concordance/easy_concordance.py samples \
    -g genotyping_norm_sorted_filtered.bcf \
    -s sequencing_norm_sorted_filtered.bcf \
    -l chrX \
    -o concordance_autosomes
```

## Concordance estimation
### Comparison matrix
The concordance estimation is based on a Cohen's Kappa.
To compute it, a comparison matrix is needed. In the case of an estimation of concordance per probes, there will be a matrix per probes, which will be populated with the comparison of the sequencing and genotyping data of the probe from each sample. In the case of a concordance per samples, each sample will be associated with a matrix, in which the comparison genotyping/sequencing of all the probes of this sample will be recorded.

| | 00 | 01 | 11 | NA | Total |
| --- | --- | --- | --- | --- | --- |
| 00 | N00 | N01 | N02 | N0. | T0. |
| 01 | N10 | N11 | N12 | N1. | T1. |
| 11 | N20 | N21 | N22 | N2. | T2. |
| NA | N.0 | N.1 | N.2 | N.. | Tna. |
| Total | T.0 | T.1 | T.2 | T.na | TM |

In columns, genotyping information, in lines, sequencing information. TM is the total number of comparison made in the matrix.
When comparing genotyping and sequencing probe, the count of the corresponding emplacement in the matrix is increased by one. Once the matrix is fully populated, the total of each line and column is made, as well as the grand total of the matrix.

Ex: For a given probe of a given sample, with a genotyping of 0/1 and a sequencing of 1/1, we will increment by 1 the emplacement N21. In an other case with a genotyping of 0/0 and a NA sequencing, the increased emplacement will be N.0. 

### Parameters
Once the matrix is full, we can compute the parameters and Cohen's kappa.

$$p_{o} = \frac{\sum_{i}{N_{ii}}}{TM}$$

$$p_{e} = \frac{\sum_{i}{T_{i.}T_{.i}}}{TM^2}$$

$$kappa = \frac{p_{o} - p_{e}}{1 - p_{e}}$$

One of them is calculated with all the data, including missing ones (kappa_na), the other excluding them (kappa). When excluding missing data, the totals for each line, column and the whole matrix are of course updated.

## Output files
The results directory will contain 2 sub-directories
### Figures
#### Concordance per probes
- probes_distribution_kappa.png: Distribution of the probes kappa.
- probes_distribution_kappa_na.png: Distribution of the probes kappa_na.
- probes_kappa_repartition_window_11snp.png: Distribution of the probes kappa on the genome, using a slidding window of 11 probes.
- probes_kappa_na_repartition_window_11snp.png: Distribution of the probes kappa_na on the genome, using a slidding window of 11 probes.

#### Concordance per samples
- samples_%_missing_gen_vs_kappa.png: Plotting of the comparison between percentage of missing genotyping and kappa score.
- samples_%_missing_gen_vs_kappa_na.png: Plotting of the comparison between percentage of missing genotyping and kappa_na score.
- samples_%_missing_seq_vs_kappa.png: Plotting of the comparison between percentage of missing sequencing and kappa score.
- samples_%_missing_seq_vs_kappa_na.png: Plotting of the comparison between percentage of missing sequencing and kappa_na score.
- samples_distribution_kappa.png: Distribution of the kappa score.
- samples_distribution_kappa_na.png: Distribution of the kappa_na score.
- samples_gen_heterozygosis_vs_kappa.png: Plotting of the comparison between proportion of genotyping heterozygosis and kappa score.
- samples_gen_heterozygosis_vs_kappa_na.png: Plotting of the comparison between proportion of genotyping heterozygosis and kappa_na score.
- samples_seq_heterozygosis_vs_kappa.png: Plotting of the comparison between proportion of sequencing heterozygosis and kappa score.
- samples_seq_heterozygosis_vs_kappa_na.png: Plotting of the comparison between proportion of sequencing heterozygosis and kappa_na score.

### Results
- concordance_probes.tsv: Concordance results per probes.
- concordance_samples.tsv: Concordance results per samples.

## Automatic checks
This tool only check that it compare the right variant and sample, otherwise there is no check nor correction on the quality of the data used for computation.
This is not optimised for that, and would not be as efficient as dedicated tool for that purpose.

## Pre-processing sequencing and genotyping files
The final genotyping and sequencing file must both be normalized into biallelic records, contain the same sites and samples, and be in the same order.
It also is recommanded to remove any QUAL, FILTER, INFO or FORMAT (besides Genotype) from the files. Those information are not needed for computation, and removing them, especially in the case of numerous sample, can significanly reduce both the files sizes and computation times.

Here are some command example using bcftools for formatting your data to the correct format.
### Annotation cleaning
```bash
  bcftools annotate\
    -x QUAL,FILTER,INFO,FORMAT\
    --threads 12\
    {genotypes.bcf}
```

### Normalization
```bash
  bcftools norm \
    -c s\
    -f {reference.fasta} \
    --threads {n_threads}\
    -O b\
    -o tmp.bcf\
    {raw_genotypes.bcf}

  bcftools index -f tmp.bcf

  bcftools sort \
    -m 45G \
    -O b \
    -o {genotypes_norm_sorted.bcf} \
    tmp.bcf

  bcftools index -f {genotypes_norm_sorted.bcf}
```

### Filtering
```bash
  bcftools annotate \
    -a {sequencing_norm_sorted_filtered.bcf} \
    -c CHROM,POS,REF,ALT \
    --regions-overlap pos \
    --mark-sites +Sequenced \
    -h <(echo '##INFO=<ID=Sequenced,Number=0,Type=Flag,Description=Sequenced>') \
    --threads 2 \
    {genotypes_norm_sorted.bcf} \
    | bcftools view \
      -e INFO/Sequenced=0 \
      --threads 2 \
      -O b \
      -o tmp.bcf

  bcftools index -f tmp.bcf

  bcftools sort\
    -O b \
    -o {genotypes_norm_sorted_filtered.bcf} \
    tmp.bcf
```
