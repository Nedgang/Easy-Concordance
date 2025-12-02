##########
# IMPORT #
##########
from bin.utils.kappa_utils import estimate_kappa_na, estimate_kappa
from bin.utils.matrix_utils import compute_matrix_totals

import itertools
import pysam


########
# MAIN #
########
def compute_samples_concordance(
    geno_path: str,
    seq_path: str,
    output_dir: str,
    probes_limit,
    limit_chr,
    list_samples: list,
    exclude_chr: list,
):
    # Input files
    file_geno = pysam.VariantFile(geno_path)
    file_seq = pysam.VariantFile(seq_path)
    # Check if there is a sample list, otherwise use genotyping file samples
    if list_samples == []:
        list_samples = [sample for sample in file_geno.header.samples]
    # For each sample, a matrix
    samples_matrix = create_return_dicts(list_samples)
    i = 0
    for genotypage, sequencage in zip(file_geno.fetch(), file_seq.fetch()):
        i += 1
        if probes_limit is not None and i > probes_limit:
            break
        elif limit_chr is not None and genotypage.chrom == limit_chr:
            break
        elif (
            genotypage.chrom == sequencage.chrom
            and genotypage.pos == sequencage.pos
            and genotypage.ref == sequencage.ref
            and genotypage.alts == sequencage.alts
        ):
            if genotypage.chrom not in exclude_chr:
                samples_matrix = fill_individual_matrix(
                    genotypage, sequencage, list_samples, samples_matrix
                )
        else:
            print("Error, genotyping and sequencing not sorted identically.")
            print(
                "Genotyping probe:",
                genotypage.chrom,
                genotypage.pos,
                genotypage.ref,
                genotypage.alts,
            )
            print(
                "Sequencing probe:",
                sequencage.chrom,
                sequencage.pos,
                sequencage.ref,
                sequencage.alts,
            )
            raise ValueError("Genotyping and sequencing not in the same order.")

    # Output file
    output_file = open(output_dir + "/concordance_samples.tsv", "w")
    # Header
    output_file.write(
        "indiv\t"
        + "N00\tN01\tN02\tN0.\tT0.\t"
        + "N10\tN11\tN12\tN1.\tT1.\t"
        + "N20\tN21\tN22\tN2.\tT2.\t"
        + "N.0\tN.1\tN.2\tN..\tTna.\t"
        + "T.0\tT.1\tT.2\tT.na\tTM\t"
        + "Po_na\tPe_na\tkappa_na\tPo\tPe\tkappa\n"
    )
    for sample in samples_matrix:
        output_file.write(write_samples_output_line(sample, samples_matrix[sample]))
    output_file.close()
    file_geno.close()
    file_seq.close()


#############
# FUNCTIONS #
#############
def create_return_dicts(list_sample_geno: list) -> dict:
    """
    Return a dict of empty matrix ready to be supplied for each sample id, and an
    a dict with a counter for each heterozygote sequencing (difficult to track with
    matrix in case of multiall)
    """
    result_dict = {}
    # All 5x8 matrix emplacement, line by line.
    # Genotypage 00, 01, 11, 02, 12, 22, NA, Total
    for sample_id in list_sample_geno:
        # Genotypage 00, 01, 11, NA, Total
        result_dict[sample_id] = [
            [0, 0, 0, 0, 0],  # Sequençage 00
            [0, 0, 0, 0, 0],  # 01
            [0, 0, 0, 0, 0],  # 11
            [0, 0, 0, 0, 0],  # NA
            [0, 0, 0, 0, 0],
        ]  # Total
    return result_dict


def fill_individual_matrix(
    genotypage: pysam.libcbcf.VariantRecord,
    sequencage: pysam.libcbcf.VariantRecord,
    list_sample_geno: list,
    samples_matrix: dict,
) -> dict:
    """
    For each sample, increment in the samples_matrix the matrix emplacement corresponding
    to the concordance between genotyping and sequencing.
    Increment also if needed the sequencing heterozygote counter.
    Return the updated dicts.
    """
    for sample in list_sample_geno:
        sample_genotypage = dict(genotypage.samples[sample])["GT"]
        sample_sequencage = dict(sequencage.samples[sample])["GT"]
        matrix_sample = samples_matrix[sample]
        if None not in sample_sequencage:
            seq_a = sample_sequencage[0]
            seq_b = sample_sequencage[1]
            coord_line = seq_a + seq_b
        else:
            coord_line = -2
        if sample_genotypage[0] is not None:
            coord_col = sample_genotypage[0] + sample_genotypage[1]
        else:
            coord_col = -2
        # When correct place as been determined, increment it.
        matrix_sample[coord_line][coord_col] += 1
    return samples_matrix


def write_samples_output_line(sample: str, matrix: list) -> str:
    """
    Format the return line with all the infos.
    """
    matrix = compute_matrix_totals(matrix)
    # Estimating kappa_na
    components_kappa_na = estimate_kappa_na(matrix)
    # Estimating kappa
    components_kappa = estimate_kappa(matrix)
    return (
        sample
        + "\t"
        + "\t".join([str(val) for val in itertools.chain.from_iterable(matrix)])
        + "\t"
        + "\t".join(
            (
                str(components_kappa_na["Po"]),
                str(components_kappa_na["Pe"]),
                str(components_kappa_na["kappa"]),
                str(components_kappa["Po"]),
                str(components_kappa["Pe"]),
                str(components_kappa["kappa"]),
            )
        )
        + "\n"
    )
