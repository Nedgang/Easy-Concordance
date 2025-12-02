##########
# IMPORT #
##########
from bin.utils.kappa_utils import estimate_kappa_na, estimate_kappa
from bin.utils.matrix_utils import compute_probe_matrix

import itertools
import pysam


########
# MAIN #
########
def compute_probes_concordance(
    geno_path: str,
    seq_path: str,
    output_dir: str,
    probes_limit: int,
    limit_chr: str,
    list_samples: list,
    exclude_chr: list,
):
    # Input files
    file_geno = pysam.VariantFile(geno_path)
    file_seq = pysam.VariantFile(seq_path)
    # Output files
    output_file = open(output_dir + "/concordance_probes.tsv", "w")
    # Check if there is a sample list, otherwise use genotyping file samples
    if list_samples == []:
        list_samples = [sample for sample in file_geno.header.samples]
    # Writing output header
    output_file.write(
        "chrom\tpos\tid\tref\talt\t"
        + "N00\tN01\tN02\tN0.\tT0.\t"
        + "N10\tN11\tN12\tN1.\tT1.\t"
        + "N20\tN21\tN22\tN2.\tT2.\t"
        + "N.0\tN.1\tN.2\tN..\tTna.\t"
        + "T.0\tT.1\tT.2\tT.na\tTM\t"
        + "Po_na\tPe_na\tkappa_na\tPo\tPe\tkappa\n"
    )
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
                output_file.write(
                    write_probes_output_line(genotypage, sequencage, list_samples)
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

    output_file.close()
    file_geno.close()
    file_seq.close()


#############
# FUNCTIONS #
#############
def write_probes_output_line(
    genotypage: pysam.libcbcf.VariantRecord,
    sequencage: pysam.libcbcf.VariantRecord,
    list_samples: list,
) -> str:
    """
    Format the return line with all the infos, and return it.
    """
    # Writing informations
    return_line = "\t".join(
        (
            str(genotypage.chrom),
            str(genotypage.pos),
            str(genotypage.id),
            str(sequencage.ref),
            str(sequencage.alts),
        )
    )
    # Writing matrix
    concordance_matrix = compute_probe_matrix(genotypage, sequencage, list_samples)
    # Estimating kappa_na
    components_kappa_na = estimate_kappa_na(concordance_matrix)
    # Estimating kappa
    components_kappa = estimate_kappa(concordance_matrix)
    return (
        return_line
        + "\t"
        + "\t".join(
            [str(val) for val in itertools.chain.from_iterable(concordance_matrix)]
        )
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
