#!/usr/bin/python

##########
# IMPORT #
##########
import argparse
import os

from bin.compute_probes_concordance import compute_probes_concordance
from bin.compute_samples_concordance import compute_samples_concordance
from bin.draw_probes_graph import draw_probes_graph
from bin.draw_samples_graph import draw_samples_graph

##########
# PARSER #
##########
parser = argparse.ArgumentParser(prog="python easy_concordance.py")
# Creation of the subparsers for the different actions
subparsers = parser.add_subparsers(
    title="Actions",
    dest="action",
    required=True,
    help="o Probes: estimate concordance per probes, and draw graphs. o Samples: estimate concordance per samples, and draw graphs. o Draw_probes: draw graphs for concordance by probes. o Draw_samples: draw graphs for concordance by samples.",
)
# Parser for the concordance estimation per probes
probes_parser = subparsers.add_parser("probes")
probes_parser.add_argument(
    "-g",
    "--genotyping",
    dest="genotyping_file",
    required=True,
    help="Path to genotyping VCF file.",
    type=str,
)
probes_parser.add_argument(
    "-s",
    "--sequencing",
    dest="sequencing_file",
    required=True,
    help="Path to sequencing VCF file.",
    type=str,
)
probes_parser.add_argument(
    "-o",
    "--output_dir",
    help="Path to directory to fill with results (default = \
                           working directory).",
    default=".",
    type=str,
)
probes_parser.add_argument(
    "-i",
    "--list_interest_samples",
    help="Text file containing the list of samples to use for concordance computation. Format: a sample per line.",
    type=str,
)
probes_parser.add_argument(
    "-l",
    "--limit_chr",
    help="Chromosome at which to stop compute concordance. This chromosome, and the ones after it, won't be taken in account in the calculations. Faster than exclude_chr.",
    type=str,
)
probes_parser.add_argument(
    "-e",
    "--exclude_chr",
    help="Chromosome to exclude from computation. Possible to pass multiple, separated by comma. Ex: chr12,chr13",
    default=[],
    type=str,
)
probes_parser.add_argument(
    "-t",
    "--test",
    dest="probes_limit",
    help="Test run on 1000 probes.",
    const=1000,
    action="store_const",
)
# Parser for the concordance estimation per samples
samples_parser = subparsers.add_parser("samples")
samples_parser.add_argument(
    "-g",
    "--genotyping",
    dest="genotyping_file",
    required=True,
    help="Path to genotyping VCF file.",
    type=str,
)
samples_parser.add_argument(
    "-s",
    "--sequencing",
    dest="sequencing_file",
    required=True,
    help="Path to sequencing VCF file.",
    type=str,
)
samples_parser.add_argument(
    "-o",
    "--output_dir",
    help="Path to directory to fill with results (default= working directory).",
    default=".",
    type=str,
)
samples_parser.add_argument(
    "-i",
    "--list_interest_samples",
    help="Text file containing the list of samples to use for concordance computation. Format: a sample per line.",
    type=str,
)
samples_parser.add_argument(
    "-l",
    "--limit_chr",
    help="Chromosome at which to stop compute concordance.\
                            This chromosome, and the ones after it, won't be taken in\
                            account in the calculations. Faster than exclude_chr.",
    type=str,
)
samples_parser.add_argument(
    "-e",
    "--exclude_chr",
    help="Chromosome to exclude from computation. Possible to\
                            pass multiple, separated by comma. Ex: chr12,chr13",
    default=[],
    type=str,
)
samples_parser.add_argument(
    "-t",
    "--test",
    dest="probes_limit",
    help="Test run on 1000 probes.",
    const=1000,
    action="store_const",
)
# Parser for drawing probes graph
draw_probes_parser = subparsers.add_parser("draw_probes")
draw_probes_parser.add_argument(
    "-c",
    "--concordance_dir",
    required=True,
    help="Path to dir with concordance file",
    type=str,
)
draw_probes_parser.add_argument(
    "-o",
    "--output_dir",
    help="Path to directory to fill with graphs (default =\
                                working directory).",
    default=".",
    type=str,
)
# Parser for drawing samples graph
draw_samples_parser = subparsers.add_parser("draw_samples")
draw_samples_parser.add_argument(
    "-c",
    "--concordance_dir",
    required=True,
    help="Path to dir with concordance file",
    type=str,
)
draw_samples_parser.add_argument(
    "-o",
    "--output_dir",
    help="Path to directory to fill with graphs (default \
                                 = working directory).",
    default=".",
    type=str,
)
# Instantiation
args = parser.parse_args()

##########
# SCRIPT #
##########
output_dir = args.output_dir
action = args.action

if action == "probes":
    if args.list_interest_samples is not None:
        with open(args.list_interest_samples) as file:
            list_samples = [sample.strip() for sample in file]
    else:
        list_samples = []
    os.system("mkdir -p " + output_dir)
    os.system("mkdir -p " + output_dir + "/results")
    os.system("mkdir -p " + output_dir + "/figures")
    compute_probes_concordance(
        args.genotyping_file,
        args.sequencing_file,
        output_dir + "/results/",
        args.probes_limit,
        args.limit_chr,
        list_samples,
        args.exclude_chr,
    )
    draw_probes_graph(
        output_dir + "/results/concordance_probes.tsv", output_dir + "/figures"
    )

elif action == "samples":
    if args.list_interest_samples is not None:
        with open(args.list_interest_samples) as file:
            list_samples = [sample.strip() for sample in file]
    else:
        list_samples = []
    os.system("mkdir -p " + output_dir)
    os.system("mkdir -p " + output_dir + "/results")
    os.system("mkdir -p " + output_dir + "/figures")
    compute_samples_concordance(
        args.genotyping_file,
        args.sequencing_file,
        output_dir + "/results/",
        args.probes_limit,
        args.limit_chr,
        list_samples,
        args.exclude_chr,
    )
    draw_samples_graph(
        output_dir + "/results/concordance_samples.tsv", output_dir + "/figures"
    )

elif action == "draw_probes":
    draw_probes_graph(args.concordance_dir + "/concordance_probes.tsv", output_dir)

elif action == "draw_samples":
    draw_samples_graph(args.concordance_dir + "/concordance_samples.tsv", output_dir)

else:
    raise ValueError("Action asked not available:", action)
