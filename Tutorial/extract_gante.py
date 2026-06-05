import dendropy
import sys
import argparse


parser = argparse.ArgumentParser(
    description="Extract Gante et al. Nexus gene trees into a DAFT-readable TXT file."
)

parser.add_argument(
    "path",
    help="Path to DAFT_utils folder, or parent folder. Example: ../"
)

parser.add_argument(
    "input_file",
    help="Input Nexus tree file. Example: raxml.trees"
)

parser.add_argument(
    "output_file",
    help="Output TXT file. Example: Gante.txt"
)

args = parser.parse_args()


sys.path.append(args.path + "/DAFT_utils")
sys.path.append(args.path + "/DAFT_utils/reconcILS")

from reconcILS import *
from utils_reconcILS import *


trees = dendropy.TreeList.get(
    path=args.input_file,
    schema="nexus"
)


gene_trees = []

for ili, trea in enumerate(trees, 1):

    reco = reconcils()
    red = readWrite.readWrite()

    tree = str(trea.as_string(schema="newick"))
    tree = tree.replace('e-', '0')

    if len(tree) == 0:
        continue

    tr = red.parse_bio(red.to_newick(red.parse_bio(tree)))

    gene_trees.append(tr.to_newick())


with open(args.output_file, "w") as f:
    for gt in gene_trees:
        f.write(gt.strip() + "\n")