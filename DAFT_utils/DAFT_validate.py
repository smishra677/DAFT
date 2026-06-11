import sys
import csv
import hashlib
import numpy as np

sys.path.append("./DAFT_utils/reconcILS")
sys.path.append("./DAFT_utils/")

from utils_reconcILS import *
from reconcILS import *
from pathlib import Path
from collections import Counter

from ete3 import Tree


reco = reconcils()
red = readWrite.readWrite()
Il = ILS.ILS()


try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    RICH_AVAILABLE = True
except Exception:
    RICH_AVAILABLE = False


class daft_validate:
    def __init__(self):
        self.species_tree_newick = None
        self.gene_treefile = None
        self.gene_tree_strings = None
        self.input_hash = None
        self.taxon_report_path = None
        self.rooted_gene_treefile = None
        self.console=None
        self.log=None
        self.verbose=None


    def log_message(self,message, log):
        if self.verbose!=0:
            print(message, flush=True)
        
        print(message, file=log, flush=True)


    def log_only(self,message, log):
        print(message, file=log, flush=True)


    def log_header(self,title, log):
        self.log_message("", log)
        self.log_message("=" * 80, log)
        self.log_message(title, log)
        self.log_message("-" * 80, log)


    def log_table_only(self,title, columns, rows, log):
        self.log_only("", log)
        self.log_only("=" * 80, log)
        self.log_only(title, log)
        self.log_only("-" * 80, log)
        self.log_only("\t".join([str(x) for x in columns]), log)

        for row in rows:
            self.log_only("\t".join([str(x) for x in row]), log)


    def print_panel(self,title, message, log, console=None, style="green"):
        if RICH_AVAILABLE and console:
            console.print(Panel(message, title=title, border_style=style))

            self.log_only("", log)
            self.log_only("=" * 80, log)
            self.log_only(title, log)
            self.log_only("-" * 80, log)
            self.log_only(message, log)

        else:
            self.log_header(title, log)
            self.log_message(message, log)
        self.log_only("*" * 80, log)


    def print_table(self,title, columns, rows, log, console=None):
        if RICH_AVAILABLE and console:
            table = Table(title=title)

            for col in columns:
                table.add_column(str(col))

            for row in rows:
                row_style = None

                if len(row) >= 2:
                    status = str(row[1]).upper()

                    if status == "OK":
                        row_style = "green"
                    elif status == "WARNING":
                        row_style = "yellow"
                    elif status == "FAIL" or status == "FAILED" or status == "ERROR":
                        row_style = "red"

                if row_style:
                    table.add_row(*[str(x) for x in row], style=row_style)
                else:
                    table.add_row(*[str(x) for x in row])

            console.print(table)
            self.log_table_only(title, columns, rows, log)

        else:
            self.log_header(title, log)
            self.log_message("\t".join([str(x) for x in columns]), log)

            for row in rows:
                self.log_message("\t".join([str(x) for x in row]), log)
        self.log_only("*" * 80, log)
        self.log_only("*" * 80, log)
        


    def print_startup_report(self,process,parser, species_tree_source, gene_treefile, output_prefix, log, console=None):
        log_path = output_prefix + f"_{process}_log.txt"
        #print(log_path)
        #exit()

        self.print_panel(
            process,
            f"Starting new {process} run\n\n" +
            "Output prefix: " + str(output_prefix) + "\n" +
            "Log file: " + str(Path(log_path).resolve()),
            log,
            console,
            "green",
        )

        rows = []

        for key, value in vars(parser).items():
            rows.append([key, value])

        self.print_table(
            "Run options",
            ["Option", "Value"],
            rows,
            log,
            console,
        )

        gt_path = Path(gene_treefile)

        rows = [
            ["Species tree source", species_tree_source],
            ["Gene tree file", gt_path.resolve()],
            ["Gene tree file exists", gt_path.exists()],
            ["Gene tree file readable", gt_path.is_file()],
            ["Output prefix", output_prefix],
            ["DAFT log file", Path(log_path).resolve()],
        ]

        self.print_table(
            "Input files",
            ["Item", "Value"],
            rows,
            log,
            console,
        )


    def read_gene_tree_strings(self,gene_treefile):
        gene_tree_strings = []

        with open(gene_treefile, "r", encoding="utf-8", newline="") as f:
            reader = csv.reader(f)

            for row in reader:
                if not row:
                    continue

                gt_string = ",".join([x.strip() for x in row if x.strip()]).strip()

                if not gt_string:
                    continue

                if gt_string.lower() == "gt":
                    continue

                gene_tree_strings.append(gt_string)

        return gene_tree_strings


    def read_gene_tree_data(self,gene_treefile):
        gene_tree_strings = self.read_gene_tree_strings(gene_treefile)
        return np.array([[gt] for gt in gene_tree_strings])


    def check_newick_string(self,newick_string, label):
        newick_string = newick_string.strip()

        if not newick_string:
            raise ValueError(label + " is empty.")

        if not newick_string.endswith(";"):
            raise ValueError(label + " does not end with ';'.")

        if newick_string.count("(") != newick_string.count(")"):
            raise ValueError(label + " has unbalanced parentheses.")

        if newick_string.count(";") != 1:
            raise ValueError(label + " should contain exactly one semicolon.")


    def parse_tree(self,newick_string, label):
        self.check_newick_string(newick_string, label)

        try:
            tree = Tree(newick_string, format=1)
            tree_rc = red.parse(newick_string)
        except Exception as err:
            raise ValueError(label + " could not be parsed as Newick: " + str(err))

        if len(tree.get_leaf_names()) == 0:
            raise ValueError(label + " contains no taxa.")

        return tree,tree_rc


    def check_duplicate_taxa(self,tree, label,ignore_duplicate):
        taxa = tree.get_leaf_names()
        counts = Counter(taxa)
        duplicates = []

        for taxon, count in counts.items():
            if count > 1:
                duplicates.append(taxon)

        if len(duplicates)>0 and ignore_duplicate==0:
            raise ValueError(label + " contains duplicate taxa: " + ", ".join(sorted(duplicates)))
        elif len(duplicates)>0 and ignore_duplicate!=0:
            return False

        return taxa


    def check_binary_rooted(self,tree, label):
        if len(tree.children) != 2:
            raise ValueError(label + " root has " + str(len(tree.children)) + " children; expected 2.")

        for node in tree.traverse():
            if not node.is_leaf() and len(node.children) != 2:
                raise ValueError(label + " contains a polytomy/non-binary node with " + str(len(node.children)) + " children.")


    def root_split(self,tree):
        if len(tree.children) != 2:
            raise ValueError("Cannot calculate root split because the root does not have 2 children.")

        left = frozenset(tree.children[0].get_leaf_names())
        right = frozenset(tree.children[1].get_leaf_names())

        return frozenset([left, right])


    def restricted_root_split(self,species_split, gene_taxa):
        restricted = []

        for side in species_split:
            side_taxa = frozenset(side & gene_taxa)

            if side_taxa:
                restricted.append(side_taxa)

        return frozenset(restricted)


    def check_root_consistency(self,species_tree, gene_tree, label):
        species_split = self.root_split(species_tree)
        gene_split = self.root_split(gene_tree)

        gene_taxa = set(gene_tree.get_leaf_names())
        expected_split = self.restricted_root_split(species_split, gene_taxa)

        if len(expected_split) != 2:
            return False

        if gene_split != expected_split:
            raise ValueError(label + " is not rooted consistently with the species tree.")

        return True


    def parse_outgroup_arg(self,rooting):
        if rooting is None:
            return []

        return [x.strip() for x in rooting.split(",") if x.strip()]


    def root_tree_on_outgroup(self,tree, outgroup_taxa, label):
        if not outgroup_taxa:
            raise ValueError("No outgroup taxa were provided for " + label + ".")

        taxa = set(tree.get_leaf_names())
        missing = []

        for taxon in outgroup_taxa:
            if taxon not in taxa:
                missing.append(taxon)

        if missing:
            raise ValueError(label + " is missing requested outgroup taxa: " + ", ".join(sorted(missing)))

        if len(outgroup_taxa) == 1:
            outgroup_node = tree & outgroup_taxa[0]
        else:
            outgroup_node = tree.get_common_ancestor(outgroup_taxa)

            found = set(outgroup_node.get_leaf_names())
            requested = set(outgroup_taxa)

            if found != requested:
                raise ValueError(label + ": requested outgroup taxa are not monophyletic.")

        tree.set_outgroup(outgroup_node)

        return tree


    def write_rooted_gene_trees(self,gene_tree_strings, rooted_file):
        with open(rooted_file, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["gt"])

            for gt in gene_tree_strings:
                writer.writerow([gt])


    def check_gene_taxa_against_species(self,gene_taxa, species_taxa, label):
        extra_taxa = set(gene_taxa) - set(species_taxa)
        missing_taxa = set(species_taxa) - set(gene_taxa)

        if extra_taxa:
            raise ValueError(label + " contains taxa not in the species tree: " + ", ".join(sorted(extra_taxa)))

        return missing_taxa


    def build_taxon_report(self,species_taxa, gene_tree_taxa_counts):
        all_taxa = sorted(species_taxa | set(gene_tree_taxa_counts.keys()))
        rows = []

        for taxon in all_taxa:
            rows.append({
                "original_taxon": taxon,
                "in_species_tree": taxon in species_taxa,
                "gene_tree_count": gene_tree_taxa_counts.get(taxon, 0),
                "daft_internal_name": taxon,
            })

        return rows


    def write_taxon_report(self,rows, output_prefix):
        taxon_report_path = output_prefix + "_taxon_report.csv"

        with open(taxon_report_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "original_taxon",
                    "in_species_tree",
                    "gene_tree_count",
                    "daft_internal_name",
                ],
            )

            writer.writeheader()
            writer.writerows(rows)

        return taxon_report_path


    def check_internal_name_collisions(self,rows):
        internal_to_original = {}

        for row in rows:
            original = row["original_taxon"]
            internal = row["daft_internal_name"]

            if internal not in internal_to_original:
                internal_to_original[internal] = []

            internal_to_original[internal].append(original)

        collisions = {}

        for internal, originals in internal_to_original.items():
            if len(originals) > 1:
                collisions[internal] = originals

        if collisions:
            details = []

            for internal, originals in collisions.items():
                details.append(internal + ": " + ", ".join(originals))

            raise ValueError("DAFT internal taxon-name collision detected: " + "; ".join(details))


    def compute_input_hash(self,species_tree_newick, gene_tree_strings):
        h = hashlib.sha256()

        for gt in gene_tree_strings:
            h.update(gt.strip().encode("utf-8"))
            h.update(b"\n")

        h.update(species_tree_newick.strip().encode("utf-8"))
        h.update(b"\n")

        return h.hexdigest()


    def print_tree_summary(self,species_taxa, gene_tree_strings, species_root_split, log, console=None):
        root_sides = list(species_root_split)
        rows = []

        rows.append(["Species tree taxa", len(species_taxa)])
        rows.append(["Gene trees", len(gene_tree_strings)])

        if len(root_sides) == 2:
            side1 = ", ".join(sorted(root_sides[0]))
            side2 = ", ".join(sorted(root_sides[1]))

            rows.append(["Species root clade 1", side1])
            rows.append(["Species root clade 2", side2])

            if len(root_sides[0]) <= len(root_sides[1]):
                outgroup_guess = side1
            else:
                outgroup_guess = side2

            rows.append(["Assumed outgroup candidate", outgroup_guess])

        self.print_table(
            "Tree summary",
            ["Item", "Value"],
            rows,
            log,
            console,
        )


    def print_validation_summary(self,validation_rows, log, console=None):
        self.print_table(
            "Validation checks",
            ["Check", "Status", "Details"],
            validation_rows,
            log,
            console,
        )


    def validateme(self,process,sp, gene_treefile, output_prefix, parser,  rooting=None, forced=0, allow_inconsistent_rooting=0, ignore_duplicate=0,verbose=0,hashed=0):
        log_path = output_prefix + f"_{process}_log.txt"
        
        
        log = None
        self.verbose=verbose

        if verbose!=0:
            if RICH_AVAILABLE:
                console = Console()
                self.console=console
            else:
                console = None
        else:
            console = None

        try:
            log = open(log_path, "w", encoding="utf-8")
            self.log=log
            if hashed!=0:
                self.print_panel(
                    "djiNNI message",
                    "input validation already done in DAFT .\n"
                    f"Please check {output_prefix}_DAFT_log.txt for input validation report",
                    log,
                    console,
                    "yellow",
                )
                return self
            

            self.print_startup_report(
                process=process,
                parser=parser,
                species_tree_source="inline --sp argument",
                gene_treefile=gene_treefile,
                output_prefix=output_prefix,
                log=log,
                console=console,
            )

            validation_rows = []

            outgroup_taxa = self.parse_outgroup_arg(rooting)

            if rooting and not forced:
                self.print_panel(
                    "DAFT warning",
                    "--rooting was provided, but --forced 1 was not used.\n"
                    "DAFT will not reroot trees unless rerooting is explicitly forced.",
                    log,
                    console,
                    "yellow",
                )

            species_tree,sp_tree_rec = self.parse_tree(sp, "Species tree")
            validation_rows.append(["Species tree parsed", "OK", ""])

            species_taxa_list = self.check_duplicate_taxa(species_tree, "Species tree",0)
            species_taxa = set(species_taxa_list)

            validation_rows.append(["Species duplicate taxa", "OK", "No duplicate taxa found"])

            gene_tree_strings = self.read_gene_tree_strings(gene_treefile)

            if not gene_tree_strings:
                raise ValueError("No gene trees were found in " + str(gene_treefile) + ".")

            original_gene_tree_count = len(gene_tree_strings)
            clean_gene_tree_strings = []
            skipped_duplicate_trees = []

            gene_trees = []
            gene_taxa_counts = Counter()
            missing_taxa_warnings = []

            for i, gt_string in enumerate(gene_tree_strings, start=1):
                gt,gt_tree_rec = self.parse_tree(gt_string, "Gene tree " + str(i))
                taxa = self.check_duplicate_taxa(gt, "Gene tree " + str(i),ignore_duplicate)

                if not taxa:
                    skipped_duplicate_trees.append("Gene tree " + str(i))
                    self.print_panel(
                    "DAFT warning",
                    "Gene tree with duplicate detected.\n"
                    "Skipping because --ignore_duplication 1 was used.\n"
                    f"DAFT will use skip input trees{gt_string} ",
                    log,
                    console,
                    "yellow")
                    self.log_message("WARNING: Gene tree " + str(i) + " skipped because it contains duplicate taxa.", log)
                    continue

                missing_taxa = self.check_gene_taxa_against_species(taxa,species_taxa,"Gene tree " + str(i),)

                if missing_taxa:
                    missing_taxa_warnings.append("Gene tree " + str(i) + " is missing species-tree taxa: " + ", ".join(sorted(missing_taxa)))

                for taxon in taxa:
                    gene_taxa_counts[taxon] += 1

                clean_gene_tree_strings.append(gt_string)
                gene_trees.append(gt)

            gene_tree_strings = clean_gene_tree_strings

            if not gene_trees:
                raise ValueError("No usable gene trees remain after duplicate-taxon filtering.")

            if skipped_duplicate_trees:
                validation_rows.append(["Duplicate gene trees", "WARNING", str(len(skipped_duplicate_trees)) + " skipped"])

                clean_file = output_prefix + "_duplicate_removed_gene_trees.csv"
                self.write_rooted_gene_trees(gene_tree_strings, clean_file)
                gene_treefile = clean_file

            else:
                validation_rows.append(["Duplicate gene trees", "OK", "No duplicate taxa found"])

            validation_rows.append(["Gene trees parsed", "OK", str(len(gene_trees)) + " usable of " + str(original_gene_tree_count)])

            if missing_taxa_warnings:
                validation_rows.append(["Missing taxa", "WARNING", str(len(missing_taxa_warnings)) + " gene trees missing taxa"])

                for warning in missing_taxa_warnings[:10]:
                    self.log_message("WARNING: " + warning, log)

                if len(missing_taxa_warnings) > 10:
                    self.log_message("WARNING: more missing-taxa warnings exist but are not printed here.", log)

            else:
                validation_rows.append(["Missing taxa", "OK", "No missing species-tree taxa"])

            rooted_file = None

            if rooting and forced:
                self.print_panel(
                    "Forced rooting requested",
                    "Outgroup: " + ", ".join(outgroup_taxa),
                    log,
                    console,
                    "yellow",
                )

                species_tree = self.root_tree_on_outgroup(species_tree,outgroup_taxa,"Species tree")

                rooted_strings = []
                rooted_gene_trees = []

                for i, gt in enumerate(gene_trees, start=1):
                    gt = self.root_tree_on_outgroup(gt,outgroup_taxa,"Gene tree " + str(i))

                    rooted_strings.append(gt.write(format=1))
                    rooted_gene_trees.append(gt)

                rooted_file = output_prefix + "_rooted_gene_trees.csv"
                self.write_rooted_gene_trees(rooted_strings, rooted_file)

                gene_tree_strings = rooted_strings
                gene_trees = rooted_gene_trees
                gene_treefile = rooted_file

                validation_rows.append(["Forced rooting", "OK", rooted_file])

            self.check_binary_rooted(species_tree, "Species tree")
            validation_rows.append(["Species tree binary/rooted", "OK", ""])

            species_root_split = self.root_split(species_tree)

            inconsistent = []
            unchecked_rooting = []

            for i, gt in enumerate(gene_trees, start=1):
                self.check_binary_rooted(gt, "Gene tree " + str(i))

                try:
                    checked = self.check_root_consistency(
                        species_tree,
                        gt,
                        "Gene tree " + str(i),
                    )

                    if checked is False:
                        unchecked_rooting.append("Gene tree " + str(i))

                except ValueError as err:
                    inconsistent.append(str(err))

            if unchecked_rooting:
                validation_rows.append([
                    "Root consistency skipped",
                    "WARNING",
                    str(len(unchecked_rooting)) + " gene trees could not be checked",
                ])

                self.log_message(
                    "WARNING: Some gene trees contain taxa from only one side of the species-tree root, "
                    "so root consistency could not be checked.",
                    log,
                )

            if inconsistent and not allow_inconsistent_rooting:
                raise ValueError(
                    "Some gene trees are rooted inconsistently with the species tree. "
                    "Use --rooting OUTGROUP --forced 1 to reroot, or "
                    "--allow_inconsistent_rooting 1 to continue without rerooting (for example, if the gene trees are already rooted as intended)."
                    "First inconsistent tree: " + inconsistent[0]
                )

            if inconsistent and allow_inconsistent_rooting:
                validation_rows.append(["Root consistency", "WARNING", str(len(inconsistent)) + " inconsistent gene trees"])

                self.print_panel(
                    "DAFT warning",
                    "Inconsistent rooting detected.\n"
                    "Continuing because --allow_inconsistent_rooting 1 was used.\n"
                    "DAFT will use the input trees unchanged.",
                    log,
                    console,
                    "yellow",
                )

            else:
                validation_rows.append(["Root consistency", "OK", "All checkable gene trees are consistent"])

            self.print_tree_summary(
                species_taxa,
                gene_tree_strings,
                species_root_split,
                log,
                console,
            )

            rows = self.build_taxon_report(
                species_taxa=species_taxa,
                gene_tree_taxa_counts=gene_taxa_counts,
            )

            self.check_internal_name_collisions(rows)
            validation_rows.append(["DAFT internal names", "OK", "No collisions found"])

            taxon_report_path = self.write_taxon_report(rows, output_prefix)
            validation_rows.append(["Taxon report", "OK", taxon_report_path])

            input_hash = self.compute_input_hash(
                species_tree.write(format=1),
                gene_tree_strings,
            )

            validation_rows.append(["djiNNI cache hash", "OK", input_hash])

            self.print_validation_summary(validation_rows, log, console)

            self.print_panel(
                "DAFT validation complete",
                "Input validation complete.\n\n" +
                "Taxon report: " + str(taxon_report_path) + "\n" +
                "djiNNI cache input hash: " + str(input_hash),
                log,
                console,
                "green",
            )

            log.close()

            self.species_tree_newick = species_tree.write(format=1)
            self.gene_treefile = gene_treefile
            self.gene_tree_strings = gene_tree_strings
            self.input_hash = input_hash
            self.taxon_report_path = taxon_report_path
            self.rooted_gene_treefile = rooted_file

            return self

        except Exception as err:
            message = (
                "Error: " + str(err) + "\n\n" +
                "Please check that:\n" +
                "- all trees are valid Newick strings\n" +
                "- all trees end with ';'\n" +
                "- all trees are binary rooted trees\n" +
                "- taxa are not duplicated within any tree\n"
                "- gene-tree taxa are present in the species tree\n" +
                "- gene trees are rooted consistently with the species tree\n" +
                "- use --rooting OUTGROUP --forced 1 only if you want DAFT to reroot\n" +
                "- use --allow_inconsistent_rooting 1 only if you want to continue unchanged\n"
                "- use --ignore_duplication 1 only if you want to ignore gene tree with duplicates\n"
            )

            if log is not None:
                self.print_panel(
                    "DAFT input validation failed",
                    message,
                    log,
                    console,
                    "red",
                )
                log.close()

            print("DAFT input validation failed.", file=sys.stderr)
            print(message, file=sys.stderr)

            sys.exit(1)