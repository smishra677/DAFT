
import argparse
import csv
import math
from pathlib import Path
import re

import dendropy


def natural_key(path):
    return [int(part) if part.isdigit() else part
            for part in re.split(r'(\d+)', path.as_posix())]


def leaf_labels(node):
    return frozenset(leaf.taxon.label for leaf in node.leaf_nodes())


def root_from_reference(tree, newick):
    references = dendropy.TreeList.get(
        data=newick, schema='newick', rooting='force-rooted',
        preserve_underscores=True,
    )
    if len(references) != 1:
        raise ValueError('Expected exactly one reference tree')
    reference = references[0]
    children = reference.seed_node.child_nodes()
    if len(children) != 2:
        raise ValueError('Reference must have exactly two root children')
    if leaf_labels(reference.seed_node) != leaf_labels(tree.seed_node):
        raise ValueError('Reference and target taxa differ')
    groups = {leaf_labels(child) for child in children}
    try:
        for node in tree.preorder_node_iter():
            if node.parent_node is not None and leaf_labels(node) in groups:
                length = node.edge.length
                if length is not None and (not math.isfinite(length) or length < 0):
                    raise ValueError('Root edge length must be finite and nonnegative')
                half = length / 2 if length is not None else None
                tree.reroot_at_edge(node.edge, length1=half, length2=half)
                return
        raise ValueError('Reference root split is absent from target tree')
    except:
        tree.resolve_polytomies(limit=2)
        return


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path, help='Folder containing IQ-TREE .treefile outputs')
    parser.add_argument('--output', '-o', type=Path, required=True, help='Output CSV path')
    parser.add_argument('--recursive', action='store_true', help='Include subfolders')
    parser.add_argument('--overwrite', action='store_true', help='Replace an existing output CSV')
    parser.add_argument('--reference-csv', type=Path, required=True,
                        help='Reference trees in a gt-column CSV, paired by row with naturally sorted files')
    args = parser.parse_args()
    if not args.folder.is_dir():
        parser.error(f'Not a directory: {args.folder}')
    if args.output.exists() and not args.overwrite:
        parser.error(f'Output exists: {args.output}. Use --overwrite to replace it.')
    paths = args.folder.rglob('*.treefile') if args.recursive else args.folder.glob('*.treefile')
    files = sorted((p for p in paths if p.is_file()), key=natural_key)
    if not files:
        parser.error(f'No .treefile files found in {args.folder}')
    if args.output.resolve() in {p.resolve() for p in files}:
        parser.error('Output cannot replace an input tree file.')

    if args.output.resolve() == args.reference_csv.resolve():
        parser.error('Output cannot replace the reference CSV.')
    try:
        with args.reference_csv.open(newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or 'gt' not in reader.fieldnames:
                raise ValueError('Reference CSV must have a gt column')
            references = []
            for index, row in enumerate(reader, start=1):
                newick = (row.get('gt') or '').strip()
                if not newick:
                    raise ValueError(f'Reference row {index} has an empty gt value')
                references.append(newick)
        #if len(references) != len(files):
        #raise ValueError(f'{len(references)} reference rows but {len(files)} tree files')
    except (OSError, ValueError, csv.Error) as exc:
        parser.error(f'{args.reference_csv}: {exc}')

    # Validate every input before opening the output file.
    rows = []
    for index, path in enumerate(files):
        try:
            trees = dendropy.TreeList.get(
                path=str(path),
                schema="newick",
                preserve_underscores=True
            )
            if len(trees) != 1:
                raise ValueError(f'Expected one tree, found {len(trees)}')
            tree = trees[0]
            root_from_reference(tree, references[index])
            rows.append([tree.as_string(schema='newick', suppress_rooting=True).strip()])
        except Exception as exc:
            context = f' (reference row {index + 1})'
            parser.error(f'{path}{context}: {exc}')

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w' if args.overwrite else 'x', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['gt'])
        writer.writerows(rows)
    print(f'Saved {len(rows)} rerooted trees to {args.output}.')


if __name__ == '__main__':
    main()
