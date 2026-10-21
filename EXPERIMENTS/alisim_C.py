import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'DATA_SET/list_A.csv')
    parser.add_argument('--output', type=Path, default=ROOT / 'DATA_SET/A_alisim_500bp_indel004')
    parser.add_argument('--iqtree', type=Path, default=ROOT / 'iqtree3')
    parser.add_argument('--length', type=int, default=300)
    parser.add_argument('--indel', type=float, default=0.04)
    parser.add_argument('--model', default='JC')
    parser.add_argument('--branch-length', type=float, default=0.01,
                        help='Length assigned to missing branches (substitutions/site).')
    parser.add_argument('--threads', type=int, default=1, help='Threads for each reconstruction.')
    parser.add_argument('--seed', type=int, default=12345)
    args = parser.parse_args()
    if (args.length < 1 or args.indel < 0 or args.branch_length <= 0
            or args.threads < 1 or args.seed < 1):
        parser.error('Lengths, threads and seed must be positive; indel must be nonnegative.')
    with args.input.open(newline='') as handle:
        reader = csv.DictReader(handle)
        if 'gt' not in (reader.fieldnames or []):
            parser.error('Input CSV must contain a gt column.')
        trees = [row['gt'].strip() for row in reader]
    if not trees or any(not tree.endswith(';') for tree in trees):
        parser.error('Input must contain nonempty semicolon-terminated Newick trees.')
    output = args.output.resolve()
   
    output.mkdir(parents=True, exist_ok=False)
    for folder in ('input_trees', 'alignments', 'trees', 'logs'):
        (output / folder).mkdir()
    config = {k: str(v.resolve()) if isinstance(v, Path) else v for k, v in vars(args).items()}
    config['tree_count'] = len(trees)
    (output / 'settings.json').write_text(json.dumps(config, indent=2) + '\n')
    executable = str(args.iqtree.resolve())

    def run(command, logfile):
        with logfile.open('w') as log:
            log.write('Command: ' + repr(command) + '\n')
            log.flush()
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)

    with (output / 'reconstructed_trees.nwk').open('w') as combined, \
            (output / 'manifest.csv').open('w', newline='') as manifest:
        writer = csv.writer(manifest)
        writer.writerow(['row', 'input_tree', 'alignment', 'reconstructed_tree', 'seed'])
        for index, tree in enumerate(trees, 1):
            name = f'gene_{index:04d}'
            # For this CSV's unquoted Newick labels, add missing edge lengths
            # before sibling/parent delimiters. The root needs no edge length.
            tree = re.sub(r'([^(),;]+|\))(?=\s*[,\)])',
                          lambda m: (m.group(1) if ':' in m.group(1) else
                                     m.group(1).rstrip() + ':' + str(args.branch_length)), tree)
            source = output / 'input_trees' / f'{name}.nwk'
            source.write_text(tree + '\n')
            aln_prefix = output / 'alignments' / name
            tree_prefix = output / 'trees' / name
            seed = args.seed + index - 1
            run([executable, '--alisim', str(aln_prefix), '-t', str(source),
                 '-m', args.model, '--length', str(args.length), '--indel',
                 f'{args.indel},{args.indel}', '--seed', str(seed), '-nt', '1'],
                output / 'logs' / f'{name}.simulation.log')
            alignment = aln_prefix.with_suffix('.phy')
            run([executable, '-s', str(alignment), '-m', args.model,
                 '-nt', str(args.threads), '--seed', str(seed), '--prefix', str(tree_prefix)],
                output / 'logs' / f'{name}.inference.log')
            inferred = tree_prefix.with_suffix('.treefile')
            combined.write(inferred.read_text().strip() + '\n')
            combined.flush()
            writer.writerow([index, source.relative_to(output), alignment.relative_to(output),
                             inferred.relative_to(output), seed])
            manifest.flush()
            print(f'[{index}/{len(trees)}] {inferred}', flush=True)
    print(f'Finished: {output}')


if __name__ == '__main__':
    try:
        main()
    except (OSError, subprocess.CalledProcessError) as exc:
        sys.exit(f'Error: {exc}. See the output logs; use a fresh --output directory to rerun.')
