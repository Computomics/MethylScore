#!/usr/bin/env python

from cytosine import Cytosine, get_input, get_samples, reference_gen 

def convert_to_genome_matrix(path, seq, chromosome, samples, format='bedgraph', verbose=False):
    """Converts per-cytosine methylation information to genome matrix format required by MethylScore

    Note that the input file needs to be preprocessed to contain the sample name in the first column!

    Args:
        path: str
            Position sorted input file
        seq: str
            Fasta sequence string
        chromosome: str
            Chromosome ID
        samples: dict
            Sample list
        format: string
            Currently supports 'bedgraph' and 'methylpy'
        verbose: bool
            Prints skipped positions to stdout
    """
    sep='\t'

    skipped_sites = set()
    line = []
    coord = []
    values = {}
    prev = None

    with get_input(path) as infile:
        with open(f'{chromosome}.genome_matrix.tsv', 'w') as outfile:
            outfile.write(f'#chr\tpos\tclass\tstrand\t{sep.join(samples)}\n')

            for l in infile:

                fields = l.split(sep)

                if format == 'bedgraph':
                    sample, chrom, pos, meth, unmeth = str(fields[0]), str(fields[1]), int(fields[3]), int(fields[5]), int(fields[6])
                elif format == "methylpy":
                    sample, chrom, pos, meth, unmeth = str(fields[0]), str(fields[1]), int(fields[2]), int(fields[5]), (int(fields[6]) - int(fields[5]))
                else:
                    raise NotImplementedError('Only bedGraph and methylpy format are supported at this time')

                if (chrom != chromosome) or (pos in skipped_sites):
                    continue

                try:
                    c = Cytosine(chrom, pos, meth, unmeth, reference=seq)
                except ValueError as e:
                    skipped_sites.add(pos)
                    if verbose:
                        print(e)
                    continue

                while pos == prev or prev is None:
                    break
                else:
                    for strain in samples:
                        line.append(values.get(strain, '0'))

                    outfile.write(f'{sep.join(coord + line)}\n')
                    line.clear()
                    coord.clear()
                    values.clear()

                prev = pos

                if not coord:
                    coord.extend([chrom, str(pos), c.context, c.strand])

                values[sample] = c.get_valuestring()

    print(f'Skipped {len(skipped_sites)} site(s) due to ambiguous contexts or polymorphisms')

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description = 'Converts per-cytosine methylation info to genome matrix format')

    parser.add_argument('--allC', type=str, help='per-cytosine information', dest='consensus')
    parser.add_argument('--format', type=str, help='input file format', dest='format', default='bedgraph')
    parser.add_argument('--fasta', type=str, help='reference genome', dest='fasta')
    parser.add_argument('--verbose', type=bool, help='print skipped sites to stdout', dest='verbose', default=False)
    parser.add_argument('--samples', type=str, help='sample file containing sample ids and paths', dest='samples', required=True)

    args = parser.parse_args()

    # read fasta file
    ref = reference_gen(args.fasta)

    # read sample file
    samples = get_samples(args.samples).keys()

    # generate genome matrix for each chromosome
    for chrom, seq in ref:
        convert_to_genome_matrix(args.consensus, seq, chrom, samples, args.format, args.verbose)