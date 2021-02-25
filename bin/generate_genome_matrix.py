#!/usr/bin/env python

import re
import gzip

revcomp = str.maketrans("ACGT","TGCA")

CG = re.compile('^CG')
CHG = re.compile('^C.G')
CHH = re.compile('^C')

def open_text(path):
    """Reads a textfile"""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    else:
        return open(path)

def get_samples(path):
    """Reads a tsv file with sample names in first column and returns a dict of samples"""
    samples = {}
    with open(path) as s:
        for line in s:
           id, path = line.split()
           samples[id] = path
    return samples

def reference_gen(path):
    """Returns an iterator over sequences in a fasta file"""
    def _get_chrom(header):
        return header.split()[0][1:]
    seq = ''
    with open(path) as fa:
        # first line is always a fasta header
        c = _get_chrom(fa.readline())
        for line in fa:
           if line.startswith('>'):              
               yield c, seq.upper()
               c = _get_chrom(line)
               seq = ''
           else:
                seq += line.rstrip('\n')
    # always yield when the file has been read through
    yield c, seq.upper()

class Cytosine( ):
    """
    A class used to represent a single Cytosine position.

    Attributes
    ----------
    pos : int
        the position on the reference
    meth : int
        number of reads supporting a methylated state
    unmeth : str
        number of reads not supporting a methylated state
    context : str
        the sequence context
    strand : str
        the strand position, either Watson (C) or Crick (G)

    Methods
    -------
    get_strand(reference)
        Gets the strand nucleotide at pos on the reference sequence

    get_context(strand)
        Gets the nucleotide context at pos.
        Raises ValueError in case of ambiguous contexts or polymorphisms

    get_rate()
        Gets the methylation rate based on methylated and unmethylated reads

    get_valuestring()
        Gets a string of the form 40/rate/meth/unmeth for the position

    """
    def __init__(self, pos, meth, unmeth, context=None, strand=None, reference=None):

        self.pos = pos
        self.meth = meth
        self.unmeth = unmeth
        if reference is not None:
            self.strand = self.get_strand(reference)
            self.context = self.get_context(reference, self.strand)
        else:
            self.strand = strand
            self.context = context

    def get_strand(self, reference):
        return reference[self.pos-1]

    def get_context(self, reference, strand):
        def _reverse_complement(triplet):
            return triplet.translate(revcomp)[::-1]

        if strand == 'C':
            triplet = reference[self.pos-1:self.pos+2]
        elif strand == 'G':
            triplet = _reverse_complement(reference[self.pos-3:self.pos])
        else:
            raise ValueError(f'position {self.pos} is not a C or G in the reference')

        if 'N' in triplet:
            raise ValueError(f'triplet {triplet} at pos {self.pos} is ambiguous, skipping')

        if CG.match(triplet):
            context = "CG"
        elif CHG.match(triplet):
            context = "CHG"
        elif CHH.match(triplet):
            context = "CHH"
        return context

    def get_rate(self):
        return float(self.meth / (self.meth + self.unmeth))*100

    def get_valuestring(self):
        return f'40/{self.get_rate():.1f}/{self.meth}/{self.unmeth}'


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

    with open_text(path) as infile:
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
                    c = Cytosine(pos, meth, unmeth, reference=seq)
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