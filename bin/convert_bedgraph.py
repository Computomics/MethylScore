#!/usr/bin/env python

from cytosine import Cytosine, get_input, reference_gen

def bedgraph2format(input, seq, chromosome, format, output, sample=None):
    """Converts bedgraph format to either MethylScore, DSS, HOME, methylKit or methylpy format

    Args:
        input: str
            output file path
        format: str
            desired output format
        output: str
            output file path
        sample: str
            sampleID
    """
    assert format in ['DSS','HOME','methylpy', 'methylKit', 'MethylScore'], 'Only MethylScore, DSS, HOME, methylKit and methylpy formats are supported at this time'
    sep='\t'
    prev = None

    with get_input(input) as bedgraph:
        with open(output, 'w') as outfile:
            for l in bedgraph:
                fields = l.split(sep)

                if format == 'MethylScore':
                    # the following code allows to treat replicates as technical ones by summing up meth. information for each position if all information is present in a single (sorted) file
                    # In the default workflow however, given that the samplesheet contains unique sampleIDs, it merely prepends the sampleID to the file
                    while prev == int(fields[2]) or prev is None:
                        break
                    else:
                        c = Cytosine(chromosome, prev, meth, unmeth, reference=seq)
                        outfile.write(sep.join(str(f) for f in [sample, c.chrom, c.pos-1, c.pos, f'{c.get_rate():.2f}', c.meth, c.unmeth]) + '\n')
                        del meth
                        del unmeth
                    try:
                        meth += int(fields[4])
                        unmeth += int(fields[5])
                    except NameError:
                        meth = int(fields[4])
                        unmeth = int(fields[5])
                    finally:
                        prev = int(fields[2])
                elif format == 'methylKit':
                    c = Cytosine(str(fields[0]), int(fields[2]), int(fields[4]), int(fields[5]), reference=seq)
                    coverage = c.get_coverage()
                    freqC = f'{c.meth/coverage*100:.2f}'
                    freqT = f'{c.unmeth/coverage*100:.2f}'
                    outfile.write(sep.join(str(f) for f in [c.chrom, c.chrom, c.pos, ('F' if c.strand == 'C' else 'R'), coverage, freqC, freqT]) + '\n')
                elif format == 'DSS':
                    c = Cytosine(str(fields[0]), int(fields[2]), int(fields[4]), int(fields[5]), reference=seq)
                    outfile.write(sep.join(str(f) for f in [c.chrom, c.pos, c.get_coverage(), c.meth]) + '\n')
                elif format == 'HOME':
                    c = Cytosine(str(fields[0]), int(fields[2]), int(fields[4]), int(fields[5]), reference=seq)
                    outfile.write(sep.join(str(f) for f in [c.chrom, c.pos, ('+' if c.strand == 'C' else '-'), c.context, c.meth, c.get_coverage()]) + '\n')
                elif format == 'methylpy':
                    c = Cytosine(str(fields[0]), int(fields[2]), int(fields[4]), int(fields[5]), reference=seq)
                    outfile.write(sep.join(str(f) for f in [c.chrom, c.pos, ('+' if c.strand == 'C' else '-'), c.triplet, c.meth, c.get_coverage(), 'NA']) + '\n')

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description = 'Converts per-cytosine bedGraph to DSS, HOME, methylKit or methylpy input format')

    parser.add_argument('--input', type=str, help='per-cytosine information', dest='consensus')
    parser.add_argument('--format', type=str, help='output file format', dest='format', required=True)
    parser.add_argument('--output', type=str, help='output file name', dest='output', required=True)
    parser.add_argument('--fasta', type=str, help='reference genome', dest='fasta')
    parser.add_argument('--sample', type=str, help='name of the sample', dest='sample')
    parser.add_argument('--verbose', type=bool, help='print skipped sites to stdout', dest='verbose', default=False)

    args = parser.parse_args()

    ref = reference_gen(args.fasta)
    
    for chrom, seq in ref:
        bedgraph2format(args.consensus, seq, chrom, args.format, args.output, args.sample)

