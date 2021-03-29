#!/usr/bin/env python

import re
import gzip

def get_input(path):
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

class Cytosine():
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

    get_triplet(reference)
        Gets the triplet starting at pos.
        Raises ValueError in case of polymorphisms

    get_context()
        Gets the nucleotide context at pos.
        Raises ValueError in case of ambiguous contexts

    get_rate()
        Gets the methylation rate based on methylated and unmethylated reads

    get_coverage()
        Gets the coverage as sum of methylated and unmethylated reads

    get_valuestring()
        Gets a string of the form 40/rate/meth/unmeth for the position

    """
    revcomp = str.maketrans("ACGT","TGCA")
    CG = re.compile('^CG')
    CHG = re.compile('^C.G')
    CHH = re.compile('^C')

    def __init__(self, chrom, pos, meth, unmeth, reference, strand=None, triplet=None, context=None):

        self.chrom = chrom
        self.pos = pos
        self.meth = meth
        self.unmeth = unmeth
        self.strand = self.get_strand(reference)
        self.triplet = self.get_triplet(reference)
        self.context = self.get_context()

    def get_strand(self, reference):
        return reference[self.pos-1]

    def get_triplet(self, reference):
        def _reverse_complement(triplet):
            return triplet.translate(Cytosine.revcomp)[::-1]
        if self.strand == 'C':
            return reference[self.pos-1:self.pos+2]
        elif self.strand == 'G':
            return _reverse_complement(reference[self.pos-3:self.pos])
        else:
            raise ValueError(f'position {self.pos} is not a C or G in the reference')

    def get_context(self):
        if 'N' in self.triplet:
            raise ValueError(f'triplet {self.triplet} at pos {self.pos} is ambiguous')

        if Cytosine.CG.match(self.triplet):
            return "CG"
        elif Cytosine.CHG.match(self.triplet):
            return "CHG"
        elif Cytosine.CHH.match(self.triplet):
            return "CHH"

    def get_rate(self):
        return float(self.meth / (self.meth + self.unmeth))*100

    def get_coverage(self):
        return int(self.meth + self.unmeth)

    def get_valuestring(self):
        return f'40/{self.get_rate():.1f}/{self.meth}/{self.unmeth}'