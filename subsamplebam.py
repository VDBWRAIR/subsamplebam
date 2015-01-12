#!/usr/bin/env python

# Handle SIGPIPE signals correctly
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

# Read given file with read names and hash them for quick lookup
# Read stdin and filter based on samtools view input and spit to stdout

import sys
import argparse
import random
import subprocess
import shlex
import contextlib
import os.path

from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'bamfile'
    )

    parser.add_argument(
        'reffile'
    )

    parser.add_argument(
        'regionstr',
        default=None,
        help='Region string to do subsampleing on ref with'
    )

    parser.add_argument(
        '--subsample',
        type=int,
        default=1000,
        help='What depth to try and subsample to'
    )

    parser.add_argument(
        '--output',
        default='-',
        help='Where to store the resulting filtered bam[Default: stdout]'
    )

    return parser.parse_args()

class MinimumDepthNotAvailable(Exception): pass

def randomly_select_read_from_samview(bamfile, regionstr, n):
    '''
    Given a samview output select n random read names
    In place operation on selection

    :param str bamfile: bamfile path
    :param str regionstr: regionstring to operate on
    :param int n: Random subselection depth, aka, how many random rows to select from samview
    '''
    newselection = None
    with samview(bamfile, regionstr) as sview:
        sview = list(sview)
        try:
            newselection = set(random.sample(sview, n))
        except ValueError as e:
            raise MinimumDepthNotAvailable('Depth for {0} is only {1}'.format(regionstr, len(sview)))
            #raise MinimumDepthNotAvailable('{0} for {1}\n'.format(str(e), regionstr))
    return newselection

def reference_info(reffile):
    '''
    Return information about references in file
    {'ref1name': ref1len, ...}

    >>> i = reference_info('Den1__Cambodia__2009.fasta')
    >>> print i['Den1/GU131895_1/Cambodia/2009/Den1_1']
    10474
    '''
    refinfo = {}
    for rec in SeqIO.parse(reffile, 'fasta'):
        refinfo[rec.id] = len(str(rec.seq))
    return refinfo

def parallel_randomly_select_read_from_samview(args):
    try:
        return randomly_select_read_from_samview(*args)
    except MinimumDepthNotAvailable as e:
        sys.stderr.write(e + '\n')
        return set()

def subselect_from_bam(bamfile, subselectdepth, reffile, regionstr):
    '''
    Iterate over every base position in regionstr and subselect randomly from the reads
    that fall under it.
    '''
    import multiprocessing

    # Get names and lengths of reference
    refinfo = reference_info(reffile)
    # Parse the regionstring
    refname, startstop = regionstr.split(':')
    start, stop = startstop.split('-')
    uniquereads = set()
    rstrings = []
    for i in range(int(start), int(stop)+1):
        rstring = '{0}:{1}-{1}'.format(refname, i, i)
        #selected = randomly_select_read_from_samview(bamfile, rstring, subselectdepth)
        #uniquereads.update(selected)
        rstrings.append((bamfile, rstring, subselectdepth))

    pool = multiprocessing.Pool()
    ureads = pool.map(parallel_randomly_select_read_from_samview, rstrings)
    for uread in ureads:
        uniquereads.update(uread)

    return uniquereads

def make_subselected_bam(bamfile, uniquereads, outbam):
    '''
    Make outbam with uniquereads
    '''
    name, ext = os.path.splitext(outbam)

    # Write header of bamfile
    cmd = 'samtools view -H {0}'.format(bamfile)
    sys.stderr.write(cmd + '\n')
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        shell=True
    )

    for line in p.stdout:
        sys.stdout.write(line)
    for line in uniquereads:
        sys.stdout.write(line)

@contextlib.contextmanager
def samview(bamfile, regionstr):
    '''
    Just return iterator over samtools view bamfile regionstr
    '''
    cmd = 'samtools view {0} {1}'.format(bamfile, regionstr)
    sys.stderr.write(cmd + '\n')
    p = subprocess.Popen(
        shlex.split(cmd),
        stdout=subprocess.PIPE
    )
    yield p.stdout

def samtools_is_available():
    try:
        p = subprocess.Popen(
            ['samtools'],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )
        sout,_ = p.communicate()
        if 'Version' in sout:
            return True
        return False
    except OSError as e:
        return False
 
def main():
    if not samtools_is_available():
        sys.stderr.write('Samtools is not installed or executable\n')
        sys.exit(1)
    args = parse_args()
    uniquereads = subselect_from_bam(args.bamfile, args.subsample, args.reffile, args.regionstr)
    make_subselected_bam(args.bamfile, uniquereads, args.output) 

if __name__ == '__main__':
    main()
