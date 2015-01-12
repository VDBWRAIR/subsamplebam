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
import os.path
import multiprocessing

from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'bamfile'
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

    return parser.parse_args()

class MinimumDepthNotAvailable(Exception): pass

def randomly_select_reads_from_samview(bamfile, regionstr, n):
    '''
    Given a samview output select n random read names
    In place operation on selection

    :param str bamfile: bamfile path
    :param str regionstr: regionstring to operate on
    :param int n: Random subselection depth, aka, how many random rows to select from samview
    '''
    newselection = None
    sview = samview(bamfile, regionstr)
    sview = list(sview)
    try:
        newselection = set(random.sample(sview, n))
    except ValueError as e:
        raise MinimumDepthNotAvailable('Depth for {0} is only {1}'.format(regionstr, len(sview)))
    return newselection

def reference_info(reffile):
    '''
    Hash reference id's with their lengths
    '''
    refinfo = {}
    for rec in SeqIO.parse(reffile, 'fasta'):
        refinfo[rec.id] = len(str(rec.seq))
    return refinfo

def parallel_randomly_select_reads_from_samview(args):
    try:
        return randomly_select_reads_from_samview(*args)
    except MinimumDepthNotAvailable as e:
        sys.stderr.write(str(e) + '\n')
        return set()

def subselect_from_bam(bamfile, subselectdepth, regionstr):
    '''
    Iterate over every base position in regionstr and subselect randomly from the reads
    that fall under it.
    '''
    # Parse the regionstring
    refname, startstop = regionstr.split(':')
    start, stop = startstop.split('-')
    uniquereads = set()
    rstrings = []
    for i in range(int(start), int(stop)+1):
        rstring = '{0}:{1}-{1}'.format(refname, i, i)
        rstrings.append((bamfile, rstring, subselectdepth))

    pool = multiprocessing.Pool()
    ureads = pool.map(parallel_randomly_select_reads_from_samview, rstrings)
    #ureads = map(parallel_randomly_select_reads_from_samview, rstrings)
    for uread in ureads:
        uniquereads.update(uread)

    return uniquereads

def make_subselected_bam(bamfile, uniquereads):
    '''
    Make outbam with uniquereads
    '''
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
    # you have to read stdout in order for return to fill in
    firstline = p.stdout.readline()
    if p.poll() not in (0,None) and firstline == '':
        raise ValueError('samtools did not return correctly')
    yield firstline
    for line in p.stdout:
        yield line

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
    uniquereads = subselect_from_bam(args.bamfile, args.subsample, args.regionstr)
    make_subselected_bam(args.bamfile, uniquereads) 
