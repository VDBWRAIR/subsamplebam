#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
import numpy as np


from Bio import SeqIO

def parse_args(wrapper=False):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'bamfile'
    )

    if not wrapper: 
        parser.add_argument(
            'refseq',
            default=None,
            help='name of the reference sequenceto do subsampleing on '
        )
        ''' reference lenght isn't necessary to use samtools view '''
#        parser.add_argument(
#            'reflength',
#            default=None,
#            type=int,
#            help='Length of the reference sequence'
#        )

    parser.add_argument(
        '--subsample',
        type=int,
        default=1000,
        help='What depth to try and subsample to'
    )

    parser.add_argument(
            '-A', '--count-orphans',
            action='store_true',
            help='Allow orphan/unpaired reads.'
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


def reads_by_pos(reads, max_pos):
       #def is_at_pos(pos): return read['pos'] == pos
    for pos in xrange(max_pos):
        yield [read for read in reads if read['pos'] == pos]

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
        sout, _ = p.communicate()
        findme = None
        if isinstance(sout, bytes):
            findme = b'Version:'
        else:
            findme = 'Version:'
        if sout.find(findme) >= 0:
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
