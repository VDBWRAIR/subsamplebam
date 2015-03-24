from __future__ import print_function
import sys
#import pandas as pd 
import numpy as np
import subsamplebam
import subprocess as sp
import re


''' bugs '''

def get_raw_reads(bamfile, refseq='', idx=None):
 return list(subsamplebam.samview(bamfile,  refseq))
 #pass

def get_next_reads(reads, pos):
   return [read for read in reads if read['pos'] == pos]

'''
Only add a read if the current depth is less than the target depth.
Unfortunately this doesn't work because you can be under the target depth at a position that has no alignments starting there.
''' 

def parse_alignment(raw_view):
   fields = raw_view.split()
   qname, pos, seq = fields[0], int(fields[3]), fields[9]
   return Alignment(raw_view, qname, pos, seq)

def add_lists():
   for j, depth in enumerate(current_depths[i:]): next_depths_array[j] += depth
   current_depths = next_depths

'''
0. Store sequences at indices with (sequence, QName) so that they can be uniquely identified.
1. set overlap_target to under_index.
2. get all sequences which could overlap under_index, and which have not already been selected. (how get these, by using samtools view again? or store in memory? only need a UNIQUE identifier and the sequence length)
3. choose sequences a. randomly b. minimize overflow (the sequence closest to under_index) this may result in a bias towards reads beginning at a given index.  c. minimize overflow probabilistically. d. the sequence which has the longest length after under_index. This will deal with the situation where you picked a short sequence which didn't cover a coming index, before you knew you wouldn't need it.
'''
class DepthMatrix(object):
    def __init__(self, ref_length, min_depth):
        self.seq_matrix = [[]]
        self.depth_array = np.zeros(ref_length)
        self.min_depth = min_depth

    def get_candidate_sequences(self, under_index):
        #have this only look backwards until reach  under-covered index 
        sub_matrix = self.seq_matrix[:under_index+1]
        #flatten the list
        prev_seqs = [seq for row in sub_matrix for seq in row]
        return filter(lambda seq: seq.overlap  >= under_index and not seq.picked, prev_seqs)
    
    def yield_greatest_overlaps(self, under_index, num_needed):
        candidate_sequences = self.get_candidate_sequences(under_index) 
        matches = 0
        #import ipdb; ipdb.set_trace()
        while matches < num_needed and candidate_sequences: 
            # could instead keep sequences in order sorted by overlap 
            farthest_overlap_seq = max(candidate_sequences, key=lambda seq: seq.overlap)
            candidate_sequences.remove(farthest_overlap_seq)
            farthest_overlap_seq.pick()
            matches += 1
            yield farthest_overlap_seq

    def backtrack(self, pos, needed_depth):
        next_sequences = list(self.yield_greatest_overlaps(pos, needed_depth))
        if not next_sequences:  # no matches found
            return False
        ''' This doesn't work because it returns a flat list of depths starting at index 0, where in reality these reads
        Can start at different positions behind the current position.'''
        min_pos = min(seq.pos for seq in next_sequences)
        next_depths = self.get_depths(next_sequences, min_pos)
        #self.depth_array[pos:pos+len(next_depths)] += next_depths[:len(self.depth_array)-pos]
        self.depth_array += next_depths
        #sys.stderr.write("depth {0} at pos {1}\n".format(str(len(next_sequences)), str(pos)))
        if  len(next_sequences) > needed_depth:
            raise Exception("Too many sequences retured")
        return len(next_sequences) == needed_depth

    ''' 
    Needlessly creates a depth_array of equal size to self.depth_array 
    '''
    def get_depths(self, reads, min_pos): 
       depths = np.zeros(len(self.depth_array))
       for seq in reads: 
           ''' Even if a sequence would overlap past the length of depth-array, we don't include that in depth-array. 
           the numpy arrays must be the same size in order to add them.
           '''
           i, j = seq.pos, min(len(depths), seq.overlap +1)
           depths[i:j] += np.array([1 for overlap in xrange(i, j)])
       return depths

    ''' for debugging'''
    def print_matrix(self):
        for row in self.seq_matrix:
            sys.stderr.write( [seq.pos for seq in row])

    def make_seq_matrix(self, bamfile, refseq): 
        all_alignments = [parse_alignment(string) for string in get_raw_reads(bamfile, refseq)]
        #map(parse_alignment, get_raw_reads(bamfile, refseq))
        max_pos = max([seq.pos for seq in all_alignments])
        #initialize sequence matrix as a 2D array
        # weird thing about pos 0 
        #self.seq_matrix = [ filter( lambda seq: seq.pos == i, all_alignments) for i in xrange(1, max_pos)]
        self.seq_matrix = [ [seq for seq in all_alignments if seq.pos == i] for i in xrange(max_pos)]
                #filter( lambda seq: seq.pos == i, all_alignments) for i in xrange(1, max_pos)]
        return

# need to handle case where fasta file does not span entire reference genome;
# depth array should be limited in that case to the largest overlap
    def main(self):
        #import ipdb; ipdb.set_trace()
        for pos, depth in enumerate(self.depth_array):
            if depth < self.min_depth:
                needed_depth = self.min_depth - depth
                depth_met = self.backtrack(pos, needed_depth)
                #sys.stderr.write("depth was {0}met.\n".format("" if depth_met else "NOT "))

class Alignment(object):

    def __init__(self, string, qname, pos, seq):
        self.qname, self.pos, self.seq = qname, pos, seq
        self.picked = False
        self.string = string

    @property
    def seq_length(self):
        return len(self.seq)

    @property
    def overlap(self):
        return self.pos + self.seq_length

    def pick(self):
        self.picked = True
    def __str__(self):
        return "pos: {0}, overlap {1} \n seq {2} \n {3}".format(self.pos, self.overlap, self.seq, self.string)

def get_num_alignments( bamfile):
   #ref_len_str = !samtools view -c $bamfile
   #return int(ref_len_str[0])
   result = sp.check_output(['samtools view -c {0}'.format(bamfile)], shell=True)
   return int(result)

def get_ref_length(bamfile): 
    view = sp.check_output(["samtools view -H {0}".format(bamfile)], shell=True)
    return int(re.findall(r'.*LN:(.+)', view)[0])

def get_first_ref_seq(bamfile): 
    view = sp.check_output(["samtools view -H {0}".format(bamfile)], shell=True)
    return re.findall(r'.*SN:(.*)\t', view)[0]

def main():
    try: 
        bamfile, min_depth = sys.argv[1], int(sys.argv[2])
    except IndexError:
        #bamfile, min_depth =  '/home/AMED/michael.panciera/projects/samtools_primer/tutorial/alignments/sim_reads_aligned.sorted.bam', 1000
        #bamfile = '/home/AMED/michael.panciera/projects/780/780.bam'
        raise IndexError
    sys.stderr.write("bamfile: {0}, min_depth: {1}\n".format(bamfile, min_depth))
    ref_length = get_ref_length(bamfile)
    refid = get_first_ref_seq(bamfile)
    region_str = ':1-'.join([refid, str(ref_length)])
    matrix = DepthMatrix(ref_length, min_depth) 
    matrix.make_seq_matrix(bamfile, region_str)
    matrix.main()
    sampled_seqs = [seq.string for row in matrix.seq_matrix for seq in row if seq.picked]
    subsamplebam.make_subselected_bam(bamfile, sampled_seqs) 
    import ipdb; ipdb.set_trace()
    return matrix.depth_array

if __name__ == "__main__":
    main()
