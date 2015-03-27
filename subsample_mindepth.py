from __future__ import print_function
import sys
import numpy as np
import subsamplebam
import subprocess as sp
import re

''' Python3 compatibility '''
try: 
    xrange 
except NameError: 
    xrange = range

class CommonEqualityMixin(object): 
        def __eq__(self, other):
            return (type(other) is type(self) and self.__dict__ == other.__dict__)

        def __ne__(self, other):
            return not self.__eq__(other)

'''
1. get all sequences which could overlap under_index, and which have not already been selected. (how get these, by using samtools view again? or store in memory? only need a UNIQUE identifier and the sequence length)
2. choose the sequence which has the longest length after under_index. This will deal with the situation where you picked a short sequence which didn't cover a coming index, before you knew you wouldn't need it.
'''

def get_raw_reads(bamfile, regionstr=""): 
   '''
   :param str bamfile: path to input file
   :param str regionstr: the seq id to select, followed by the region to look up, seperated by colon. i.e.:
   <refid>:1-1000
   :return a list of sequences from the output of `samtools view` 
   '''
   return list(subsamplebam.samview(bamfile,  regionstr))

def parse_alignment(raw_view): 
   '''
   :param str raw_view: a single sequence (line as displayed by `samtools view`
   :return: An Alignment object including the QNAME, POS, SEQ and whether or not it's an orphan.
   '''
   fields = raw_view.split()
   flag = int(fields[1])
   qname, pos, seq = fields[0], int(fields[3]), fields[9]
   return Alignment(raw_view, qname, pos, seq, flag)

def get_alignments(bamfile, regionstr): 
        return [parse_alignment(string) for string in get_raw_reads(bamfile, regionstr)]

class DepthMatrix(CommonEqualityMixin): 
    '''
    Keep track of the picked alignments (seq_matrix) and .depth_array, representing the coverage at each sequence position.
    '''

    def __init__(self, min_depth, allow_orphans=False):
        self.seq_matrix = [[]]
        #self.depth_array = np.zeros(ref_length)
        self.min_depth = min_depth 
        self.allow_orphans = allow_orphans


    def allow(self, seq):
        if self.allow_orphans:
            return True
        else:
            return not seq.orphan

    def get_candidate_sequences(self, under_index): 
        '''
        :param int under_index: the current position in the matrix which is under-covered.
        Candidates have not been used already, overlap the current index, and are not "orphan" or "anomalous" reads.
        orphaned/anomalous reads don't get spotted by default mpileup command. see https://github.com/VDBWRAIR/ngs_mapper/issues/112
        '''
        #TODO: have this only look backwards until reach under-covered index 
        start = max(under_index - self.max_seq_length, 0)
        sub_matrix = self.seq_matrix[start:under_index+1]
        '''flatten the matrix'''
        prev_seqs = [seq for row in sub_matrix for seq in row]
        return filter(lambda seq: seq.overlap  >= under_index  and self.allow(seq) and not seq.picked, prev_seqs)

    def yield_greatest_overlaps(self, under_index, num_needed): 
        '''
        @side-effect: sequences are set to "picked" when they are yielded by this function.
        :param int under_index: the position on the reference that is under minimum depth.
        :param int num_needed: the amount the overlap is under the minimum depth
        :return (yield) the sequence with the greatest overlap, which is a candidate.
        '''
        candidate_sequences = self.get_candidate_sequences(under_index) 
        matches = 0
        while matches < num_needed and candidate_sequences: 
            # could instead keep sequences in order sorted by overlap 
            farthest_overlap_seq = max(candidate_sequences, key=lambda seq: seq.overlap)
            candidate_sequences.remove(farthest_overlap_seq)
            farthest_overlap_seq.pick()
            matches += 1
            yield farthest_overlap_seq
    def pickreads(self, pos, needed_depth):
        ''' 
        1. get all sequences which could overlap under_index, and which have not already been selected. (how get these, by using samtools view again? or store in memory? only need a UNIQUE identifier and the sequence length)
        2. choose the sequence which has the longest length after under_index. This will deal with the situation where you picked a short sequence which didn't cover a coming index, before you knew you wouldn't need it.
        :param int pos: the current position which needs to be covered.
        :param int needed_depth: The amount of coverage needed to reach the minimum depth.
        :return True if the minimum depth was met at position pos, False otherwise.
        '''
        next_sequences = list(self.yield_greatest_overlaps(pos, needed_depth))
        if not next_sequences:  # no matches found
            return False
        ''' This doesn't work because it returns a flat list of depths starting at index 0, where in reality these reads
        Can start at different positions behind the current position.'''
        next_depths = self.get_depths(next_sequences)
        self.depth_array += next_depths
        #sys.stderr.write("depth {0} at pos {1}\n".format(str(len(next_sequences)), str(pos)))
        if  len(next_sequences) > needed_depth:
            raise Exception("Too many sequences retured")
        return len(next_sequences) == needed_depth

    def get_depths(self, reads): 
       ''' 
       :param list reads: Alignment objects that have been picked.
       :return: a numpy array representing the total depths or overlap of these alignments.
       '''
       #TODO:Needlessly creates a depth_array of equal size to self.depth_array. Also needlessly slow.
       depths = np.zeros(len(self.depth_array))
       #import ipdb; ipdb.set_trace()
       for seq in reads: 
           ''' Even if a sequence would overlap past the length of depth-array, we don't include that in depth-array. 
           the numpy arrays must be the same size in order to add them.  '''
           i, j = seq.pos, min(len(depths), seq.overlap)# +1)
           depths[i:j] += np.array([1 for overlap in xrange(i, j)])
       return depths

    def make_seq_matrix(self, bamfile, regionstr): 
        '''
        :param str bamfile: path to file
        :param str regionstr: reference sequence and length as listed in .bam header i.e.  <refname>:1-1000
        Parse a bam file using sam tools, and store the alignments as a 2d matrix, where each row is a posiiton in the reference sequence.  
        '''
        #import ipdb; ipdb.set_trace()
        all_alignments = get_alignments(bamfile, regionstr)
        max_pos = max([seq.pos for seq in all_alignments])
        max_overlap = max([seq.overlap for seq in all_alignments])
        self.depth_array = np.zeros(max_overlap)
        '''initialize sequence matrix as a 2D array in order of position.'''
        self.seq_matrix = [ [seq for seq in all_alignments if seq.pos == i] for i in xrange(max_pos+1)]
        self.max_seq_length = max([seq.seq_length for seq in all_alignments])

    def minimize_depths(self):
        ''' 
        Trim self.seq_matrix to minimize coverage overflow.
        For each position in the reference sequence, pick reads until the minimum depth is met if possible.
        '''
        #import ipdb; ipdb.set_trace()
        for pos, depth in enumerate(self.depth_array):
            if depth < self.min_depth:
                needed_depth = self.min_depth - depth
                depth_met = self.pickreads(pos, needed_depth)
                #sys.stderr.write("depth was {0}met.\n".format("" if depth_met else "NOT "))


class Alignment(CommonEqualityMixin):
    ''' Store info for the alignment, created from parsing the results of `samtools view`. See parse_alignment method. '''
    
    def __init__(self, string, qname, pos, seq, flag):
        self.qname, self.pos, self.seq = qname, pos, seq
        self.picked = False
        self.string = string
        self.orphan = self.is_orphan(flag) 

    ''' return true if the pair is not mapped. This is necessary to maintain consistency with `samtools mpileup`'''
    def is_orphan(self, flag): 
        return not flag & 0x2

    @property
    def seq_length(self):
        return len(self.seq)

    @property
    def overlap(self):
        return self.pos + self.seq_length

    def pick(self):
        self.picked = True

    def __repr__(self):
        #return "pos: {0}, overlap {1} \n seq {2} \n {3}".format(self.pos, self.overlap, self.seq, self.string)
        return "pos: {0}, overlap {1} ".format(self.pos, self.overlap)

def flatten_and_filter_matrix(matrix):
    return [seq.string for row in matrix for seq in row if seq.picked]

def main(): 
    '''
    Parse the sam view, initialize the DepthMatrix, and trim the sequence matrix using minimize_depths.
    Then, ouptut the results (as a sam view) to stdout.
    '''
    #TODO: need to handle case where fasta file does not span entire reference genome; depth array should be limited in that case to the largest overlap
    args = subsamplebam.parse_args(wrapper=False)
    sys.stderr.write(str(args)+'\n') 
#    ref_length = int(args.regionstr.split(':')[-1].split('-')[-1])
#    region_str = args.regionstr.split(':')[0]
    matrix = DepthMatrix(args.subsample, allow_orphans=args.count_orphans) 
    matrix.make_seq_matrix(args.bamfile, args.refseq)
    matrix.minimize_depths()
    ''' Flatten the matrix '''
    sampled_seqs = flatten_and_filter_matrix(matrix.seq_matrix)
    '''Print results to stdout as a sam view (.sam file)'''
    subsamplebam.make_subselected_bam(args.bamfile, sampled_seqs) 

if __name__ == "__main__":
    main()
