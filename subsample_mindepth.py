from __future__ import print_function
import sys
import numpy as np
import subsamplebam
import subprocess as sp
import re

'''
1. get all sequences which could overlap under_index, and which have not already been selected. (how get these, by using samtools view again? or store in memory? only need a UNIQUE identifier and the sequence length)
2. choose the sequence which has the longest length after under_index. This will deal with the situation where you picked a short sequence which didn't cover a coming index, before you knew you wouldn't need it.
'''

'''
@param bamfile: input file
@param regionstr: the seq id to select, followed by the region to look up, seperated by colon. i.e.:
<refid>:1-1000
@return a list of sequences from the output of `samtools view` 
'''
def get_raw_reads(bamfile, regionstr=""):
   return list(subsamplebam.samview(bamfile,  regionstr))

'''
@param raw_view: a single sequence (line as displayed by `samtools view`
@return: An Alignment object including the QNAME, POS, SEQ and whether or not it's an orphan.
'''
def parse_alignment(raw_view):
   fields = raw_view.split()
   flag = int(fields[1])
   qname, pos, seq = fields[0], int(fields[3]), fields[9]
   return Alignment(raw_view, qname, pos, seq, flag)

'''
Keep track of the picked alignments (seq_matrix) and .depth_array, representing the coverage at each sequence position.
'''
class DepthMatrix(object):
    def __init__(self, ref_length, min_depth, allow_orphans=False):
        self.seq_matrix = [[]]
        self.depth_array = np.zeros(ref_length)
        self.min_depth = min_depth 
        self.allow_orphans = allow_orphans

    '''
    Candidates have not been used already, overlap the current index, and are not "orphan" or "anomalous" reads.
    orphaned/anomalous reads don't get spotted by default mpileup command. see https://github.com/VDBWRAIR/ngs_mapper/issues/112
    '''
    def get_candidate_sequences(self, under_index):
        #TODO: have this only look backwards until reach under-covered index 
        start = min(under_index - self.max_seq_length, 0)
        sub_matrix = self.seq_matrix[start:under_index+1]
        '''flatten the matrix'''
        prev_seqs = [seq for row in sub_matrix for seq in row]
        return filter(lambda seq: seq.overlap  >= under_index and not seq.picked and not (seq.orphan and not self.allow_orphans), prev_seqs)
    '''
    @side-effect: sequences are set to "picked" when they are yielded by this function.
    @param under_index: the position on the reference that is under minimum depth.
    @param num_needed: the amount the overlap is under the minimum depth
    @yield: the sequence with the greatest overlap, which is a candidate.
    '''
    def yield_greatest_overlaps(self, under_index, num_needed):
        candidate_sequences = self.get_candidate_sequences(under_index) 
        matches = 0
        while matches < num_needed and candidate_sequences: 
            # could instead keep sequences in order sorted by overlap 
            farthest_overlap_seq = max(candidate_sequences, key=lambda seq: seq.overlap)
            candidate_sequences.remove(farthest_overlap_seq)
            farthest_overlap_seq.pick()
            matches += 1
            yield farthest_overlap_seq
    ''' 
    2. get all sequences which could overlap under_index, and which have not already been selected. (how get these, by using samtools view again? or store in memory? only need a UNIQUE identifier and the sequence length)
    2. choose the sequence which has the longest length after under_index. This will deal with the situation where you picked a short sequence which didn't cover a coming index, before you knew you wouldn't need it.
    @param pos: the current position which needs to be covered.
    @param needed_depth: The amount of coverage needed to reach the minimum depth.
    @return True if the minimum depth was met at position pos, False otherwise.
    '''
    def pickreads(self, pos, needed_depth):
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

    ''' 
    @param reads: Alignment objects that have been picked.
    @return: a numpy array representing the total depths or overlap of these alignments.
    '''
    def get_depths(self, reads): 
       #TODO:Needlessly creates a depth_array of equal size to self.depth_array. Also needlessly slow.
       depths = np.zeros(len(self.depth_array))
       for seq in reads: 
           ''' Even if a sequence would overlap past the length of depth-array, we don't include that in depth-array. 
           the numpy arrays must be the same size in order to add them.  '''
           i, j = seq.pos, min(len(depths), seq.overlap +1)
           depths[i:j] += np.array([1 for overlap in xrange(i, j)])
       return depths

    ''' for debugging'''
    def print_matrix(self):
        for row in self.seq_matrix:
            sys.stderr.write( [seq.pos for seq in row])

    '''
    Parse a bam file using sam tools, and store the alignments as a 2d matrix, where each row is a posiiton in the reference sequence.  
    '''
    def make_seq_matrix(self, bamfile, refseq): 
        all_alignments = [parse_alignment(string) for string in get_raw_reads(bamfile, refseq)]
        max_pos = max([seq.pos for seq in all_alignments])
        '''initialize sequence matrix as a 2D array in order of position.'''
        self.seq_matrix = [ [seq for seq in all_alignments if seq.pos == i] for i in xrange(max_pos)]
        self.max_seq_length = max([seq.seq_length for seq in all_alignments])

    ''' 
    Trim self.seq_matrix to minimize coverage overflow.
    For each position in the reference sequence, pick reads until the minimum depth is met if possible.
    '''
    def minimize_depths(self):
        #import ipdb; ipdb.set_trace()
        for pos, depth in enumerate(self.depth_array):
            if depth < self.min_depth:
                needed_depth = self.min_depth - depth
                depth_met = self.pickreads(pos, needed_depth)
                #sys.stderr.write("depth was {0}met.\n".format("" if depth_met else "NOT "))


''' Store info for the alignment, created from parsing the results of `samtools view`. See parse_alignment method. '''
class Alignment(object):

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

    def __str__(self):
        return "pos: {0}, overlap {1} \n seq {2} \n {3}".format(self.pos, self.overlap, self.seq, self.string)

''' Not used '''
def get_num_alignments( bamfile):
   #ref_len_str = !samtools view -c $bamfile
   result = sp.check_output(['samtools view -c {0}'.format(bamfile)], shell=True)
   return int(result)

''' 
@return the length of the first reference sequence. 
'''
#TODO: Make this accept the sequence ID so it can return the correct length, not just the first one.
def get_ref_length(bamfile): 
    view = sp.check_output(["samtools view -H {0}".format(bamfile)], shell=True)
    return int(re.findall(r'.*LN:(.+)', view)[0])
'''
@return the ref idea of the first sequence listed in the .sam header. 
'''
#TODO: Replace this function with a command-line argument to allow parallel running for each ref sequence.
def get_first_ref_seq(bamfile): 
    view = sp.check_output(["samtools view -H {0}".format(bamfile)], shell=True)
    return re.findall(r'.*SN:(.*)\t', view)[0]

'''
Parse the sam view, initialize the DepthMatrix, and trim the sequence matrix using minimize_depths.
Then, ouptut the results (as a sam view) to stdout.
'''
def main():
    #TODO: expect the region string as a command-line argument, or parse the header to get all sequences.
    #TODO: need to handle case where fasta file does not span entire reference genome; depth array should be limited in that case to the largest overlap
    try:
        args = sys.argv[1:]
        sys.stderr.write(str(args))
        if "-A" in args:
            allow_orphans = True
            args.remove("-A")
            bamfile, min_depth = args[0], int(args[1])
        else:
            allow_orphans = False
            #bamfile, min_depth = sys.argv[1], int(sys.argv[2])
    except IndexError:
        raise IndexError("Not enough arguments")
    sys.stderr.write("bamfile: {0}, min_depth: {1}\n".format(bamfile, min_depth))
    ref_length = get_ref_length(bamfile)
    refid = get_first_ref_seq(bamfile)
    region_str = ':1-'.join([refid, str(ref_length)])
    matrix = DepthMatrix(ref_length, min_depth, allow_orphans=allow_orphans) 
    matrix.make_seq_matrix(bamfile, region_str)
    matrix.minimize_depths()
    ''' Flatten the matrix '''
    sampled_seqs = [seq.string for row in matrix.seq_matrix for seq in row if seq.picked]
    '''Print results to stdout as a sam view (.sam file)'''
    subsamplebam.make_subselected_bam(bamfile, sampled_seqs) 

if __name__ == "__main__":
    main()
