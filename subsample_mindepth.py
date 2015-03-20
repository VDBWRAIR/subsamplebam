from __future__ import print_function
import sys
#import pandas as pd 
import numpy as np
import subsamplebam

''' bugs '''

def get_raw_reads(bamfile, refseq=None, idx=None):
 #  if not (refseq or idx):
 #      res = !samtools view $bamfile
 #  else:
 #     region_str = "{0}:{1}-{1}".format(refseq, idx)
 #     res = !samtools view $bamfile "$region_str"
 #  return res
 return list(subsamplebam.samview(bamfile, ''))
 #pass

def add_arrays(a, b, add_idx):
   bigger, smaller = max(a, b, key=len), min(a, b, key=len)
   result = bigger.copy()
   result[add_idx:add_idx+len(smaller)] += smaller
   return result


def get_next_reads(reads, pos):
   #def is_at_pos(pos): return read['pos'] == pos
   return [read for read in reads if read['pos'] == pos]


    
'''
Only add a read if the current depth is less than the target depth.
Unfortunately this doesn't work because you can be under the target depth at a position that has no alignments starting there.
'''



def parse_alignment(raw_view):
   fields = raw_view.split()
   #return dict(zip(['qname', 'pos', 'seq'], [fields[0], int(fields[3]), fields[9]]))
   qname, pos, seq = fields[0], int(fields[3]), fields[9]
   return Alignment(qname, pos, seq)
   #return dict(zip(['qname', 'pos', 'seq'], [fields[0], int(fields[3]), fields[9]]))




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
        #return [filter(lambda seq: seq.overlap  >= under_index and not seq.picked, seqs) for seqs in self.seq_matrix[:under_index+1] ]
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
        print("depth {0} at pos {1}".format(str(len(next_sequences)), str(pos)))
        if  len(next_sequences) > needed_depth:
            raise Exception("Too many sequences retured")
        return len(next_sequences) == needed_depth
    ''' 
    Needlessly creates a depth_array of equal size to self.depth_array 
    '''
    def get_depths(self, reads, min_pos): 
       #def seq_len(alignment): return len(alignment['seq'])
       #return [len(filter(lambda x: seq_len(x) > i, reads)) for i in xrange(max(map(seq_len, reads))) ]
       #return [len(filter(lambda seq: seq.seq_length > i, reads)) for i in xrange(max([seq.seq_length for seq in reads])) ]
       #return [len(seq) for i in xrange(max([seq.seq_length for seq in reads])) ]
       #return [len(filter(lambda x: x.seq_length > i, reads)) for i in xrange(max([s.seq_length for s in reads])) ] 
       depths = np.zeros(len(self.depth_array))
       for seq in reads: 
           ''' Even if a sequence would overlap past the length of depth-array, we don't include that in depth-array. 
           the numpy arrays must be the same size in order to add them.
           '''
           i, j = seq.pos, min(len(depths), seq.overlap +1)
           depths[i:j] += np.array([1 for overlap in xrange(i, j)])
       return depths

    def print_matrix(self):
        for row in self.seq_matrix:
            print( [seq.pos for seq in row])

    def make_seq_matrix(self, bamfile): 
        all_alignments = map(parse_alignment, get_raw_reads(bamfile))
        max_pos = max([seq.pos for seq in all_alignments])#max(all_alignments, key=lambda seq: seq.pos)
        #initialize sequence matrix as a 2D array
        # weird thing about pos 0 
        seq_matrix = [ filter( lambda seq: seq.pos == i, all_alignments) for i in xrange(1, max_pos)]
        self.seq_matrix = seq_matrix
        ''' does htis work lol?'''
        return
        #self.seq_matrix = [ [seq for seq in all_alignments] for i in xrange(max_pos) if seq.pos == i]
        #self.print_matrix()

# need to handle case where fasta file does not span entire reference genome;
# depth array should be limited in that case to the largest overlap
    def main(self):
        for pos, depth in enumerate(self.depth_array):
            #print( self.depth_array )
            if depth < self.min_depth:
                needed_depth = self.min_depth - depth
                depth_met = self.backtrack(pos, needed_depth)
                print("depth was {0}met.".format("" if depth_met else "NOT "))


def get_ref_length( bamfile):
   #ref_len_str = !samtools view -c $bamfile
   #return int(ref_len_str[0])
   return 1000


def main():
    try: 
        bamfile, min_depth = sys.argv[1], int(sys.argv[2])
    except IndexError:
        bamfile, min_depth =  '/Users/wovenhead/clones/samtools_primer/tutorial/alignments/sim_reads_aligned.sorted.bam', 58
        bamfile = '/Users/wovenhead/clones/ngs_mapper/tdir/780/780.bam'
    print("bamfile: {0}, min_depth: {1}".format(bamfile, min_depth))
    ref_length = get_ref_length(bamfile)
    matrix = DepthMatrix(ref_length, min_depth) 
    matrix.make_seq_matrix(bamfile)
    matrix.main()
    print(matrix.depth_array)
    return matrix.depth_array

class Alignment(object):

    def __init__(self, qname, pos, seq):
        self.qname, self.pos, self.seq = qname, pos, seq
        self.picked = False

    @property
    def seq_length(self):
        return len(self.seq)

    @property
    def overlap(self):
        return self.pos + self.seq_length

    def pick(self):
        self.picked = True


if __name__ == "__main__":
    main()
