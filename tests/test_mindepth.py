from __future__ import print_function
import unittest
import sys 
from compat import (
        unittest,
        mock
        )
#sys.path.append('/Users/wovenhead/clones/ngs_mapper/ngs_mapper/')
import numpy as np
import subsample_mindepth as sub
from subsample_mindepth import Alignment, DepthMatrix
from numpy.ma.testutils import assert_equal
from argparse import Namespace 

from past.builtins import map, xrange, filter
#from six import StringIO
#try:
#        from StringIO import StringIO
#except ImportError:
#        from io import StringIO
#
from io import BytesIO
import os.path


THISD = os.path.dirname(os.path.abspath(__file__))

def mock_args(wrapper=False): 
        bamfile = os.path.join(THISD, 'ecoli.bam')
        return Namespace(reflength=1000, subsample=40, count_orphans=True, bamfile=bamfile, refseq="gi|110640213|ref|NC_008253.1|")


def mock_get_raw_reads( bamfile, regionstr=''):
     return ''' gi|110640213|ref|NC_008253.1|_2_451_1:0:0_1:0:0_12b/1	0	gi|110640213|ref|NC_008253.1|	0	42	70M	*	0	0	GCT 222 AS:i:-3	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:45A24	YT:Z:UU
10640213|ref|NC_008253.1|_4_510_1:0:0_2:0:0_114/1	0	gi|110640213|ref|NC_008253.1|	0	42	70M	*	0	0	TTT	222	AS:i:-3	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:26A43	YT:Z:UU
10640213|ref|NC_008253.1|_5_452_2:0:0_1:0:0_1d/1	0	gi|110640213|ref|NC_008253.1|	0	42	70M	*	0	0	TT	22	AS:i:-6	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:31C23C14	YT:Z:UU
10640213|ref|NC_008253.1|_6_397_0:0:0_0:0:0_2c1/1	0	gi|110640213|ref|NC_008253.1|	0	42	70M	*	0	0	TTTC	222	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YT:Z:UU
10640213|ref|NC_008253.1|_7_500_3:0:0_1:0:0_250/1	0	gi|110640213|ref|NC_008253.1|	1	40	70M	*	0	0	TCATTC	222222	AS:i:-9	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:9C27A5A26	YT:Z:UU
10640213|ref|NC_008253.1|_9_557_0:0:0_2:0:0_13d/1	0	gi|110640213|ref|NC_008253.1|	1	42	70M	*	0	0	AT	222	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YT:Z:UU
10640213|ref|NC_008253.1|_10_554_1:0:0_1:0:0_2a/1	0	gi|110640213|ref|NC_008253.1|   2		42	70M	*	0	0	TTCTGA	222222	AS:i:-3	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:15G54	YT:Z:UU
10640213|ref|NC_008253.1|_11_515_1:0:0_1:0:0_335/1	0	gi|110640213|ref|NC_008253.1|	2	42	70M	*	0	0	TCT	222	AS:i:-3	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:10A59	YT:Z:UU'''.split('\n')
    
def short_mock_get_raw_reads( a, b):
    "X 	0	X	1	42	70M	*	0	0	ATT	222	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:70	YT:Z:UU"

def short_mock_get_alignments( a, b): 
    return [Alignment('', '', 0, 'A'*1, 0x2),
            Alignment('', '', 0, 'A'*2, 0x2),
            Alignment('', '', 0, 'A'*3, 0x2)]


class SimpleTest(unittest.TestCase):

    #@mock.patch("subsample_mindepth.get_raw_reads", side_effect=short_mock_get_raw_reads)
    @mock.patch("subsample_mindepth.get_alignments", side_effect=short_mock_get_alignments)
    def setUp(self, get_raw_reads_function): 
        #self.raw_reads = self.mock_get_raw_reads("", "")
        self.reads = short_mock_get_alignments('', '')
        self.raw_reads = mock_get_raw_reads('', '')
        self.matrix = DepthMatrix(3, allow_orphans=True) 
        self.matrix.make_seq_matrix("", "")

    @mock.patch("subsample_mindepth.get_raw_reads", side_effect=mock_get_raw_reads)
    def complexSetUp(self, func):
        self.matrix.make_seq_matrix('', '')
        pass

    

    def test_simple_get_depths(self):
        expected_depths = np.array([3, 2, 1])
        result = self.matrix.get_depths(self.reads) 
        #print(result)
        assert_equal(expected_depths, result)

    def test_get_candidate_sequences(self):
        # use the reads from setup and at position 0
        expected = self.reads
        result = self.matrix.get_candidate_sequences(0)
        self.assertEquals(expected, result)

    def test_get_candidate_sequences_no_candidates(self):
        result = self.matrix.get_candidate_sequences(5)
        self.assertEquals([], result)


    def test_yield_greatest_overlaps(self):
        expected = [Alignment('', '', 0, 'A'*1, 0x2)]


    def test_parse_alignment(self):
        expected = Alignment(self.raw_reads[0], "gi|110640213|ref|NC_008253.1|_2_451_1:0:0_1:0:0_12b/1", 0, "GCT", 0)
        result = sub.parse_alignment(self.raw_reads[0]) 
        self.assertEquals(expected, result)

    @mock.patch("subsample_mindepth.get_raw_reads", side_effect=mock_get_raw_reads)
    def test_make_seq_matrix(self, get_reads_function):
        expected_lengths = [4, 2, 2]
        self.matrix.make_seq_matrix("mockecoli.sam", "gi|110640213|ref|NC_008253.1|")
        actual_lengths = map(len, self.matrix.seq_matrix)
        self.assertEquals(expected_lengths, actual_lengths)

    def test_get_raw_reads(self):
        pass

    def test_yield_greatest_overlaps_ecoli(self):
        self.complexSetUp()
        expected = [sub.parse_alignment(self.raw_reads[-2]),  sub.parse_alignment(self.raw_reads[4])]
        for seq in expected:
            seq.pick()
        result = list(self.matrix.yield_greatest_overlaps(3, 2))
        self.assertEquals(expected, result)
        #expected = self.matrix.seq_matrix[2][0] 
        pass

    def test_pickreads(self):
        pass

    def test_minimize_depths_ecoli(self):
        self.complexSetUp()
        seqs = map(sub.parse_alignment, self.raw_reads) 
        expected_matrix =  [
                [seqs[0],  seqs[1],  seqs[3]],
                [seqs[4]],
                [seqs[6], seqs[7]]
                ]
        expected_list = [  seqs[0], seqs[1], seqs[3], seqs[4], seqs[6], seqs[7] ]
        for seq in expected_list:
            seq.pick()
        expected_depths = np.array([3, 4, 6, 4, 3, 2, 2, 1])
        self.matrix.min_depth = 3
        self.matrix.minimize_depths()
        actual_matrix = self.matrix.seq_matrix
        for row in actual_matrix:
            for seq in row:
                if not seq.picked:
                    row.remove(seq)

        self.assertEquals(expected_matrix, actual_matrix)
        assert_equal(expected_depths, self.matrix.depth_array)
        self.assertEquals([seq.string for seq in expected_list], sub.flatten_and_filter_matrix(self.matrix.seq_matrix)) 




    @mock.patch("subsamplebam.parse_args", side_effect=mock_args)
    @mock.patch('subsamplebam.sys.stdout', new_callable=BytesIO)
    def test_main(self, mock_stdout, func):
        sub.main() 

        samfile = os.path.join(THISD, 'test40.sam')
        expected = open(samfile, 'r').read().strip()
        result = str(mock_stdout.getvalue().strip().decode('UTF-8'));
        self.assertEquals(expected, result)







