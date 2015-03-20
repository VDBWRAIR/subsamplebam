from __future__ import print_function
import unittest
import sys
sys.path.append('/Users/wovenhead/clones/ngs_mapper/ngs_mapper/')

from  subsample_mindepth import *

class SimpleTest(unittest.TestCase):

    def setUp(self): 
        self.reads = [Alignment('', 0, 'A'*1),
                Alignment('', 0, 'A'*2),
                Alignment('', 0, 'A'*3)]
    def test_empty_get_depths(self):
        expected_depths = [3, 2, 1]
        result = get_depths(self.reads) 
        print(result)
        self.assertEquals(expected_depths, result)

    def test_get_candidate_sequences(self):
        # use the reads from setup and at position 0
        expected =  [Alignment('', 0, 'A'*1)]






