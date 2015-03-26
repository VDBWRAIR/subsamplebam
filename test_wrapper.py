import unittest
import wrapper
from wrapper import subsample_reference, subsample_file
class TestWrapper(unittest.TestCase):
    def setUp(self):
        self.ecoli = {'ref' : 'gi|110640213|ref|NC_008253.1|' , 
                'length'  : 4938920}

    def test_parse_header(self): 
        bamfile = 'ecoli.bam'
        expected = {'gi|110640213|ref|NC_008253.1|' : 4938920}
        result = wrapper.parse_header(bamfile)
        self.assertEqual(expected, result)

    def test_parse_header_multiple(self):
        bamfile = '780.bam'
        expected_seqs = 8
        result = len(wrapper.parse_header(bamfile))
        self.assertEqual(expected_seqs, result)



    def check_output_bam(self, refname, depth, stdout=None): 
        expected_depth_file = open('expected/{0}.min.{1}.bam'.format(refname, str(depth)), 'r').read()
        filename = "{0}.min.{1}.bam".format(refname, str(depth))
        if not stdout: stdout = open(filename).read()
        actual = stdout
        #actual = stdout or open(filename).read()
        self.assertEqual(expected, actual)

    def test_subsample_reference_over_max_depth(self): 
        depth = 1000
        refname, length = 'gi|110640213|ref|NC_008253.1|' , 4938920
        stdout = subsample_reference(refname, length, 'ecoli.bam', depth, allow_orphans=True)
        self.check_output_bam(refname, depth, stdout=stdout)

    def test_subsample_reference_minimize_depth(self): 
        depth = 30
        refname, length = self.ecoli['ref'], self.ecoli['length']
        stdout = subsample_reference(refname, length, 'ecoli.bam', depth, allow_orphans=True)
        self.check_output_bam(refname, depth, stdout=stdout)

    def test_subsample_reference_minimize_depth(self): 
        depth = 200
        refname, length ='H3N2/CY074915/Managua/2010/HA_4', 1731
        bamfile = '../../780.bam'
        stdout = subsample_reference(refname, length, bamfile, depth, allow_orphans=False)
        self.check_output_bam(refname, depth, stdout=stdout)

    def test_subsample_file_ecoli(self):
        depth, bamfile  = 40, 'ecoli.bam'
        subsample_file(bamfile, depth=depth, allow_orphans=True)
        self.check_output_bam(self.ecoli['ref'], depth)
        

        




