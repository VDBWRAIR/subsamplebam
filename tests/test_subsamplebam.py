from os.path import *
from io import StringIO

from nose.plugins.attrib import attr

from compat import (
    unittest,
    mock
)

import subsamplebam

THISD = dirname(abspath(__file__))

@mock.patch('subsamplebam.sys')
@mock.patch('subsamplebam.subprocess.Popen')
class TestSamtoolsIsAvailable(unittest.TestCase):
    def test_missing_samtools_executable_returns_false(self, mpopen, msys):
        mpopen.side_effect = OSError('samtools missing')
        r = subsamplebam.samtools_is_available()
        self.assertFalse(r)

    def test_samtools_returns_incorrect_string_returns_false(self, mpopen, msys):
        mpopen.return_value.communicate.return_value = ('stdout','')
        r = subsamplebam.samtools_is_available()
        self.assertFalse(r)

    def test_samtools_executes_and_has_version_returns_true(self, mpopen, msys):
        mpopen.return_value.communicate.return_value = ('Version: 0.1.19+','')
        r = subsamplebam.samtools_is_available()
        self.assertTrue(r)

@mock.patch('subsamplebam.sys')
@mock.patch('subsamplebam.subprocess.Popen')
class TestSamview(unittest.TestCase):
    def test_returns_samtools_view_output_as_subprocess_stdout(self, mpopen, msys):
        samlines = ['line1', 'line2', 'line3']
        mock_popen = mock.MagicMock()
        mock_popen.poll.return_value = None
        mock_popen.stdout.readline.return_value = samlines[0]
        mock_popen.stdout.__iter__.return_value = samlines[1:]

        mpopen.return_value = mock_popen

        r = subsamplebam.samview('foo.bam', 'chr1:1-100')
        cmd = 'samtools view foo.bam chr1:1-100'
        for eline, rline in zip(samlines, r):
            self.assertEqual(eline, rline)
        mpopen.assert_called_once_with(cmd.split(), stdout=-1)
        msys.stderr.write.assert_called_with(cmd + '\n')

    def test_samtools_returns_non0_raises_valueerror(self, mpopen, msys):
        mock_popen = mock.MagicMock()
        mock_popen.poll.return_value = 1
        mock_popen.stdout.readline.return_value = ''
        mock_popen.stdout.__iter__.return_value = []

        mpopen.return_value = mock_popen

        r = subsamplebam.samview('foo.bam', 'chr1:1-100')
        self.assertRaises(ValueError, next, r)

    def test_only_one_line_samtools_view_no_error(self, mpopen, msys):
        mock_popen = mock.MagicMock()
        mock_popen.poll.return_value = 0
        mock_popen.stdout.readline.return_value = 'line1'
        mock_popen.stdout.__iter__.return_value = []

        mpopen.return_value = mock_popen

        r = subsamplebam.samview('foo.bam', 'chr1:1-100')
        lines = list(r)
        self.assertEqual('line1', lines[0])

    def test_0_line_samtools_view_no_error(self, mpopen, msys):
        mock_popen = mock.MagicMock()
        mock_popen.poll.return_value = 0
        mock_popen.stdout.readline.return_value = ''
        mock_popen.stdout.__iter__.return_value = []

        mpopen.return_value = mock_popen

        r = subsamplebam.samview('foo.bam', 'chr1:1-100')
        lines = list(r)
        self.assertEqual([''], lines)

@mock.patch('subsamplebam.multiprocessing.Pool')
@mock.patch('subsamplebam.sys')
@mock.patch('subsamplebam.subprocess.Popen')
class TestMakeSubselectBam(unittest.TestCase):
    def setUp(self):
        self.uniquereads = (
            'read1  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
        )
        self.headers = (
            '@HD VN:1.3  SO:coordinate',
            '@SQ SN:chr1 LN:10',
            '@RG ID:MiSeq    SM:foo CN:None PL:ILLUMINA',
        )

    def test_contains_bam_header_and_samviewrows(self, mpopen, msys, mpool):
        mpool.return_value.map = map
        mpopen.return_value.stdout = self.headers
        r = subsamplebam.make_subselected_bam('foo.bam', self.uniquereads)
        msys.stdout.write.assert_has_calls([
            mock.call(self.headers[0]),
            mock.call(self.headers[1]),
            mock.call(self.headers[2]),
            mock.call(self.uniquereads[0])
        ])
    
    def test_uniquereads_empty_sequence(self, mpopen, msys, mpool):
        mpopen.return_value.stdout = self.headers
        r = subsamplebam.make_subselected_bam('foo.bam', self.uniquereads)
        ecalls = []
        for c in self.headers:
            ecalls.append(mock.call(c))
        msys.stdout.write.assert_has_calls(ecalls)

class TestReferenceInfo(unittest.TestCase):
    def setUp(self):
        self.mock_seqrecords = [
            mock.Mock(
                id='chr1',
                seq='ATGCATGCAT',
                description='chr1 descr'
            ),
            mock.Mock(
                id='chr2',
                seq='ATGC',
                description='chr2 descr'
            )
        ]

    @mock.patch('subsamplebam.SeqIO.parse')
    def test_returns_dictionary_by_id_values_are_lengths(self, m_seqio_parse):
        m_seqio_parse.return_value = self.mock_seqrecords
        r = subsamplebam.reference_info('reference.fasta')
        for rec in self.mock_seqrecords:
            seqlen = len(rec.seq)
            self.assertEqual(seqlen, r[rec.id])

@mock.patch('subsamplebam.SeqIO.parse')
@mock.patch('subsamplebam.multiprocessing.Pool')
@mock.patch('subsamplebam.subprocess.Popen')
class TestRandomlySelectsReads(unittest.TestCase):
    def setUp(self):
        self.samview_lines = (
            'read1  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read2  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read3  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read4  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read5  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read6  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read7  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read8  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read9  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
            'read10  1  chr1    1   60  10M    =   1  1    TTTCGAATC    FFFFFFFFF    NM:i:3  AS:i:231    XS:i:0  RG:Z:MiSeq',
        )

    def test_randomly_selects_correct_amount_of_items(self, mpopen, mpool, mparse):
        mpopen.return_value.poll.return_value = None
        mpopen.return_value.stdout.__iter__.return_value = self.samview_lines[1:]
        mpopen.return_value.stdout.readline.return_value = self.samview_lines[0]
        mpool.return_value.map = map

        with mock.patch('subsamplebam.sys') as msys:
            r = subsamplebam.subselect_from_bam(
                'foo.bam', 10, 'chr1:1-10',
            )

            self.assertEqual(set(self.samview_lines), r)

    def test_handles_missing_depth_error(self, mpopen, mpool, mparse):
        mpopen.return_value.poll.return_value = None
        mpopen.return_value.stdout.__iter__.return_value = self.samview_lines[1:]
        mpopen.return_value.stdout.readline.return_value = self.samview_lines[0]
        mpool.return_value.map = map

        with mock.patch('subsamplebam.sys') as msys:
            r = subsamplebam.subselect_from_bam(
                'foo.bam', 11, 'chr1:1-10',
            )
            self.assertEqual(set([]), r)
            for i in range(1,11):
                msys.stderr.write.assert_has_call('Depth for chr1:{0}-{0} is only 10\n'.format(i))

@mock.patch('subsamplebam.multiprocessing.Pool')
@mock.patch('subsamplebam.sys')
@mock.patch('subsamplebam.argparse.ArgumentParser')
class TestMain(unittest.TestCase):
    def test_correct_args(self, margparser, msys, mpool):
        mpool.return_value.map = map
        args = mock.Mock()
        margparser.return_value.parse_args.return_value = args
        
        args.bamfile = join(THISD, 'test.bam')
        args.subsample = 10
        args.regionstr = 'chr1:1-10'
        subsamplebam.main()
