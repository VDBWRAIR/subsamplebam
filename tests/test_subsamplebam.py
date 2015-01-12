from nose.plugins.attrib import attr

from compat import (
    unittest,
    mock
)

import subsamplebam

@mock.patch('subsamplebam.subprocess.Popen')
class TestSamtoolsIsAvailable(unittest.TestCase):
    def test_missing_samtools_executable_returns_false(self, mpopen):
        mpopen.side_effect = OSError('samtools missing')
        r = subsamplebam.samtools_is_available()
        self.assertFalse(r)

    def test_samtools_returns_incorrect_string_returns_false(self, mpopen):
        mpopen.return_value.communicate.return_value = ('stdout','')
        r = subsamplebam.samtools_is_available()
        self.assertFalse(r)

    def test_samtools_executes_and_has_version_returns_true(self, mpopen):
        mpopen.return_value.communicate.return_value = ('Version: 0.1.19+','')
        r = subsamplebam.samtools_is_available()
        self.assertTrue(r)
