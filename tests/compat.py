import sys

try:
    import unittest2 as unittest
    unittest # silence pyflakes
except ImportError:
    import unittest
    unittest # silence pyflakes

try:
    import mock
    mock # silence pyflakes
except ImportError:
    import unittest.mock as mock
    mock # silence pyflakes

PY3 = sys.version_info[0] == 3
if PY3:
    builtins_name = "builtins"
else:
    builtins_name = "__builtin__"
