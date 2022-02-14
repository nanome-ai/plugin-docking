import os
import unittest
import sys

from nanome.util import Logs

docking_dir = f'{os.getcwd()}/nanome_docking/'
sys.path.append(docking_dir)

test_directory = 'tests'
suite = unittest.TestLoader().discover(test_directory)

output = unittest.TextTestRunner(verbosity=3).run(suite)

if output.wasSuccessful():
    sys.exit(0)
else:
    sys.exit(1)
