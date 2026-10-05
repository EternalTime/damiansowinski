"""The suite runs in the environment _tools/setup-env.sh builds, or not at all.

environment.py says what that environment is. Its check runs here as the suite is gathered,
before any test, and on a problem prints them all with the command that mends them and ends
the process: unittest turns an exception raised while gathering into one more failing test
among hundreds of skips, which is how a broken environment once went unnoticed.
"""
import os
import sys
import unittest

import environment

FOUND = environment.problems()
if FOUND:
    sys.stdout.flush()
    print(environment.report(FOUND), file=sys.stderr, flush=True)
    os._exit(1)


class Environment(unittest.TestCase):
    def test_the_environment_is_whole(self):
        self.assertEqual(environment.problems(), [])


if __name__ == "__main__":
    unittest.main()
