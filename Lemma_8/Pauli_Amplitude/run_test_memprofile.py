# run_test_memprofile.py

import unittest
from Lemma_8.test_prob_dist import TestProbDist

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(TestProbDist("test_stat_measures"))
    unittest.TextTestRunner().run(suite)
