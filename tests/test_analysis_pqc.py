import unittest
import numpy as np

import analysis_pqc

x = np.array([0., 1., 2., 3., 4., 5.])
y = np.array([.1, .2, .3, .4, .5, .6])

class AnalysisPQCTest(unittest.TestCase):

    def test_states(self):
        self.assertEqual(analysis_pqc.STATUS_NONE, "none")
        self.assertEqual(analysis_pqc.STATUS_PASSED, "passed")
        self.assertEqual(analysis_pqc.STATUS_FAILED, "failed")

    def test_line_regr_with_cuts(self):
        r = analysis_pqc.line_regr_with_cuts(x, y, cut_param=-.5)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_iv(self):
        r = analysis_pqc.analyse_iv(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_cv(self):
        r = analysis_pqc.analyse_cv(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_FAILED)

    def test_analyse_mos(self):
        r = analysis_pqc.analyse_mos(x, y, cut_param=-1)
        self.assertEqual(r.status, analysis_pqc.STATUS_FAILED)

    def test_analyse_gcd(self):
        r = analysis_pqc.analyse_gcd(x, y, cut_param=-1)
        self.assertEqual(r.status, analysis_pqc.STATUS_FAILED)

    def test_analyse_fet(self):
        r = analysis_pqc.analyse_fet(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_van_der_pauw(self):
        r = analysis_pqc.analyse_van_der_pauw(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_cross(self):
        r = analysis_pqc.analyse_cross(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_linewidth(self):
        r = analysis_pqc.analyse_linewidth(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_cbkr(self):
        r = analysis_pqc.analyse_cbkr(x, y, r_sheet=-1, cut_param=-1)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)
        r = analysis_pqc.analyse_cbkr(x, y, r_sheet=1, cut_param=1e-5)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_contact(self):
        r = analysis_pqc.analyse_contact(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_meander(self):
        r = analysis_pqc.analyse_meander(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_breakdown(self):
        r = analysis_pqc.analyse_breakdown(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

    def test_analyse_capacitor(self):
        r = analysis_pqc.analyse_capacitor(x, y)
        self.assertEqual(r.status, analysis_pqc.STATUS_PASSED)

if __name__ == '__main__':
    unittest.main()
