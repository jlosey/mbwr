import unittest
from context import mbwr

class TestEOSPressure(unittest.TestCase):

    def test_press1(self):
        self.assertAlmostEqual(mbwr.eos.mbwrP(2.64721E-03,7.22000E-01),1.86402E-03,delta=1.0E-06)

    def test_press2(self):
        self.assertAlmostEqual(mbwr.eos.mbwrP(1.39898E-02,8.94000E-01),1.13874E-02,delta=1.0E-06)

if __name__ == '__main__':
        unittest.main()
