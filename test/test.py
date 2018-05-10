import unittest
from context import mbwr

class TestEOSPressure(unittest.TestCase):

    def test_press1(self):
        self.assertAlmostEqual(mbwr.eos.mbwrP(2.64721E-03,7.22000E-01),1.86402E-03,delta=1.0E-06)

    def test_press2(self):
        self.assertAlmostEqual(mbwr.eos.mbwrP(1.39898E-02,8.94000E-01),1.13874E-02,delta=1.0E-06)
    def test_intEng(self):
        self.assertAlmostEqual(mbwr.eos.intrnEng(0.8,1.0),-5.533,delta=0.01)
    def test_intEng2(self):
        self.assertAlmostEqual(mbwr.eos.intrnEng(0.2,2.0),-1.306,delta=0.01)
if __name__ == '__main__':
        unittest.main()
