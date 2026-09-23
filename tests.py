import unittest
from security import SETS, bound, add
from run import inputs, ROOT

class Checks(unittest.TestCase):
    def test_log_sum(self):
        self.assertAlmostEqual(add([-1000,-1000]),-999)

    def test_design_results(self):
        expected=[131.0347706951663,194.4954373857997,260.65464487358156]
        for s,e in zip(SETS,expected):
            b,dp,dq,_=inputs(ROOT/'data/raw',s)
            r=bound(s,s['design'],b,dp,dq)
            self.assertAlmostEqual(r['security_bits'],e,places=8)
            self.assertGreater(r['security_bits']-s['target'],2)
            self.assertIsNone(bound(s,1,b,dp,dq)['outer_log2_terms']['salt_bad'])

    def test_frodo_direct_equation(self):
        # Independent ordinary-domain evaluation at N=1.
        for s in SETS:
            t,q,M=2**32,2**14,2**s['fmu']
            delta=2**s['fd']
            a=6*q/M+2**-512+1/M+q*delta+2*t*2**-s['fb']
            import math
            samples=16*s['frodo']+16*s['frodo']+64
            B=2*((a*math.exp(samples*s['D']))**(1-1/s['alpha'])+delta+q/M)
            self.assertAlmostEqual(bound(s,1,frodo=True)['security_bits'],32-math.log2(B),places=9)

if __name__=='__main__': unittest.main()
