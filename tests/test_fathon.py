import numpy as np
import fathon
import math

#####
# Time series
#####
ts1 = np.loadtxt('./tests/ts1.txt')
ts1 = fathon.toAggregated(ts1)
ts2 = np.loadtxt('./tests/ts2.txt')
ts2 = fathon.toAggregated(ts2)
ts3 = np.loadtxt('./tests/ts3.txt')
ts3 = fathon.toAggregated(ts3)

#####
# Test 1
#####
def test_dfa():
    pydfa = fathon.DFA(ts1)
    n1, F1 = pydfa.computeFlucVec(10, nMax=200, revSeg=True)
    H1, H_int1 = pydfa.fitFlucVec()

    assert math.isclose(H1, 0.7982194289592676)

#####
# Test 2
#####
def test_mfdfa():
    pymfdfa = fathon.MFDFA(ts3)
    qs = np.arange(-3, 3, 1)
    n2, F2 = pymfdfa.computeFlucVec(10, qs, nMax=200, revSeg=True)
    H2, H_int2 = pymfdfa.fitFlucVec()

    assert math.isclose(H2[2], 1.1956312585360254)

#####
# Test 3
#####
def test_multifractal_spectrum():
    pymfdfa = fathon.MFDFA(ts3)
    qs = np.arange(-3, 3, 1)
    n2, F2 = pymfdfa.computeFlucVec(10, qs, nMax=200, revSeg=True)
    H2, H_int2 = pymfdfa.fitFlucVec()
    a2, m2 = pymfdfa.computeMultifractalSpectrum()

    assert math.isclose(m2[4], 0.8445200259231695)

#####
# Test 4
#####
def test_dcca():
    pydcca = fathon.DCCA(ts1, ts2)
    n3, rho = pydcca.computeRho(10, nMax=200, nStep=2)

    assert math.isclose(rho[53], 0.48503322468665233)

#####
# Test 5
#####
def test_rho_thresholds():
    np.random.seed(42)
    pydcca = fathon.DCCA()
    n4, int1, int2 = pydcca.rhoThresholds(len(ts1), 10, 200, 10, 0.95, nStep=2)

    assert math.isclose(int1[53], 0.03131478865331007) and math.isclose(int2[53], -0.05672796198121624)
