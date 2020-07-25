import numpy as np
import fathon
from fathon import fathonUtils as fu
import math
import os

# FUNCTIONALITY TESTS
# -------------------
# - Results of the Ihlen's "Introduction to Multifractal Detrended Fluctuation
#   Analysis in Matlab" paper are reproduced
#   (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3366552/)
# - MATLAB code to reproduce the paper results can be found at:
#   https://it.mathworks.com/matlabcentral/fileexchange/38262-multifractal-detrended-fluctuation-analyses

TESTS_PATH = os.path.dirname(os.path.abspath(__file__))

#####
# Time series to test functionality
# (extracted from Ihlen's MATLAB package "Multifractal detrended fluctuation analyses")
#####
# white noise
wn = np.loadtxt(os.path.join(TESTS_PATH, 'mat/whitenoise.txt'))
# monofractal time series
mn = np.loadtxt(os.path.join(TESTS_PATH, 'mat/monofractal.txt'))
# multifractal time series
mf = np.loadtxt(os.path.join(TESTS_PATH, 'mat/multifractal.txt'))
#####

# defining parameters and scales as in Ihlen's paper
nMin = 16
nMax = 1024
nScales = 19
exponents = np.linspace(np.log2(nMin), np.log2(nMax), nScales)
scales = np.round(np.power(2.0, exponents)).astype(int)
q_list = [-5, -3, -1, 1, 3, 5]

def get_idxs(vec, scales):
    idxs = []
    for s in scales:
        idxs.append(np.where(vec==s)[0][0])
    return idxs

#####
# Functionality test 1
# It tests DFA for the white noise time series
# H should be 0.45
#####
def test_mat_dfa_wn():
    w_dfa = fathon.DFA(fu.toAggregated(wn))
    #n_w, F_w = w_dfa.computeFlucVec(nMin=nMin, nMax=nMax, revSeg=False, nStep=1, polOrd=1)
    n_w, F_w = w_dfa.computeFlucVec(scales, revSeg=False, polOrd=1)
    idxs = get_idxs(n_w, scales)
    n_w = n_w[idxs]
    F_w = F_w[idxs]
    H, _ = np.polyfit(np.log2(n_w), np.log2(F_w), 1)
    assert math.isclose(H, 0.45, rel_tol=1e-2, abs_tol=0)

#####
# Functionality test 2
# It tests DFA for the monofractal time series
# H should be 0.72
#####
def test_mat_dfa_mn():
    mn_dfa = fathon.DFA(fu.toAggregated(mn))
    #n_mn, F_mn = mn_dfa.computeFlucVec(nMin=nMin, nMax=nMax, revSeg=False, nStep=1, polOrd=1)
    n_mn, F_mn = mn_dfa.computeFlucVec(scales, revSeg=False, polOrd=1)
    idxs = get_idxs(n_mn, scales)
    n_mn = n_mn[idxs]
    F_mn = F_mn[idxs]
    H, _ = np.polyfit(np.log2(n_mn), np.log2(F_mn), 1)
    assert math.isclose(H, 0.72, rel_tol=1e-2, abs_tol=0)

#####
# Functionality test 3
# It tests DFA for the multifractal time series
# H should be 0.77
#####
def test_mat_dfa_mf():
    mf_dfa = fathon.DFA(fu.toAggregated(mf))
    #n_mf, F_mf = mf_dfa.computeFlucVec(nMin=nMin, nMax=nMax, revSeg=False, nStep=1, polOrd=1)
    n_mf, F_mf = mf_dfa.computeFlucVec(scales, revSeg=False, polOrd=1)
    idxs = get_idxs(n_mf, scales)
    n_mf = n_mf[idxs]
    F_mf = F_mf[idxs]
    H, _ = np.polyfit(np.log2(n_mf), np.log2(F_mf), 1)
    assert math.isclose(H, 0.77, rel_tol=1e-2, abs_tol=0)

#####
# Functionality test 4
# It tests MFDFA for the white noise time series
# Hq should be [0.4583, 0.4555, 0.4546, 0.4515, 0.4445, 0.4340]
#####
def test_mat_mfdfa_wn():
    w_mfdfa = fathon.MFDFA(fu.toAggregated(wn))
    #n_w, F_w = w_mfdfa.computeFlucVec(nMin=nMin, qList=q_list, nMax=nMax, revSeg=False, nStep=1, polOrd=1)
    n_w, F_w = w_mfdfa.computeFlucVec(scales, qList=q_list, revSeg=False, polOrd=1)
    idxs = get_idxs(n_w, scales)
    n_w = n_w[idxs]
    F_w_vec = np.zeros((len(q_list), len(idxs)))
    for i in range(len(q_list)):
        F_w_vec[i] = F_w[i, idxs]
    Hq = []
    for i in range(len(q_list)):
        Hq.append(np.polyfit(np.log2(n_w), np.log2(F_w_vec[i]), 1)[0])
    np.testing.assert_allclose(Hq, [0.4583, 0.4555, 0.4546, 0.4515, 0.4445, 0.4340], rtol=1e-4, atol=0)

#####
# Functionality test 5
# It tests MFDFA for the monofractal time series
# Hq should be [0.7542, 0.7392, 0.7301, 0.7240, 0.7149, 0.7023]
#####
def test_mat_mfdfa_mn():
    mn_mfdfa = fathon.MFDFA(fu.toAggregated(mn))
    #n_mn, F_mn = mn_mfdfa.computeFlucVec(nMin=nMin, qList=q_list, nMax=nMax, revSeg=False, nStep=1, polOrd=1)
    n_mn, F_mn = mn_mfdfa.computeFlucVec(scales, qList=q_list, revSeg=False, polOrd=1)
    idxs = get_idxs(n_mn, scales)
    n_mn = n_mn[idxs]
    F_mn_vec = np.zeros((len(q_list), len(idxs)))
    for i in range(len(q_list)):
        F_mn_vec[i] = F_mn[i, idxs]
    Hq = []
    for i in range(len(q_list)):
        Hq.append(np.polyfit(np.log2(n_mn), np.log2(F_mn_vec[i]), 1)[0])
    np.testing.assert_allclose(Hq, [0.7542, 0.7392, 0.7301, 0.7240, 0.7149, 0.7023], rtol=1e-4, atol=0)

#####
# Functionality test 6
# It tests MFDFA for the multifractal time series
# Hq should be [1.4477, 1.3064, 1.0823, 0.8846, 0.6606, 0.5174]
#####
def test_mat_mfdfa_mf():
    mf_mfdfa = fathon.MFDFA(fu.toAggregated(mf))
    #n_mf, F_mf = mf_mfdfa.computeFlucVec(nMin=nMin, qList=q_list, nMax=nMax, revSeg=False, nStep=1, polOrd=1)
    n_mf, F_mf = mf_mfdfa.computeFlucVec(scales, qList=q_list, revSeg=False, polOrd=1)
    idxs = get_idxs(n_mf, scales)
    n_mf = n_mf[idxs]
    F_mf_vec = np.zeros((len(q_list), len(idxs)))
    for i in range(len(q_list)):
        F_mf_vec[i] = F_mf[i, idxs]
    Hq = []
    for i in range(len(q_list)):
        Hq.append(np.polyfit(np.log2(n_mf), np.log2(F_mf_vec[i]), 1)[0])
    np.testing.assert_allclose(Hq, [1.4477, 1.3064, 1.0823, 0.8846, 0.6606, 0.5174], rtol=1e-4, atol=0)

#---------------------------------------------------------------------------------------------------------

# REGRESSION TESTS
# ----------------
# - the element of the array tested in tests 2-5,7 has no
#   particular meaning, it was just chosen randomly.
# - the numbers to which results are tested against have
#   been produced by fathon and have many significant figures
#   due to the fact that the precision of the testing machine
#   is not known.

#####
# Time series for regression tests
#####
# co2 residuals after yearly cycle removal
ts1 = np.loadtxt(os.path.join(TESTS_PATH, 'ts1.txt'))
ts1 = fu.toAggregated(ts1)
# other co2 residuals after yearly cycle removal
ts2 = np.loadtxt(os.path.join(TESTS_PATH, 'ts2.txt'))
ts2 = fu.toAggregated(ts2)
# multifractal data
ts3 = np.loadtxt(os.path.join(TESTS_PATH, 'ts3.txt'))
ts3 = fu.toAggregated(ts3)
#####

#####
# Regression test 1
# It tests if the Hurst exponent of `ts1` is correct
#####
def test_dfa_lin_step():
    pydfa = fathon.DFA(ts1)
    winSizes = fu.linRangeByStep(10, 200)
    #n1, F1 = pydfa.computeFlucVec(10, nMax=200, revSeg=True)
    n1, F1 = pydfa.computeFlucVec(winSizes, revSeg=True)
    H1, H_int1 = pydfa.fitFlucVec()

    assert math.isclose(H1, 0.7982194289592676)

#####
# Regression test 1.1
# It tests if the Hurst exponent of `ts1` is correct
#####
def test_dfa_lin_count():
    pydfa = fathon.DFA(ts1)
    winSizes = fu.linRangeByCount(10, 200)
    n1, F1 = pydfa.computeFlucVec(winSizes, revSeg=True)
    H1, H_int1 = pydfa.fitFlucVec()

    assert math.isclose(H1, 0.7982194289592676)

#####
# Regression test 1.2
# It tests if the Hurst exponent of `ts1` is correct
#####
def test_dfa_pow_step():
    pydfa = fathon.DFA(ts1)
    winSizes = fu.powRangeByStep(2, 7, step=2)
    n1, F1 = pydfa.computeFlucVec(winSizes, revSeg=True)
    H1, H_int1 = pydfa.fitFlucVec()
    
    assert math.isclose(H1, 0.7316763359204684)

#####
# Regression test 1.3
# It tests if the Hurst exponent of `ts1` is correct
#####
def test_dfa_pow_count():
    pydfa = fathon.DFA(ts1)
    winSizes = fu.powRangeByCount(1, 5, count=4, base=3)
    n1, F1 = pydfa.computeFlucVec(winSizes, revSeg=True)
    H1, H_int1 = pydfa.fitFlucVec()
    
    assert math.isclose(H1, 0.8465301940909294)

#####
# Regression test 1.3
# It tests if the Hurst exponent of `ts1` is correct
# using different log bases
#####
def test_dfa_pow_count():
    pydfa = fathon.DFA(ts1)
    winSizes = fu.powRangeByCount(1, 5, count=4, base=3)
    n1, F1 = pydfa.computeFlucVec(winSizes, revSeg=True)
    H1, H_int1 = pydfa.fitFlucVec()
    H2, H_int2 = pydfa.fitFlucVec(logBase=2)
    H3, H_int3 = pydfa.fitFlucVec(logBase=10)
    H4, H_int4 = pydfa.fitFlucVec(logBase=1.4)
    
    assert (math.isclose(H1, H2) and math.isclose(H1, H3) and math.isclose(H1, H4)
            and math.isclose(H2, H3) and math.isclose(H2, H4) and math.isclose(H3, H4))

#####
# Regression test 2
# It tests if the generalised Hurst exponent of `ts3`
# for q = -1 is correct
#####
def test_mfdfa():
    pymfdfa = fathon.MFDFA(ts3)
    qs = np.arange(-3, 3, 1)
    winSizes = fu.linRangeByStep(10, 200)
    #n2, F2 = pymfdfa.computeFlucVec(10, qs, nMax=200, revSeg=True)
    n2, F2 = pymfdfa.computeFlucVec(winSizes, qs, revSeg=True)
    H2, H_int2 = pymfdfa.fitFlucVec()

    assert math.isclose(H2[2], 1.1956312585360254)

#####
# Regression test 3
# It tests if the fifth element of the multifractal
# spectrum of `ts3` is correct
#####
def test_multifractal_spectrum():
    pymfdfa = fathon.MFDFA(ts3)
    qs = np.arange(-3, 3, 1)
    winSizes = fu.linRangeByStep(10, 200)
    #n2, F2 = pymfdfa.computeFlucVec(10, qs, nMax=200, revSeg=True)
    n2, F2 = pymfdfa.computeFlucVec(winSizes, qs, revSeg=True)
    H2, H_int2 = pymfdfa.fitFlucVec()
    a2, m2 = pymfdfa.computeMultifractalSpectrum()

    assert math.isclose(m2[4], 0.8445200259231695)

#####
# Regression test 4
# It tests if the fiftythird element of the cross-correlation
# coefficient between `ts1` and `ts2` is correct
#####
def test_dcca():
    pydcca = fathon.DCCA(ts1, ts2)
    winSizes = fu.linRangeByStep(10, 200, step=2)
    n3, rho = pydcca.computeRho(winSizes)

    assert math.isclose(rho[53], 0.48503322468665233)

#####
# Regression test 5
# It tests if the fiftythird elements of the lower and upper
# confidence levels of the cross-correlation coefficient
# for series of the same length of `ts1` and `ts2` are correct
#####
def test_rho_thresholds():
    np.random.seed(42)
    pydcca = fathon.DCCA()
    winSizes = fu.linRangeByStep(10, 200, step=2)
    n4, int1, int2 = pydcca.rhoThresholds(len(ts1), winSizes, 10, 0.95)

    assert math.isclose(int1[53], 0.03131478865331007) and math.isclose(int2[53], -0.05672796198121624)

#####
# Regression test 6
# It tests if the Hurst exponent between `ts1`
# and `ts2` is correct
#####
def test_dcca():
    pydcca = fathon.DCCA(ts1, ts2)
    winSizes = fu.linRangeByStep(10, 200, step=2)
    n, F = pydcca.computeFlucVec(winSizes, polOrd=2)
    H, H_int = pydcca.fitFlucVec()
    
    assert math.isclose(H, 0.9604004237165071)

#####
# Regression test 7
# It tests if the Hurst exponent between `ts1`
# and `ts2` is correct
#####
def test_dcca():
    pyht = fathon.HT(ts2)
    scales = [10, 100]
    ht = pyht.computeHt(scales, mfdfaPolOrd=3)
    
    #assert math.isclose(H, 0.9604004237165071)
