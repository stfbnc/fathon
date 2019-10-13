import numpy as np
import fathon
import math

#####
ts1 = np.loadtxt('ts1.txt')
ts2 = np.loadtxt('ts2.txt')
ts3 = np.loadtxt('ts3.txt')
ts1 = fathon.toAggregated(ts1)
ts2 = fathon.toAggregated(ts2)
ts3 = fathon.toAggregated(ts3)
#####
# Test 1
#####
pydfa = fathon.DFA(ts1)
n1, F1 = pydfa.computeFlucVec(10, nMax=200, revSeg=True)
H1, H_int1 = pydfa.fitFlucVec()

print('')
if(math.isclose(H1, 0.7982194289592676)):
    print('Test 1 passed')
else:
    print('Test 1 not passed')
print('')
#####
# Test 2
#####
pymfdfa = fathon.MFDFA(ts3)
qs = np.arange(-3, 3, 1)
n2, F2 = pymfdfa.computeFlucVec(10, qs, nMax=200, revSeg=True)
H2, H_int2 = pymfdfa.fitFlucVec()

print('')
if(math.isclose(H2[2], 1.1956312585360254)):
    print('Test 2 passed')
else:
    print('Test 2 not passed')
print('')
#####
# Test 3
#####
a2, m2 = pymfdfa.computeMultifractalSpectrum()

print('')
if(math.isclose(m2[4], 0.8445200259231695)):
    print('Test 3 passed')
else:
    print('Test 3 not passed')
print('')
#####
# Test 4
#####
pydcca = fathon.DCCA(ts1, ts2)
n3, rho = pydcca.computeRho(10, nMax=200, nStep=2)

print('')
if(math.isclose(rho[53], 0.48503322468665233)):
    print('Test 4 passed')
else:
    print('Test 4 not passed')
print('')
#####
# Test 5
#####
np.random.seed(42)
n4, int1, int2 = pydcca.rhoThresholds(len(ts1), 10, 200, 10, 0.95, nStep=2)

print('')
if(math.isclose(int1[53], 0.03131478865331007) and math.isclose(int2[53], -0.05672796198121624)):
    print('Test 5 passed')
else:
    print('Test 5 not passed')
print('')
#####
