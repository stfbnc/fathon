import numpy as np
import matplotlib.pyplot as plt
import fathon
from fathon import fathonUtils as fu

print('This is fathon v{}'.format(fathon.__version__))

a = np.random.randn(10000)
b = np.random.randn(10000)

a = fu.toAggregated(a)
b = fu.toAggregated(b)

pydcca = fathon.DCCA(a, b)

winSizes = fu.linRangeByStep(20, 1000, step=50)
polOrd = 1

n, F = pydcca.computeFlucVec(winSizes, polOrd=polOrd)

H, H_intercept = pydcca.fitFlucVec()

plt.plot(np.log(n), np.log(F), 'ro')
plt.plot(np.log(n), H_intercept+H*np.log(n), 'k-', label='H = {:.2f}'.format(H))
plt.xlabel('ln(n)', fontsize=14)
plt.ylabel('ln(F(n))', fontsize=14)
plt.title('DCCA', fontsize=14)
plt.legend(loc=0, fontsize=14)
plt.show()

limits_list = np.array([[20,120], [220,870]], dtype=int)
list_H, list_H_intercept = pydcca.multiFitFlucVec(limits_list)

clrs = ['k', 'b', 'm', 'c', 'y']
stls = ['-', '--', '.-']
plt.plot(np.log(n), np.log(F), 'ro')
for i in range(len(list_H)):
    n_rng = np.arange(limits_list[i][0], limits_list[i][1]+1, 50)
    plt.plot(np.log(n_rng), list_H_intercept[i]+list_H[i]*np.log(n_rng),
             clrs[i%len(clrs)]+stls[(i//len(clrs))%len(stls)], label='H = {:.2f}'.format(list_H[i]))
plt.xlabel('ln(n)', fontsize=14)
plt.ylabel('ln(F(n))', fontsize=14)
plt.title('DCCA', fontsize=14)
plt.legend(loc=0, fontsize=14)
plt.show()

winSizes = fu.linRangeByStep(20, 1000, step=50)
polOrd = 1

n, rho = pydcca.computeRho(winSizes, polOrd=polOrd)

plt.plot(n, rho, 'ro')
plt.ylim(-1, 1)
plt.xlabel('n', fontsize=14)
plt.ylabel('$\\rho_{DCCA}$', fontsize=14)
plt.title('rhoDCCA', fontsize=14)
plt.show()

pythresh = fathon.DCCA()

L = 300
winSizes = fu.linRangeByStep(4, 100, step=1)
nSim = 100
confLvl = 0.95
polOrd = 1

n, cInt1, cInt2 = pythresh.rhoThresholds(L, winSizes, nSim, confLvl, polOrd=polOrd, verbose=True)

plt.plot(n, cInt1, 'r-')
plt.plot(n, cInt2, 'b-')
plt.ylim(-1, 1)
plt.xlabel('n', fontsize=14)
plt.show()
