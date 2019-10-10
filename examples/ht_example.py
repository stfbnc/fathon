import numpy as np
import matplotlib.pyplot as plt
import fathon

print('This is fathon v{}'.format(fathon.__version__))

a = np.random.randn(10000)
b = np.random.randn(10000)

a = fathon.toAggregated(a)
b = fathon.toAggregated(b)

pyht = fathon.HT(a)

scales = [100, 200, 1000]
mfdfaPolOrd = 1
polOrd = 1

ht = pyht.computeHt(scales, mfdfaPolOrd=mfdfaPolOrd, polOrd=polOrd)

plt.rc('font', size=14)
plt.figure(figsize=(10, 6))
w = 3 if len(scales) >= 3 else len(scales)
h = np.ceil(len(scales)/3)
for i, scale in enumerate(scales):
    plt.subplot(h, w, i+1)
    plt.plot(np.arange(1, len(ht[i, 0:len(a)-scale+1])+1), ht[i, 0:len(a)-scale+1],
             'r-', label='scale = {}'.format(scale))
    plt.xlabel('window number', fontsize=14)
    plt.ylabel('$h_t$', fontsize=14)
    plt.legend(loc=0, fontsize=14)
plt.subplots_adjust(hspace=0.6, wspace=0.3)
plt.show()
