MFDFA
=====

.. currentmodule:: fathon

.. autoclass:: MFDFA
   :show-inheritance:

   .. automethod:: computeFlucVec
   .. automethod:: computeMassExponents
   .. automethod:: computeMultifractalSpectrum
   .. automethod:: fitFlucVec

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   import fathon

   #time series
   a = np.random.randn(10000)

   #zero-mean cumulative sum
   a = fathon.toAggregated(a)

   #initialize mfdfa object
   pymfdfa = fathon.MFDFA(a)
   #compute fluctuation function and generalized Hurst exponents
   n, F = pymfdfa.computeFlucVec(10, np.arange(-3, 4, 0.1), nMax=2000, revSeg=True, nStep=1, polOrd=1)
   list_H, list_H_intercept = pymfdfa.fitFlucVec()

   #compute mass exponents
   tau = pymfdfa.computeMassExponents()

   #compute multifractal spectrum
   alpha, mfSpect = pymfdfa.computeMultifractalSpectrum()

