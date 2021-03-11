MFDCCA
======

.. currentmodule:: fathon

.. autoclass:: MFDCCA
   :show-inheritance:

   .. automethod:: computeFlucVec
   .. automethod:: computeMassExponents
   .. automethod:: computeMultifractalSpectrum
   .. automethod:: fitFlucVec
   .. automethod:: saveObject

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   import fathon
   from fathon import fathonUtils as fu

   #time series
   a = np.random.randn(10000)
   b = np.random.randn(10000)

   #zero-mean cumulative sum
   a = fu.toAggregated(a)
   b = fu.toAggregated(b)

   #initialize mfdfa object
   pymfdcca = fathon.MFDCCA(a, b)
   #compute fluctuation function and generalized Hurst exponents
   wins = fu.linRangeByStep(10, 2000)
   n, F = pymfdcca.computeFlucVec(wins, np.arange(-3, 4, 0.1), revSeg=True, polOrd=1)
   list_H, list_H_intercept = pymfdcca.fitFlucVec()

   #compute mass exponents
   tau = pymfdcca.computeMassExponents()

   #compute multifractal spectrum
   alpha, mfSpect = pymfdcca.computeMultifractalSpectrum()

