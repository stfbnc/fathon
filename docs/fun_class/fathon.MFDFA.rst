MFDFA
=====

.. currentmodule:: fathon

.. autoclass:: MFDFA
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

   #zero-mean cumulative sum
   a = fu.toAggregated(a)

   #initialize mfdfa object
   pymfdfa = fathon.MFDFA(a)
   #compute fluctuation function and generalized Hurst exponents
   wins = fu.linRangeByStep(10, 2000)
   n, F = pymfdfa.computeFlucVec(wins, np.arange(-3, 4, 0.1), revSeg=True, polOrd=1)
   list_H, list_H_intercept = pymfdfa.fitFlucVec()

   #compute mass exponents
   tau = pymfdfa.computeMassExponents()

   #compute multifractal spectrum
   alpha, mfSpect = pymfdfa.computeMultifractalSpectrum()

