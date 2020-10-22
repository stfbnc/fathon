DCCA
====

.. currentmodule:: fathon

.. autoclass:: DCCA
   :show-inheritance:

   .. automethod:: computeFlucVec
   .. automethod:: computeRho
   .. automethod:: fitFlucVec
   .. automethod:: multiFitFlucVec
   .. automethod:: rhoThresholds
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

   #initialize non-empty dcca object
   pydcca = fathon.DCCA(a, b)
   #compute fluctuation function and Hurst exponent
   wins = fu.linRangeByStep(20, 100, step=50)
   n, F = pydcca.computeFlucVec(wins, polOrd=1)
   H, H_intercept = pydcca.fitFlucVec()

   #compute Hurst exponent in different ranges
   limits_list = np.array([[20,120], [220,870]], dtype=int)
   list_H, list_H_intercept = pydcca.multiFitFlucVec(limits_list)

   #compute rho index
   wins = fu.linRangeByStep(20, 100, step=50)
   n, rho = pydcca.computeRho(wins, polOrd=1)

   #initialize empty dcca object
   pythresh = fathon.DCCA()
   #compute confidence levels
   wins = fu.linRangeByStep(4, 100)
   n, cInt1, cInt2 = pythresh.rhoThresholds(300, wins, 100, 0.95, polOrd=1)

