DFA
===

.. currentmodule:: fathon

.. autoclass:: DFA
   :show-inheritance:

   .. automethod:: computeFlucVec
   .. automethod:: fitFlucVec
   .. automethod:: multiFitFlucVec

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   import fathon

   #time series
   a = np.random.randn(10000)

   #zero-mean cumulative sum
   a = fathon.toAggregated(a)

   #initialize dfa object
   pydfa = fathon.DFA(a)
   #compute fluctuation function and Hurst exponent
   n, F = pydfa.computeFlucVec(10, nMax=2000, revSeg=True, nStep=1, polOrd=3)
   H, H_intercept = pydfa.fitFlucVec()

   #compute Hurst exponent in different ranges
   limits_list = np.array([[15,2000], [200,1000]], dtype=int)
   list_H, list_H_intercept = pydfa.multiFitFlucVec(limits_list)

