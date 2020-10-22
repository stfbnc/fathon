HT
==

.. currentmodule:: fathon

.. autoclass:: HT
   :show-inheritance:

   .. automethod:: computeHt
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

   #initialize ht object
   pyht = fathon.HT(a)
   #compute time-dependent Hurst exponent
   ht = pyht.computeHt([100, 200, 1000], mfdfaPolOrd=1, polOrd=1)

