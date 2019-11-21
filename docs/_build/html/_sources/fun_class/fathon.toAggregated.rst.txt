toAggregated
============

.. currentmodule:: fathon

.. autofunction:: toAggregated

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   import fathon

   #time series
   a = np.random.randn(10000)

   #zero-mean cumulative sum
   a = fathon.toAggregated(a)
