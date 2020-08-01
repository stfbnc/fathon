toAggregated
============

.. currentmodule:: fathon.fathonUtils

.. autofunction:: toAggregated

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   from fathon import fathonUtils as fu

   #time series
   a = np.random.randn(10000)

   #zero-mean cumulative sum
   a = fu.toAggregated(a)
