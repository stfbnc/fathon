subtractMean
============

.. currentmodule:: fathon

.. autofunction:: subtractMean

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   import fathon

   #time series
   a = np.random.randn(10000)

   #zero-mean time series
   a = fathon.subtractMean(a)
