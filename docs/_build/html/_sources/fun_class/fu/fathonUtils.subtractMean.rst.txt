subtractMean
============

.. currentmodule:: fathon.fathonUtils

.. autofunction:: subtractMean

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   from fathon import fathonUtils as fu

   #time series
   a = np.random.randn(10000)

   #zero-mean time series
   a = fu.subtractMean(a)
