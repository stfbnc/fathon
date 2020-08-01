powRangeByStep
==============

.. currentmodule:: fathon.fathonUtils

.. autofunction:: powRangeByStep

Usage examples
^^^^^^^^^^^^^^

.. code-block:: python

   from fathon import fathonUtils as fu

   #elements given by `base` raised to linearly separated exponents
   a = fu.powRangeByStep(10, 1000, step=2, base=3)
