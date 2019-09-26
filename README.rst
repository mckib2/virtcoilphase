Virtual Reference Coil
======================

Python implementation of the algorithm described in [1]_ for finding
optimal phase map estimates using a constructed virtual coil.

Installation
============

Should be as easy as a pip install:

.. code-block:: bash

    pip install virtcoilphase

Usage
=====

Examples can be called from the command line, for example:

.. code-block:: bash

    python -m virtcoilphase.examples.simple

Given a number of coil images, you can find the absolute phase
estimate like this:

.. code-block:: python

    from virtcoilphase import virtcoilphase
    phase = virtcoilphase(x, coil_axis=-1)

References
==========
.. [1] Parker, Dennis L., et al. "Phase reconstruction from
       multiple coil data using a virtual reference coil."
       Magnetic resonance in medicine 72.2 (2014): 563-569.
