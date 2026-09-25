.. highlight:: shell

============
Installation
============


Stable release
--------------

To install the SEAMM Forcefield Utilities, run this command in your terminal:

.. code-block:: console

    $ pip install seamm-ff-util

This installs every dependency, including the compiled ones (RDKit and Open Babel,
via molsystem), from their PyPI wheels. Use `pip`_ or uv; a conda-forge package also
exists but lags the PyPI release, and the two should not be mixed in one environment.

.. _pip: https://pip.pypa.io

From sources
------------

The sources for the seamm_ff_util can be downloaded
from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/molssi-seamm/seamm_ff_util

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/molssi-seamm/seamm_ff_util/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/molssi-seamm/seamm_ff_util
.. _tarball: https://github.com/molssi-seamm/seamm_ff_util/tarball/master
