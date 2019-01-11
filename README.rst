======
uvplot
======
A simple package to make nice plots of deprojected interferometric visibilities, often called **uvplots**.
It can be installed inside the `NRAO CASA package <https://casa.nrao.edu/>`_ (see instructions below) and has functionalities to export visibilities from the MS Table format to ASCII.


.. image:: https://travis-ci.org/mtazzari/uvplot.svg?branch=master
    :target: https://travis-ci.org/mtazzari/uvplot

.. image:: https://img.shields.io/github/release/mtazzari/uvplot/all.svg
    :target: https://github.com/mtazzari/uvplot/releases
    
.. image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
    :target: https://www.gnu.org/licenses/lgpl-3.0

.. image:: https://zenodo.org/badge/105298533.svg
   :target: https://zenodo.org/badge/latestdoi/105298533
   
|

The current version implements the basic plotting functionality.

Features on the road map:
    - better handling of multiple spectral windows during visibilities export;
    - new functionality to import visibilities from ASCII to MS Table.

If you are interested, have feature requests, or encounter issues, consider creating an `Issue <https://github.com/mtazzari/uvplot/issues>`_ or writing me an `email  <marco.tazzari@gmail.com>`_. I am happy to have your feedback!


Attribution
-----------
If you use **uvplot** for your publication, please cite the `Zenodo reference <https://zenodo.org/badge/latestdoi/105298533>`_ ::

    @misc{uvplot_mtazzari,
      author       = {Marco Tazzari},
      title        = {mtazzari/uvplot: v0.1.1},
      month        = oct,
      year         = 2017,
      doi          = {10.5281/zenodo.1003113},
      url          = {https://doi.org/10.5281/zenodo.1003113}
    }


License
-------
**uvplot** is free software licensed under the LGPLv3 License. For more details see the LICENSE.
© Copyright 2018-2019 Marco Tazzari.

Documentation
-------------
Check out the `documentation <https://mtazzari.github.io/uvplot/>`_.

Changelog
---------
- **v0.2.4**: support for binary (.npz) uvtables; introduce COLUMNS format.
- **v0.2.3**: a dedicated `documentation <https://mtazzari.github.io/uvplot/>`_ website.
- **v0.2.2**: a new export visibilities option in UVTable.plot(), automatically mask empty uv-bins, bugfixes.
- **v0.2.0**: a new `export_uvtable` function to export visibilities from an MS to an ASCII table.


Installation
------------

**uvplot** works on `Python` >=2.7 and >=3.6 and can be installed with:

.. code-block :: bash

    pip install uvplot

To make **uvplot** available in CASA, run from the shell:

.. code-block :: bash

    casa-pip install uvplot

where `casa-pip` is a tool that can be downloaded `here <https://github.com/radio-astro-tools/casa-python>`_ .

**uvplot** has been tested on CASA versions >= 4.7.0.

Example
-------
This is an example plot:

.. image:: static/uvplot.png
   :width: 60 %
   :alt: example uv plot
   :align: center

created with uvplot:

.. code-block:: py

    import numpy as np
    from uvplot import UVTable, arcsec
    from uvplot import COLUMNS_V0       # use uvplot >= 0.2.4

    wle = 0.88e-3         # Observing wavelength         [m]

    dRA = 0.3 * arcsec    # Delta Right Ascension offset [rad]
    dDec = 0.07 * arcsec  # Delta Declination     offset [rad]
    inc = np.radians(73.) # Inclination    [rad]
    PA = np.radians(59)   # Position Angle [rad]

    uvbin_size = 30e3     # uv-distance bin [wle]

    uv = UVTable(filename='uvtable.txt', wle=wle, columns=COLUMNS_V0)
    uv.apply_phase(dRA, dDec)
    uv.deproject(inc, PA)

    uv_mod = UVTable(filename='uvtable_mod.txt', wle=wle, COLUMNS_V0)
    uv_mod.apply_phase(dRA=dRA, dDec=dDec)
    uv_mod.deproject(inc=inc, PA=PA)

    axes = uv.plot(label='Data', uvbin_size=uvbin_size)
    uv_mod.plot(label='Model', uvbin_size=uvbin_size, axes=axes, yerr=False, linestyle='-', color='r')

    axes[0].figure.savefig("uvplot.png")

From version v0.2.4 it is necessary to provide the `columns` parameter
when reading an ASCII uvtable. The `columns` parameter can be specified
either as a parameter to the `UVTable()` command, or as the **2nd** line
in the ASCII file. The available `columns` formats are:

.. code-block:: bash

    FORMAT          COLUMNS                                                 COLUMNS_LINE (copy-paste as 2nd line in the ASCII file)
    COLUMNS_V0      ['u', 'v', 'Re', 'Im', 'weights']                       '# Columns      u v Re Im weights'
    COLUMNS_V1      ['u', 'v', 'Re', 'Im', 'weights', 'freqs', 'spws']      '# Columns      u v Re Im weights freqs spws'
    COLUMNS_V2      ['u', 'v', 'V', 'weights', 'freqs', 'spws']             '# Columns      u v V weights freqs spws'

To import an ASCII uvtable with 5 columns with uvplot < 0.2.4:

.. code-block:: py

    from uvplot import UVTable
    uvt = UVTable(filename='uvtable.txt', format='ascii', columns=COLUMNS_V0)


and with uvplot >= 0.2.4:

.. code-block:: py

    from uvplot import UVTable
    from uvplot import COLUMNS_V0  # ['u', 'v', 'Re', 'Im', 'weights']
    uvt = UVTable(filename='uvtable.txt', format='ascii', columns=COLUMNS_V0)

