======
uvplot
======
A simple package to make nice plots of deprojected interferometric visibilities, often called **uvplot** s.


.. image:: https://travis-ci.org/mtazzari/uvplot.svg?branch=master
    :target: https://travis-ci.org/mtazzari/uvplot

.. image:: https://img.shields.io/github/release/mtazzari/uvplot/all.svg
    :target: https://github.com/mtazzari/uvplot/releases
    
.. image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
    :target: https://www.gnu.org/licenses/lgpl-3.0


|

The current version works on `Python` >2.7 and >3.6 and implements the basic plotting functionality. More features are to come in the future.

If you are interested, or have feature requests, or encounter issues, consider creating an Issue or writing me an `email  <mtazzari@ast.cam.ac.uk>`_.

This is an example plot:

.. image:: static/uvplot.png
   :width: 60 %
   :alt: example uv plot
   :align: center

created with:

.. code-block:: py

    import numpy as np
    from uvplot import UVTable, arcsec

    wle = 0.88e-3         # Observing wavelength         [m]

    dRA = 0.3 * arcsec    # Delta Right Ascension offset [rad]
    dDec = 0.07 * arcsec  # Delta Declination     offset [rad]
    inc = np.radians(73.) # Inclination    [rad]
    PA = np.radians(59)   # Position Angle [rad]

    uvbin_size = 30e3     # uv-distance bin [wle]

    uv = UVTable(filename='uvtable.txt', wle=wle)
    uv.apply_phase(dRA, dDec)
    uv.deproject(inc, PA)

    uv_mod = UVTable(filename='uvtable_mod.txt', wle=wle)
    uv_mod.apply_phase(dRA=dRA, dDec=dDec)
    uv_mod.deproject(inc=inc, PA=PA)

    axes = uv.plot(label='Data', uvbin_size=uvbin_size)
    uv_mod.plot(label='Model', uvbin_size=uvbin_size, axes=axes, yerr=False, linestyle='-', color='r')

    axes[0].figure.savefig("uvplot.png")


License
-------
**uvplot** is free software licensed under the LGPLv3 License. For more details see the LICENSE.

© Copyright 2017 Marco Tazzari.
