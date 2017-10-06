#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

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
