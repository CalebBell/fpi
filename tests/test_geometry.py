# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''

from fpi import *
import numpy as np
from numpy.testing import assert_allclose
import pytest

close = assert_allclose


def test_geometry():
    psi = sphericity(10., 2.)
    close(psi, 0.767663317071005)

    a_r = aspect_ratio(.2, 2)
    close(a_r, 0.1)

    f_circ = circularity(1.5, .1)
    close(f_circ, 1884.9555921538756)

    A = A_cylinder(0.01, .1)
    close(A, 0.0032986722862692833)

    V = V_cylinder(0.01, .1)
    close(V, 7.853981633974484e-06)

    A = A_hollow_cylinder(0.005, 0.01, 0.1)
    close(A, 0.004830198704894308)

    V = V_hollow_cylinder(0.005, 0.01, 0.1)
    close(V, 5.890486225480862e-06)

    A =  A_multiple_hole_cylinder(0.01, 0.1, [(0.005, 1)])
    close(A, 0.004830198704894308)

    V = V_multiple_hole_cylinder(0.01, 0.1, [(0.005, 1)])
    close(V, 5.890486225480862e-06)