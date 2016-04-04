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

from __future__ import division
from math import pi

__all__ = ['sphericity', 'aspect_ratio', 'circularity', 'A_cylinder',
'V_cylinder', 'A_hollow_cylinder', 'V_hollow_cylinder',
'A_multiple_hole_cylinder', 'V_multiple_hole_cylinder']


def sphericity(A, V):
    r'''Returns the sphericity of a particle of surface area `A` and volume
    `V`.

    .. math::
        \Psi = \frac{\text{A of sphere with } V_p }  {{A}_p}
        = \frac{\pi^{\frac{1}{3}}(6V_p)^{\frac{2}{3}}}{A_p}

    Parameters
    ----------
    A : float
        Surface area of particle, [m^2]
    V : float
        Volume of particle, [m^3]

    Returns
    -------
    Psi : float
        Sphericity [-]

    Examples
    --------
    >>> sphericity(10., 2.)
    0.767663317071005

    References
    ----------
    .. [1] Rhodes, Martin J., ed. Introduction to Particle Technology. 2E.
    Chichester, England ; Hoboken, NJ: Wiley, 2008.
    '''
    Psi = pi**(1/3.)*(6*V)**(2/3.)/A
    return Psi


def aspect_ratio(Dmin, Dmax):
    r'''Returns the aspect ratio of a shape with minimum and maximum dimension,
    `Dmin` and `Dmax`.

    .. math::
        A_R = \frac{D_{min}}{D_{max}}

    Parameters
    ----------
    Dmin : float
        Minimum dimension, [m]
    Dmax : float
        Maximum dimension, [m]

    Returns
    -------
    a_r : float
        Aspect ratio [-]

    Examples
    --------
    >>> aspect_ratio(.2, 2)
    0.1
    '''
    a_r = Dmin/Dmax
    return a_r


def circularity(A, P):
    r'''Returns the circularity of a shape with area `A` and perimeter `P`.

    .. math::
        f_{circ} = \frac {4 \pi A} {P^2}

    Parameters
    ----------
    A : float
        Area of the shape, [m^2]
    P : float
        Perimeter of the shape, [m]

    Returns
    -------
    f_circ : float
        Circularity of the shape [-]

    Examples
    --------
    >>> circularity(1.5, .1)
    1884.9555921538756
    '''
    f_circ = 4*pi*A/P**2
    return f_circ


def A_cylinder(D, L):
    r'''Returns the surface area of a cylinder.

    .. math::
        A = \pi D L + 2\cdot \frac{\pi D^2}{4}

    Parameters
    ----------
    D : float
        Diameter of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    A : float
        Surface area [m]

    Examples
    --------
    >>> A_cylinder(0.01, .1)
    0.0032986722862692833
    '''
    cap = pi*D**2/4*2
    side = pi*D*L
    A = cap + side
    return A


def V_cylinder(D, L):
    r'''Returns the volume of a cylinder.

    .. math::
        V = \frac{\pi D^2}{4}L

    Parameters
    ----------
    D : float
        Diameter of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    >>> V_cylinder(0.01, .1)
    7.853981633974484e-06
    '''
    V = pi*D**2/4*L
    return V


def A_hollow_cylinder(Di, Do, L):
    r'''Returns the surface area of a hollow cylinder.

    .. math::
        A = \pi D_o L + \pi D_i L  + 2\cdot \frac{\pi D_o^2}{4}
        - 2\cdot \frac{\pi D_i^2}{4}

    Parameters
    ----------
    Di : float
        Diameter of the hollow in the cylinder, [m]
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    A : float
        Surface area [m]

    Examples
    --------
    >>> A_hollow_cylinder(0.005, 0.01, 0.1)
    0.004830198704894308
    '''
    side_o = pi*Do*L
    side_i = pi*Di*L
    cap_circle = pi*Do**2/4*2
    cap_removed = pi*Di**2/4*2
    A = side_o + side_i + cap_circle - cap_removed
    return A


def V_hollow_cylinder(Di, Do, L):
    r'''Returns the volume of a hollow cylinder.

    .. math::
        V = \frac{\pi D_o^2}{4}L - L\frac{\pi D_i^2}{4}

    Parameters
    ----------
    Di : float
        Diameter of the hollow in the cylinder, [m]
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    >>> V_hollow_cylinder(0.005, 0.01, 0.1)
    5.890486225480862e-06
    '''
    V = pi*Do**2/4*L - pi*Di**2/4*L
    return V


def A_multiple_hole_cylinder(Do, L, holes):
    r'''Returns the surface area of a cylinder with multiple holes.
    Calculation will naively return a negative value or other impossible
    result if the number of cylinders added is physically impossible.
    Holes may be of different shapes, but must be perpendicular to the
    axis of the cylinder.

    .. math::
        A = \pi D_o L + 2\cdot \frac{\pi D_o^2}{4} +
        \sum_{i}^n \left( \pi D_i L  - 2\cdot \frac{\pi D_i^2}{4}\right)

    Parameters
    ----------
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]
    holes : list
        List of tuples containing (diameter, count) pairs of descriptions for
        each of the holes sizes.

    Returns
    -------
    A : float
        Surface area [m]

    Examples
    --------
    >>> A_multiple_hole_cylinder(0.01, 0.1, [(0.005, 1)])
    0.004830198704894308
    '''
    side_o = pi*Do*L
    cap_circle = pi*Do**2/4*2
    A = cap_circle + side_o
    for Di, n in holes:
        side_i = pi*Di*L
        cap_removed = pi*Di**2/4*2
        A = A + side_i*n - cap_removed*n
    return A


def V_multiple_hole_cylinder(Do, L, holes):
    r'''Returns the solid volume of a cylinder with multiple cylindrical holes.
    Calculation will naively return a negative value or other impossible
    result if the number of cylinders added is physically impossible.

    .. math::
        V = \frac{\pi D_o^2}{4}L - L\frac{\pi D_i^2}{4}

    Parameters
    ----------
    Do : float
        Diameter of the exterior of the cylinder, [m]
    L : float
        Length of the cylinder, [m]
    holes : list
        List of tuples containing (diameter, count) pairs of descriptions for
        each of the holes sizes.

    Returns
    -------
    V : float
        Volume [m^3]

    Examples
    --------
    >>> V_multiple_hole_cylinder(0.01, 0.1, [(0.005, 1)])
    5.890486225480862e-06
    '''
    V = pi*Do**2/4*L
    for Di, n in holes:
        V -= pi*Di**2/4*L*n
    return V

