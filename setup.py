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

from distutils.core import setup

classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Intended Audience :: Education',
    'Intended Audience :: Manufacturing',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Natural Language :: English',
    'Operating System :: MacOS',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: BSD',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: Implementation :: CPython',
    'Topic :: Education',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Topic :: Scientific/Engineering :: Physics',
]

long_description = '''fpi is open-source software for engineers and technicians
working in the fields of chemical or mechanical engineering. It includes
modules for various fluid-particle interaction calculations. fpi runs on
all operating systems which support Python, is quick to install, and is free
of charge. fpi is designed to be easy to use while still providing powerful
functionality. If you need to perform some fpi calculations, give
fpi a try.'''


setup(
  name = 'fpi',
  packages = ['fpi'],
  license='GPL3',
  version = '0.1.1',
  description = 'Fluid-particle Interaction component of Chemical Engineering Design Library (ChEDL)',
  author = 'Caleb Bell',
  long_description = long_description,
  platforms=["Windows", "Linux", "Mac OS", "Unix"],
  author_email = 'Caleb.Andrew.Bell@gmail.com',
  url = 'https://github.com/CalebBell/fpi',
  download_url = 'https://github.com/CalebBell/fpi/tarball/0.1.1',
  keywords = ['chemical engineering', 'fluid mechanics', 'mechanical engineering'],
  classifiers = classifiers,
)
