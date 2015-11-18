STGFM
=====

Fortran implementation of formulas for symplectic tracking based on generating
functions applied to multipoles, developed at Helmholtz-Zentrum Berlin 
by Malte Titze, Johannes Bahrdt and Godehard Wuestefeld.

* Windows
* Linux
* Mac OS X

Documentation
-------------
The code is not optimized for numeric applications or meant as a starting point for
higher-order implementations (because it will get very clumsy this way).
It should only test the principles. Especially the Fourier decomposition for the higher-order
generating functions and the cylinder coordinates for small offsets
produce systematic errors. The basic principles however do not require such restrictions and do not have
such errors.  Further documentation has to be done, contact the authors in case of problems.

Pre-requisites
--------------
Working fortran compiler.

Downloading
-----------
### Latest (HEAD) ###

Via git:

    git clone https://github.com/mtitze/stgfm

Installation
------------

    gfortran STGFM.for -o STGFM.exe

License
-------
*
Copyright (c) 2014, 2015 by Johannes Bahrdt (johannes.bahrdt@helmholtz-berlin.de)
and Malte Titze (malte.titze@cern.ch)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*

*The URAD subroutine was created by Michael Scheer and
is distributed under the terms of the licence found in urad.f*

Contributors
------------
### [GitHub Contributor Graph](https://github.com/mtitze/stgfm/graphs/contributors) ###

### Current Maintainers: ###
* Malte Titze

### Original Authors: ###
* Johannes Bahrdt
* Malte Titze

### Contributors not included in github history ###
* Johannes Bahrdt
* Michael Scheer


