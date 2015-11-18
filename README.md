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
*STGFM is free software. You may copy, distribute, and modify it under
the terms of the License contained in the file LICENSE distributed
with this package. The URAD subroutine was created by Michael Scheer and
is distributed under the licence found in urad.f*

Contributors
------------
### [GitHub Contributor Graph](https://github.com/mtitze/stgfm) ###

### Current Maintainers: ###
* Malte Titze

### Original Authors: ###
* Malte Titze
* Johannes Bahrdt

### Contributors not included in github history ###
* Johannes Bahrdt
* Michael Scheer

