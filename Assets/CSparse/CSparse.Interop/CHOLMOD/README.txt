Contents
--------

This archive contains shared libraries for SuiteSparse (version 5.3.0).

Included projects:
    libamd        - AMD, CAMD, CCOLAMD, COLAMD
    libcholmod    - CHOLMOD (with NGPL option)
    libmetis      - METIS

Dependencies:
    libamd        >
    libcholmod    > libamd, libmetis
    libmetis      >


Platform Toolset: Visual Studio 2017 (v141)
   Configuration: Release

See https://github.com/wo80/vs-suitesparse/ for additional information.
See https://github.com/stNoh/vs-suitesparse/tree/LGPL for minor changes on LGPL.
Date: 2019-08-20


License
-------

See LICENSE.txt.

Original LICENSE.txt in SuiteSparse-5.3.0 included several licenses
which are not actually used to compile binaries in this directory,
so they were excluded from LICENSE.txt now.
