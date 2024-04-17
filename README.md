# The Improved Akima ('86) Method
A python-wrapped modern Fortran implementation of the Improved Akima Interpolation Method. The `akima` interpolation found in scipy or Matlab is the original algorithm from 1970.
Akima wrote a technical document in 1986 (linked below) that highlighted weaknesses in the original method and provided proposed improvements.

# Installation
Requires: Numpy and a Fortran compiler supported by f2py and OpenMP;
  also requires matplotlib for tests

Installation:

`python setup.py install`

(to change default fortran compiler you can use e.g.
 `python setup.py build config_fc --fcompiler=g95`)

Citation:
Akima, H., 1986: A method of univariate interpolation that has
    the accuracy of a third-degree polynomial. [NTIA Report 86-208. 76 pp.](https://its.ntia.gov/publications/details.aspx?pub=2231)

Legal:
Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation. For published articles, please cite the
"Citation" listed above.
THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

