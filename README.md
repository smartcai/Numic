Numic 
======

Numic is short for NUmerical Matrix In C.

The design philosophy is to build up a clean but not huge matrix computations
package, so that any project in the C language can incorporate it as part.

The major components include operations on matrix and vectors such as BLAS,
Householder and Givens transformations, QR, SVD and bidiagonal decompositions,
PCA, and some optimizations algorithms based on them.

./src
---

- libmatrix.h - The interface.
- libmatrix.c - The major C code.
- example - Some test code.
- prototype - Experimental ideas in Octave/Matlab.


./doc
-----

The document.

Major Reference
---------------

Matrix Computations by Gene H. Golub, Charles F. Van Loan

Numerical Linear Algebra by Lloyd N. Trefethen, David Bau III
