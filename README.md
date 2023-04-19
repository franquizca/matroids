# matroids
Python implementation of the class matroid

Implementation based on the ground set and independent subsets definition of matroid stated in "Matroid theory" by James G. Oxley.

Dependencies: sympy, rustworkx, numpy
Optional dependencies: graphviz (for the drawLattice method of the AbstractMatroid class)

The main file for the program is matroids.py. It implements the class AbstractMatroid and contains three functions to initialize matroids:

  unifromMatroid(n,r) which returns the matroid where any subset of the ground set {0,...,n-1} with at most r elements is independent.
  starMatroid(n) which returns the matroid associated to the line arrangement of n lines passing through a point and one line not containing said point.
  linearMatroid(m) which returns the linear matroid determined by the matrix m.

License: you can change or modify the contents of this repository under GNU General Public License v3.0.
