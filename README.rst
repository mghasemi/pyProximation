========================
pyProximation Package
========================

This package was originally written to solve systems of integro-differential
equations via collocation method in arbitrary number of variables.
The current implementation of the method is based on a finite dimensional
orthogonal system of functions. Therefore additional modules were required 
to achieve this goal.

The modules are listed below:

measure
========================

Consists of one class called `Measure` which performs measure theoretic 
computations such as:
	- calculate the measure of a set according to the given measure;
	- calculate the integral of a function w.r.t. the measure;
	- calculate the norm-p of a given function;
	- generate samples from the support of the measure according to the distribution.

orthsys
========================
Consists of a class called `OrthSystem` which produces orthogonal systems
of functions according to a basis provided as a priori. It takes a measure
on a given set and a finite set of independent functions such as polynomials
or Fourier series, and generates an orthonormal basis w.r.t. the measure.
Then, given a function, it also computes coefficients of the corresponding
series.

interpolation
========================
The interpolation module consists of `Interpolation` class. This class
provides multi-variate interpolation.

collocation
========================

The class `Collocation`, takes a list of symbolic integro-differential 
equations, an orthonormal system of functions and some collocation points 
(if provided) then finds approximate solutions of the system.

subregion
========================
The `SubRegion` class partitions the region into sub-regions, solve
the system of Integro-differential equations on each and glue them
together.

Graphics
========================

This module provides basic command to illustrate 2D and 3D plots of 
symbolic or numeric functions based on `matplotlib`. The 3D plots also
include interactive ones based on `mayavi`.


For more detals refer to the `documentation <http://pyproximation.readthedocs.io/>`_.

To build the documentation run::

	$ sphinx-build -b html ./doc/ ./doc/_build

or from ``./doc`` run ``$ make html``.

To generate a pdf file run ``$ make latexpdf``. This requires `latex` to be available.
