=============================
Symbolic toolkits
=============================

The 'pyProximation' library performs symbolic computations as well as numerical calculations.
To perform numerical calculations, it relies on the excelent libraries 'NumPy' and 'SciPy'.
There are various toolkits that are able to handel symbolic computations and each one is suitable 
for certain purposes. 'pyProximation' takes advantage of three different toolkits:

	1. `sympy`
	2. `sage`
	3. `symengine`

SymPy
=============================
`SymPy <http://sympy.org/>`_ is a Python library for symbolic computation. It provides computer algebra capabilities either 
as a standalone application, as a library to other applications, or live on the web as `SymPy Live <http://live.sympy.org/>`_ 
or `SymPy Gamma <http://www.sympygamma.com/>`_. SymPy is trivial to install and to inspect because is written entirely in Python 
with few dependencies. This ease of access combined with a simple and extensible code base in a well known language make SymPy a 
computer algebra system with a relatively low barrier to entry.

SymPy includes features ranging from basic symbolic arithmetic to calculus, algebra, discrete mathematics and quantum physics. 
It is capable of formatting the result of the computations as LaTeX code.

SymPy is free software and is licensed under New BSD License. The lead developers are Ondřej Čertík and Aaron Meurer.

Usage
----------------------------
Sympy can be imported by::

	from sympy import *

To introduce a symbolic variable called `x`, use::

	x = Symbol('x')

and to define a symbolic function `y` depending in symbolic variables `x`, `t`, use::

	y = Function('y')(x, t)

To introduce an equation like :math:`y(x,t)=xt`, use `Eq`::

	equation = Eq(y, x*t)

For differentiation of `y` with respect to `x` of order `n`, i.e. :math:`\frac{\partial^ny}{\partial x^n}`, enter::

	diff(y, x, n)

For integration of `y` with respect to `x`, i.e., :math:`\int y dx`, write::

	integrate(y, x)

and for definite integral :math:`\int_a^b y dx` write::

	integrate(y, (x, a, b))

Sage
============================
`SageMath <http://www.sagemath.org/>`_ (previously Sage or SAGE, System for Algebra and Geometry Experimentation) is a free 
open-source mathematics software system licensed under the GPL with features covering many aspects of mathematics, including 
algebra, combinatorics, numerical mathematics, number theory, and calculus.
It builds on top of many existing open-source packages: NumPy, SciPy, matplotlib, Sympy, Maxima, GAP, FLINT, R and many more. 
Access their combined power through a common, Python-based language or directly via interfaces or wrappers.

Usage
----------------------------
When present, Sage can be imported by::

	from sage.all import *

To introduce a symbolic variable called `x`, use::

	x = var('x')

and to define a symbolic function `y` depending in symbolic variables `x`, `t`, use::

	y = function('y')(x, t)

To introduce an equation like :math:`y(x,t)=xt`, use `==`::

	equation = (y == x*t)

For differentiation of `y` with respect to `x` of order `n`, i.e. :math:`\frac{\partial^ny}{\partial x^n}`, enter::

	diff(y, x, n)

For integration of `y` with respect to `x`, i.e., :math:`\int y dx`, write::

	integral(y, x)

and for definite integral :math:`\int_a^b y dx` write::

	integral(y, (x, a, b))

SymEngine
============================
`SymEngine <https://github.com/symengine/symengine>`_ is a standalone fast C++ symbolic manipulation library. 
Optional thin wrappers allow usage of the library from other languages, e.g.:

	+ C wrappers allow usage from C, or as a basis for other wrappers;
	+ Python wrappers allow easy usage from Python and integration with SymPy and Sage;
	+ Ruby wrappers;
	+ Julia wrappers;
	+ Haskell wrappers.

It is licensed under MIT license.

Usage
----------------------------
SymEngine can be imported by::

	from symengine import *

To introduce a symbolic variable called `x`, use::

	x = Symbol('x')

and to define a symbolic function `y` depending in symbolic variables `x`, `t`, use::

	y = function_symbol('y', x, t)

`SymEngine` does not have specific command for defining an equation, so work with equations, one should always equalities with 0.
To introduce an equation like :math:`y(x,t)=xt`, use::

	equation = y - x*t

For differentiation of `y` with respect to `x`, i.e. :math:`\frac{\partial y}{\partial x}`, enter::

	diff(y, x)

To achieve a certain order, repeat the above code as many times as necessary.

The current version of `SymEngine` does not support integration.