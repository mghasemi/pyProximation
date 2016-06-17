=====================
Introduction
=====================

This is a brief documentation for using *pyProximation*.
pyProximation was originally written to solve systems of integro-differential
equations via collocation method in arbitrary number of variables.
The current implementation of the method is based on a finite dimensional
orthogonal system of functions.

Requirements and dependencies
===============================

Since this is a python package, so clearly one needs have python available.
pyProximation relies on the following packages:

	+ for numerical calculation:
		- `numpy <http://www.numpy.org/>`_,
		- `scipy <https://www.scipy.org/>`_ ,
	+ for symbolic computation, either:
		- `sympy <http://www.sympy.org/>`_ or,
		- `sage <http://www.sagemath.org/>`_ or,
		- `symengine <https://github.com/symengine/symengine>`_,
	+ for visualization:
		- `mathplotlib <http://matplotlib.org/>`_,
		- `mayavi <http://mayavi.sourceforge.net/>`_.

For the symbolic computation, pyProximation need the presence of either `sympy` or `sage`.

Download
================

`pyProximation` can be obtained from `https://github.com/mghasemi/pyProximation <https://github.com/mghasemi/pyProximation>`_.

Installation
=========================

To install `pyProximation`, run the following in terminal::

	sudo python setup.py install

Documentation
--------------------------
The documentation of `pyProximation` is prepaqred via `sphinx <http://www.sphinx-doc.org/>`_.

License
=======================
`pyProximation` is distributed under MIT license:

MIT License
------------------

	Copyright (c) 2016 Mehdi Ghasemi

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.