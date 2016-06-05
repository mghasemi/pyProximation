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
		- `sage <http://www.sagemath.org/>`_,
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