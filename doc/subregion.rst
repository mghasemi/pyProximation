=============================
Distributing over subregions
=============================

Sometimes the system of equations is too complicated and it might be very costly to produce a reasonably good 
estimation by raising the size of the function basis. A possible solution is to divide the domain into various
smaller subregions and try to solve the system for those smaller subregions. Use the partial solutions to force
appropriate boundary condition to the remaining subregions in order to get consistent solutions on all subregions.
This task is done by use of the `SubRegion` class.

.. note::
	Currently, this only works with symbolic package `sympy`.


