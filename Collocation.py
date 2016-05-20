#******************************************************************************#
#       Copyright (C) 2016 Mehdi Ghasemi <mehdi.ghasemi@usask.ca>     	       #
#                                                                              #
#  Distributed under the terms of the GNU General Public License (GPL)         #
#  as published by the Free Software Foundation; either version 2 of           #
#  the License, or (at your option) any later version                          #
#                  http://www.gnu.org/licenses/                                #
#******************************************************************************#

from sage.all import *

def DetSymEnv():
	Env = []
	from sys import modules
	if 'sympy' in modules:
		Env.append('sympy')
	if 'sage' in modules:
		Env.append('sage')
	return Env

if 'sage' in DetSymEnv():
	from sage.all import *
else:
	import sympy

class OrthogonalSystem:
	"""
	`OrthogonalSystem` class produces an orthogonal system of functions
	according to a suggested basis of functions and a given measure 
	supported on a given region.

	This basically performs a 'Gram-Schmidt' method to extract the 
	orthogonal basis. The inner product is obtained by integration of 
	the product of functions with respect to the given measure (more
	accurately, the distribution).

	To initiate an instance of this class one should provide a list of
	symbolic variables `variables` and the range of each variable as a
	list of lists `var_range`.
	"""

	def __init__(self, variables, var_range, env='sage'):
		"""
		To initiate an orthogonal system of functions, one should provide
		a list of symbolic variables `variables` and the range of each 
		these variables as a list of lists `var_range`.
		"""
		assert (type(variables) is list) and (type(var_range) is list), """The OrthogonalSystem class object
		requires two lists as inputs: 1) list od fymbolic variables; 2) range of each variable."""
		self.EnvAvail = DetSymEnv()
		if self.EnvAvail == []:
			raise Exception("No Symbolic tool is available.")
		elif ('sage' in self.EnvAvail) and (env=='sage'):
			self.Env = 'sage'
		else:
			self.Env = 'sympy'
		self._pdf = 1
		self.Vars = variables
		self.num_vars = len(self.Vars)
		self.Domain = var_range
		self.OriginalBasis = []
		self.OrthBase = []
		self.Numerical = False

	####### Symbolic functions from symbolic tools
	def Integral(self, f, x, a, b):
		"""
		Local version of integral
		"""
		if self.Env == 'sage':
			if not self.Numerical:
				return integral(f, (x, a, b))
			else:
				return numerical_integral(f, a, b)[0]
		else:
			return sympy.integrate(f, (x, a, b))

	def expand(self, expr):
		"""
		Local version of expand
		"""
		if self.Env == 'sage':
			return expand(expr)
		else:
			return sympy.expand(expr)

	def sqrt(self, expr):
		"""
		Local version of sqrt
		"""
		if self.Env == 'sage':
			return sqrt(expr)
		else:
			return sympy.sqrt(expr)

	def Neval(self, expr):
		"""
		Local version of numerical evaluation
		"""
		if self.Env == 'sage':
			return N(expr)
		else:
			return sympy.N(expr)
	############################################

	def PolyBasis(self, n):
		"""
		Generates a polynomial basis from variables consisting of all
		monomials of degree at most `n`.
		"""
		assert n>=0, "'n' must be a positive integer."
		from itertools import product
		B = []
		for o in product(range(n+1), repeat=self.num_vars):
			if sum(o) <= n:
				T_ = 1
				for idx in range(self.num_vars):
					T_ *= self.Vars[idx]**o[idx]
				B.append(T_)
		return B

	def FourierBasis(self, n):
		"""
		Generates a Fourier basis from variables consisting of all
		sin & cos functions with coefficients at most `n`.
		"""
		assert n>=0, "'n' must be a positive integer."
		from itertools import product
		B = []
		for o in product(range(n+1), repeat=self.num_vars):
			if sum(o) <= n:
				SinCos = product(range(2), repeat=self.num_vars)
				for ex in SinCos:
					T_ = 1
					for idx in range(self.num_vars):
						period = self.Domain[idx][1]-self.Domain[idx][0]
						if o[idx] != 0:
							if ex[idx] == 0:
								T_ *= cos(2*pi*o[idx]*self.Vars[idx]/period)
							else:
								T_ *= sin(2*pi*o[idx]*self.Vars[idx]/period)
					B.append(T_)
		return list(set(B))

	def SetEnv(self, env):
		"""
		Sets the symbolic tool environment.
		"""
		env = env.lower()
		if env not in self.EnvAvail:
			raise Exception("The selected symbolic tool is not available")
		else:
			self.Env = env

	def pdf(self, dm):
		"""
		To set the measure which the orthogonal system will be computed,
		simply call this method with the corresponding distribution as
		its parameter `dm`; i.e, the parameter is 'd(m)'' where `m' is 
		the original measure.
		"""
		self._pdf = dm

	def Basis(self, base_set):
		"""
		To specify a particular family of function as a basis, one should
		call this method with a list `base_set` of linearly independent 
		functions.
		"""
		self.OriginalBasis = base_set
		self.num_base = len(self.OriginalBasis)

	def FormBasis(self):
		"""
		Call this method to generate the orthogonal basis corresponding
		to the given basis via `Basis` method.
		The result will be stored in a property called `OrthBase` which
		is a list of function that are orthogonal to each other with
		respect to the measure `_pdf` over the given range `Domain`.
		"""
		for f in self.OriginalBasis:
			nf = 0
			for u in self.OrthBase:
				nf += self.project(f, u)
			nf = f - nf
			F = self.expand(nf/self.sqrt(self.inner(nf, nf)))
			self.OrthBase.append(F)

	def inner(self, f, g):
		"""
		Computes the inner product of the two parameters with respect to
		the measure `_pdf`.
		"""
		F = f*g*self._pdf
		for v in range(self.num_vars):
			F = self.Integral(F, self.Vars[v], self.Domain[v][0], self.Domain[v][1])
		return F

	def project(self, f, g):
		"""
		Finds the projection of `f` on `g` with respect to the inner 
		product induced by the measure `_pdf`.
		"""
		return g*self.inner(f, g)/self.inner(g, g)

	def Series(self, f):
		"""
		Given a function `f`, this method finds and returns the 
		coefficients of the	series that approximates `f` as a 
		linear combination of the elements of the orthogonal basis.
		"""
		cfs = []
		for b in self.OrthBase:
			cfs.append(self.Neval(self.inner(f, b)))
		return cfs

##################################################################################################

class Collocation:
	"""
	`Collocation` class tries to approximate the solutions of a system
	of partial differential equations with respect to an orthogonal
	system of functions.

	To initiate an instance of this class one needs to provide two set
	of parameters:
		1) List of independent symbolic variables `variables`;
		2) List of unknown functions to be found that depend on the
		independent variables `ufunc`.
	"""
	def __init__(self, variables, ufunc):
		self.Env = DetSymEnv()
		if self.Env == []:
			raise Exception("No Symbolic tool is available.")
		self.Vars = variables
		self.num_vars = len(variables)
		self.uFuncs = ufunc
		self.num_funcs = len(ufunc)
		self.degree = 1
		self.Cnds = []
		self.CndVals = []
		self.Coeffs = {}
		self.Apprx = []
		self.Points = []
		self.Solver = 'sage'

	def SetOrthSys(self, obj):
		"""
		To approximate the solutions of the system of pdes, the class 
		requires an orthogonal system of functions `OrthogonalSystem`.
		This method accepts such a system.
		"""
		self.OrthSys = obj
		self.degree = self.OrthSys.num_base

	def Equation(self, eq):
		"""
		To enter the system of equations, use this meyhod with a list 
		of equations as input.
		"""
		self.Eq = eq

	def Condition(self, eq, val):
		"""
		List of initial and boundary conditions.
		"""
		#if eq not in self.Cnds:
		self.Cnds.append(eq)
		self.CndVals.append(val)

	def setSolver(self, solver):
		"""
		Currently only two solvers are supported:
			1) `sage`s defult solver for rather simple system of algebraic
			equations.
			2) `scipy`s `fsolves` to handel more complex and larger systems.
		"""
		self.Solver = solver.lower()

	def CollPoints(self, pnts):
		"""
		Accepts alist of collocation point `pnts`, to form the algebraic 
		system of equations and find the coefficients of the orthogonal 
		functions from `OrthogonalSystem.OrthBase`.
		"""
		self.Points += pnts

	def collocate(self):
		"""
		Internal use: generates the system of equations for coefficients 
		using the collocation points.
		"""
		if 'sage' not in self.Env:
			from sympy import Symbol as var
		else:
			from sage.all import var

		var_syms = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
		self.CF = {var_syms[s]:[var('%s%d'%(var_syms[s], i)) for i in range(self.degree)] for s in range(self.num_funcs)}
		self.SymCF = []
		for s in self.CF:
			self.SymCF += self.CF[s]
		self.SR = {}
		self.REq = []
		for f_idx in range(self.num_funcs):
			s = var_syms[f_idx]
			T = 0
			for i in range(self.degree):
				T += self.CF[s][i]*(self.OrthSys.OrthBase[i])
			self.SR[s] = T
		for eq in self.Eq:
			f_idx = 0
			Teq = eq
			for f in self.uFuncs:
				for v in self.Vars:
					for d_ord in range(1,6):
						Teq = Teq.subs({diff(f, v, d_ord):diff(self.SR[var_syms[f_idx]], v, d_ord)})
					Teq = expand(Teq.subs({f:self.SR[var_syms[f_idx]]}))
				f_idx += 1
			if Teq not in self.REq:
				self.REq.append(Teq)
		cnd_idx = 0
		for eq in self.Cnds:
			f_idx = 0
			Teq = eq
			for f in self.uFuncs:
				for v_idx in range(self.num_vars):
					v = self.Vars[v_idx]
					Teq = Teq.subs({diff(f, v):diff(self.SR[var_syms[f_idx]], v)})
					Teq = expand(Teq.subs({f:self.SR[var_syms[f_idx]]}))
				f_idx += 1
			Teq = expand(Teq.subs({self.Vars[v]:self.CndVals[cnd_idx][v] for v in range(self.num_vars)}))
			if Teq not in self.REq:
				self.REq.append(Teq)
			cnd_idx += 1

	def PlugPoints(self):
		"""
		Internal use: plug in collocation points to elliminate independent 
		variables and keep the coefficients.
		"""
		numeric_eqs = []
		for p in self.Points:
			chg = {self.Vars[i]:p[i] for i in range(self.num_vars)}
			for eq in self.REq:
				Teq = eq.subs(chg)
				if Teq not in numeric_eqs:
					numeric_eqs.append(Teq)
			if len(numeric_eqs) >= len(self.SymCF):
				break
		#print len(numeric_eqs), len(self.SymCF)
		if len(numeric_eqs) != len(self.SymCF):
			raise Exception("Number of points and equations are not equal! Check the conditions.")
		if self.Solver == 'sage':
			sols = solve(numeric_eqs, self.SymCF, solution_dict=True)
			sols = sols[0]
			return sols
		elif self.Solver in ['scipy', 'newton_krylov']:
			import scipy 
			from scipy import optimize as opt
			def f(x):
				chng = {}
				U = self.SymCF
				n_var = len(U)
				chng = {U[i]:float(x.item(i)) for i in range(n_var)}
				EQs_ = []
				for eq in numeric_eqs:
					teq = eq.lhs()-eq.rhs()
					EQs_.append(teq.subs(chng).n())
				return EQs_
			nvars = len(self.SymCF)
			if self.Solver == 'scipy':
				sol = list(opt.fsolve(f,[0 for _ in range(nvars)]))
			elif self.Solver == 'newton_krylov':
				sol = list(opt.newton_krylov(f,[0 for _ in range(nvars)]))
			sols = {self.SymCF[i]: sol[i] for i in range(nvars)}
			return sols
	
	def Solve(self):
		"""
		Solves the collocation equations and keep a dictionary of 
		coefficients in `self.Coeffs` and returns a list of functions
		in the span of orthoginal system.
		"""
		
		CltnPnts = Sample(self.OrthSys.Vars, self.OrthSys.Domain)
		CltnPnts.setPDF(self.OrthSys._pdf)
		num = self.degree - len(self.Points)#+len(self.Cnds))
		#print self.degree, len(self.Points)+len(self.Cnds)
		if num < 0:
			raise Exception("Too many points are associated. Reduce at least %d"%(-num))
		cl_points = CltnPnts.sample(num)
		self.CollPoints(cl_points)
		print "Start..."
		self.collocate()
		print "Solving ..."
		self.Coeffs = self.PlugPoints()

		for s in self.SR:
			self.Apprx.append(self.SR[s].subs(self.Coeffs))
		return self.Apprx

##################################################################################################

class Sample():
	"""
	`Sample` accepts a distribution through `setPDF` and a list of 
	symbolic variables and range of each at initialization:
		`Sample(list of variables, list of ranges of variables)`

	Then calling the method `sample(n)`, returns `n` points from the 
	given domain which follows the given distribution.
	"""
	def __init__(self, variables, region):
		"""
		Initializes the `Sample` object with list of variables `variables`
		and a list of lists `region` to determine the domain.
		"""
		self.Vars = variables
		self.dim = len(self.Vars)
		self.Lebesque = 1
		self.Region = region
		for rng in region:
			self.Lebesque *= (rng[1] - rng[0])
		
		self.pdf = 1./self.Lebesque
		self.SubRegions = {}
		self.NumSample = {}
		self.Method = 'fast'
		return
	####### Symbolic functions from symbolic tools
	def Integral(self, f, x, a, b):
		"""
		Local version of integral
		"""
		if 'sage' in DetSymEnv():
			return integral(f, (x, a, b))
		else:
			return sympy.integrate(f, (x, a, b))
	###############################################
	def setPDF(self, p):
		"""
		Sets the probability distribution function over the domain.
		"""
		self.pdf = p

	def Measure(self, rgn):
		"""
		Approximates or Determines accurately the measure of the given 
		subregion:
			If `self.Method` is set to be equal to `fast`, the method 
			uses the average value theorem to approximate the wight of
			the region. Otherwise, it uses integration to return the 
			exact value.
		"""
		p = self.pdf
		if self.Method == 'fast':
			pnt = {}
			vol = 1
			idx = 0
			for v in self.Vars:
				pnt[v] = (rgn[idx][1]+rgn[idx][0])/2.
				vol *= (rgn[idx][1] - rgn[idx][0])
			if type(p) is type(self.Vars[0]):
				rs = N(vol*p.subs(pnt))
			else:
				rs = N(vol*p)
		else:
			idx = 0
			for v in self.Vars:
				p = self.Integral(p, v, rgn[idx][0], rgn[idx][1])
				idx += 1
			rs = N(p)
		return rs

	def NormalizePDF(self):
		"""
		Normalizes the given pdf to make sure that the measure of the
		domain sums up to 1.
		"""
		p = self.pdf
		idx = 0
		for v in self.Vars:
			p = self.Integral(p, v, self.Region[idx][0], self.Region[idx][1])
			idx += 1
		M = N(p)
		self.pdf = self.pdf / M

	def random_points(self, n, rgn):
		"""
		Returns `n` points uniformly distributed on the region `rgn`.
		"""
		# Check for inputs
		assert (n>=0), "'n' must be a non-negative integer number"
		assert type(rgn) is list, "The range of variables must be given as a list of lists."
		#
		from random import uniform
		pnts = []
		while len(pnts)<n:
			v = []
			for rng in rgn:
				v.append(uniform(rng[0], rng[1]))
			v = tuple(v)
			if v not in pnts:
				pnts.append(v)
		return pnts

	def Partition(self, m):
		"""
		Partitions the domain into at least `m` equally sliced subdomains.
		"""
		# Check for the input
		assert (m>0), "'m' must be a positive integer."
		#
		from math import ceil, pow
		from itertools import product
		n = int(ceil(pow(m, (1./self.dim))))
		delta = [(r[1]-r[0])/float(n) for r in self.Region]
		for o in product(range(n), repeat=self.dim):
			SR = [[self.Region[i][0]+o[i]*delta[i], self.Region[i][0]+(o[i]+1)*delta[i]] for i in range(self.dim)]
			self.SubRegions[o] = SR

	def SizeSubRegion(self, n):
		"""
		Approximates the number of points to be sellected uniformly from 
		any subregion such that the total number exeeds `n`
		"""
		from math import ceil
		num = max(n, len(self.SubRegions))
		for o in self.SubRegions:
			self.NumSample[o] = ceil(num*self.Measure(self.SubRegions[o]))#max(1, round(num*self.Measure(self.SubRegions[o])))

	def sample(self, num):
		"""
		Returns a list of `num` points from the given domain according
		to the pdf.
		"""
		assert num>=1, "Sample size must be a positive integer number."
		import random
		if num<= 0:
			return []
		self.Partition(num)
		self.SizeSubRegion(num)
		pnts = []
		for o in self.NumSample:
			pnts += self.random_points(self.NumSample[o], self.SubRegions[o])
		return random.sample(set(pnts), num)