class Measure:
	"""
	An instance of this class is a measure on a given set `supp`.
	The support is either 
		-- a python variable of type `set`, or
		-- a list of tuples which represents a box in euclidean space.
	"""
	def __init__(self, dom, w=None):
		"""
		Initializes a measure object according to the inputs:
			--`dom` must be either
				-- a list of 2-tuples
				-- a non-empty dictionary
			-- `w` must be a
				-- a function if `dom` defines a region
				-- left blank (None) if `dom` is a dictionary
		"""
		self.ErrorMsg = ""
		self.dim = 0
		self.card = 0
		if not self.check(dom, w):
			raise Exception(self.ErrorMsg)

	def boxCheck(self, B):
		"""
		Checks the structure of the box `B`.
		"""
		flag = True
		for interval in B:
			flag = flag and (type(interval)==tuple) and (len(interval)==2)
		return flag

	def check(self, dom, w):
		"""
		Checks the input types and their consistency
		"""
		from types import FunctionType, IntType, LongType, FloatType
		if type(dom) == list:
			if not self.boxCheck(dom):
				self.ErrorMsg = "Each member of the support's list must be a tuple of 2 elements"
				return False
			self.DomType = "box"
			self.dim = len(dom)
			self.supp = dom
			if type(w) not in [FunctionType, IntType, LongType, FloatType]:
				self.ErrorMsg = "Weight must be a `function` defined over the support or a number"
				return False
			elif type(w)==FunctionType:
				self.weight = w
			else:
				self.weight = lambda *x: w
		elif type(dom) == dict:
			if len(dom)==0:
				self.ErrorMsg = "A measure can not have an empty support"
				return False
			self.DomType = "set"
			self.card = len(dom)
			self.supp = dom.keys()
			self.weight = dom
		else:
			self.ErrorMsg = "First parameter must be either a list of 2-tupels or a dictionary"
			return False
		return True

	def measure(self, S):
		"""
		Returns the measure of the set `S`.
		"""
		m = 0
		if self.DomType == "set":
			if type(S) not in [set, list, tuple]:
				raise Exception("The type of the given set must be either `set`, `list` or `tuple`")
			for p in S:
				if p in self.supp:
					m += self.weight[p]
		else:
			if (type(S) != list) or (not self.boxCheck(S)):
				raise Exception("The given set must be a list of 2-tuples")
			from scipy import integrate
			m = integrate.nquad(self.weight, S)[0]
		return m

	def integral(self, f):
		"""
		Returns the integral of `f` with respect to the currwnt measure 
		over the support.
		"""
		from types import FunctionType
		m = 0
		if self.DomType == "set":
			if type(f) not in [dict, FunctionType]:
				raise Exception("The integrand must be a `function` or a `dict`")
			if type(f) == dict:
				for p in self.supp:
					if p in f:
						m += self.weight[p]*f[p]
			else:
				for p in self.supp:
					try:
						m += self.weight[p]*f(p)
					except:
						pass
		else:
			if type(f) != FunctionType:
				raise Exception("The integrand must be a `function`")
			from scipy import integrate
			fw = lambda *x: f(*x)*self.weight(*x)
			m = integrate.nquad(fw, self.supp)[0]
		return m

	def norm(self, p, f):
		"""
		Computes the norm-p of the `f` with respect to the current measure.
		"""
		from math import pow
		absfp = lambda *x: pow(abs(f(*x)), p)
		return pow(self.integral(absfp), 1./p)

	def sample(self, num):
		"""
		Samples from the support according to the measure.

		"""
		assert num>=1, "Sample size must be a positive integer number."
		if self.DomType == 'box':
			from math import ceil, pow
			from itertools import product
			import random
			from random import uniform
			SubRegions = {}
			NumSample = {}
			points = []
			n = int(ceil(pow(num, (1./self.dim))))
			delta = [(r[1]-r[0])/float(n) for r in self.supp]
			for o in product(range(n), repeat=self.dim):
				SubRegions[o] = [(self.supp[i][0]+o[i]*delta[i], self.supp[i][0]+(o[i]+1)*delta[i]) for i in range(self.dim)]
			numpnts = max(num, len(SubRegions))
			muSupp = self.measure(self.supp)
			for o in SubRegions:
				NumSample[o] = ceil(numpnts*self.measure(SubRegions[o]))
			for o in NumSample:
				pnts = []
				while len(pnts)<NumSample[o]:
					v = []
					for rng in SubRegions[o]:
						v.append(uniform(rng[0], rng[1]))
					v = tuple(v)
					if v not in pnts:
						pnts.append(v)
				points += pnts
			return random.sample(set(points), num)
##################################################################################################
class OrthSystem:
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

	def __init__(self, variables, var_range, env='sympy'):
		"""
		To initiate an orthogonal system of functions, one should provide
		a list of symbolic variables `variables` and the range of each 
		these variables as a list of lists `var_range`.
		"""
		assert (type(variables) is list) and (type(var_range) is list), """The OrthSystem class object
		requires two lists as inputs: (1) list od fymbolic variables; (2) range of each variable."""
		self.EnvAvail = self.DetSymEnv()
		if self.EnvAvail == []:
			raise Exception("No Symbolic tool is available.")
		elif (env in self.EnvAvail):
			self.Env = env
		else:
			raise Exception("The selected symbolic tool is not supported.")
		self.Vars = variables
		self.num_vars = len(self.Vars)
		self.Domain = var_range
		self.measure = Measure(self.Domain, 1)
		self.OriginalBasis = []
		self.OrthBase = []
		self.Numerical = False

	def DetSymEnv(self):
		"""
		Returns a list. The list consists of all symbolic tools 
		present among 'sympy' and 'sage'.
		"""
		Env = []
		from sys import modules
		if 'sympy' in modules:
			Env.append('sympy')
		if 'sage' in modules:
			Env.append('sage')
		return Env

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
		if self.Env == "sympy":
			from sympy import sin, cos, pi
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

	def SetMeasure(self, M):
		"""
		To set the measure which the orthogonal system will be computed,
		simply call this method with the corresponding distribution as
		its parameter `dm`; i.e, the parameter is 'd(m)'' where `m' is 
		the original measure.
		"""
		assert isinstance(M, Measure), "The argument must be a `Measure`."
		self.measure = M

	def Basis(self, base_set):
		"""
		To specify a particular family of function as a basis, one should
		call this method with a list `base_set` of linearly independent 
		functions.
		"""
		assert type(base_set) is list, "A list of symbolic functions is expected."
		self.OriginalBasis = base_set
		self.num_base = len(self.OriginalBasis)

	def inner(self, f, g):
		"""
		Computes the inner product of the two parameters with respect to
		the measure `_pdf`.
		"""
		if self.Env == "sympy":
			from sympy import lambdify
			F = lambdify(self.Vars, f*g, "numpy")
			m = self.measure.integral(F)
		return m

	def project(self, f, g):
		"""
		Finds the projection of `f` on `g` with respect to the inner 
		product induced by the measure `_pdf`.
		"""
		return g*self.inner(f, g)/self.inner(g, g)

	def FormBasis(self):
		"""
		Call this method to generate the orthogonal basis corresponding
		to the given basis via `Basis` method.
		The result will be stored in a property called `OrthBase` which
		is a list of function that are orthogonal to each other with
		respect to the measure `_pdf` over the given range `Domain`.
		"""
		if self.Env == 'sympy':
			from sympy import expand, sqrt
		for f in self.OriginalBasis:
			nf = 0
			for u in self.OrthBase:
				nf += self.project(f, u)
			nf = f - nf
			F = expand(nf/sqrt(self.inner(nf, nf)))
			self.OrthBase.append(F)

	def Series(self, f):
		"""
		Given a function `f`, this method finds and returns the 
		coefficients of the	series that approximates `f` as a 
		linear combination of the elements of the orthogonal basis.
		"""
		cfs = []
		for b in self.OrthBase:
			cfs.append(self.inner(f, b))
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
	def __init__(self, variables, ufunc, env='sympy'):
		Env = self.DetSymEnv()
		if Env == []:
			raise Exception("No Symbolic tool is available.")
		if env not in Env:
			raise Exception("The selected symbolic tool is not available.")
		self.Env = env # The selected tool for symbolic computations
		self.Vars = variables # Symbolic variables
		self.num_vars = len(variables) # Number of symbolic variables
		self.uFuncs = ufunc # Unknown functions
		self.num_funcs = len(ufunc) # Number of unknown functions
		self.degree = 1 # Number of elements in the orthogonal basis
		self.EQs = [] # Lists of functional equations
		self.Cnds = [] # Storage for initial and boundary conditions
		self.CndVals = [] # Storage for the values of CND
		self.Coeffs = {} 
		self.Apprx = [] # Approximate solutions to uFuncs
		self.Points = [] # Collocation points
		self.Solver = 'scipy' # The solver to find the roots

	def DetSymEnv(self):
		"""
		Returns a list. The list consists of all symbolic tools 
		present among 'sympy' and 'sage'.
		"""
		Env = []
		from sys import modules
		if 'sympy' in modules:
			Env.append('sympy')
		if 'sage' in modules:
			Env.append('sage')
		return Env

	def SetOrthSys(self, obj):
		"""
		To approximate the solutions of the system of pdes, the class 
		requires an orthogonal system of functions `OrthogonalSystem`.
		This method accepts such a system.
		"""
		assert isinstance(obj, OrthSystem), "An object of type `OrthSystem` is expected."
		self.OrthSys = obj
		self.degree = self.OrthSys.num_base

	def Equation(self, eq):
		"""
		To enter the system of equations, use this meyhod with a list 
		of equations as input.
		"""
		if type(eq) is list:
			self.EQs += eq
		else:
			self.EQs.append(eq)

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
		to be used via collocation points.
		"""
		if self.Env == 'sympy':
			from sympy import Symbol as var
			from sympy import Subs, expand, diff
		elif self.Env == 'sage':
			from sage.all import var
		# symbols for coefficients
		var_syms = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 
		'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w'
		'x', 'y', 'z']
		# Produce enough symbolic variables
		self.CF = {var_syms[s]:[var('%s%d'%(var_syms[s], i)) for i in range(self.degree)] for s in range(self.num_funcs)}
		self.SymCF = []
		for s in self.CF:
			self.SymCF += self.CF[s]
		self.SR = {}
		self.REq = []
		# loop over unknown functions to be found
		for f_idx in range(self.num_funcs):
			s = var_syms[f_idx]
			T = 0
			# loop over elements of orthogonal basis
			for i in range(self.degree):
				T += self.CF[s][i]*(self.OrthSys.OrthBase[i])
			self.SR[s] = T
		# loop over entered equations
		for eq in self.EQs:
			f_idx = 0
			Teq = eq
			for f in self.uFuncs:
				for v in self.Vars:
					if self.Env == 'sage':
						for d_ord in range(1,6):
							Teq = Teq.subs({diff(f, v, d_ord):diff(self.SR[var_syms[f_idx]], v, d_ord)})
					Teq = expand(Teq.subs({f:self.SR[var_syms[f_idx]]})).doit()
				f_idx += 1
			if Teq not in self.REq:
				self.REq.append(Teq)
		# loop over initial and boundary conditions
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
		# Plug in the collocation points to form the algebraic equations
		numeric_eqs = []
		for p in self.Points:
			chg = {self.Vars[i]:p[i] for i in range(self.num_vars)}
			for eq in self.REq:
				Teq = eq.subs(chg)
				if Teq not in numeric_eqs:
					numeric_eqs.append(Teq)
			if len(numeric_eqs) >= len(self.SymCF):
				break
		
		if len(numeric_eqs) != len(self.SymCF):
			raise Exception("Number of points and equations are not equal! Check the conditions.")
		# Solve the algebraic equations
		if self.Solver == 'sage':
			if self.Env != 'sage':
				raise Exception("Sage solver is not available in selected symbolic environment.")
			sols = solve(numeric_eqs, self.SymCF, solution_dict=True)
			sols = sols[0]
			return sols
		elif self.Solver in ['scipy', 'newton_krylov']:
			from scipy import zeros
			from scipy import optimize as opt
			if self.Env == 'sympy':
				from sympy import lambdify
				f_ = [lambdify(self.SymCF, (eq.lhs-eq.rhs), "numpy") for eq in numeric_eqs]
				def f(x):
					z = tuple(float(x.item(i)) for i in range(len(self.SymCF)))
					return [fn(*z) for fn in f_]
			elif self.Env == 'sage':
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
				sol = list(opt.fsolve(f, tuple(0. for _ in range(nvars))))
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
		
		num = self.degree - len(self.Points)
		# Check for too many points
		if num < 0:
			raise Exception("Too many points are associated. Reduce at least %d"%(-num))
		cl_points = []
		# Add enough random points to match up for variables
		if num>0:
			cl_points = self.OrthSys.measure.sample(num)
		self.CollPoints(cl_points)

		self.collocate()
		
		self.Coeffs = self.PlugPoints()

		for s in self.SR:
			self.Apprx.append(self.SR[s].subs(self.Coeffs))
		return self.Apprx