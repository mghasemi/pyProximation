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