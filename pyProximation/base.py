class Foundation:
	"""
	This class contains common features of all modules.
	"""
	def __init__(self):
		pass

	def DetSymEnv(self):
		"""
		Returns a list. The list consists of all symbolic tools 
		present among 'sympy','sage' and 'symengine'.
		"""
		Env = []
		from sys import modules
		if 'sympy' in modules:
			Env.append('sympy')
		if 'sage' in modules:
			Env.append('sage')
		if 'symengine' in modules:
			Env.append('symengine')
		return Env

	def CommonSymFuncs(self, env):
		if env == 'sympy':
			from sympy import expand, sqrt, sin, cos, pi, diff
			from sympy import Symbol as Symbol
			#from sympy import Function as SymFunc
		elif env == 'sage':
			from sage.all import expand, sqrt, sin, cos, pi, diff
			from sympy import var as Symbol
			#from sympy import function as SymFunc
		elif env == 'symengine':
			from symengine import expand, sqrt, sin, cos, pi, diff
			from sympy import Symbol as Symbol
			#from sympy import function_symbol as SymFunc

		self.expand = expand
		self.sqrt = sqrt
		self.sin = sin
		self.cos = cos
		self.pi = pi
		self.diff = diff
		self.Symbol = Symbol
