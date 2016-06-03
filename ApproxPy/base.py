class Foundation:
	"""
	This class contains common features of all modules.
	"""
	def __init__(self):
		pass

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