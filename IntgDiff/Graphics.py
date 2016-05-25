class Graphics:
	"""
	This class tends to provide basic graphic tools based on `matplotlib` 
	and `mayavi`.
	"""
	def __init__(self, env='numeric', numpoints=50):
		"""
		Accepts one optional argument `env` which determines the types
		of the function to be visualized:
			- 'numeric': is a numerical function (regular python functions)
			- 'sympy': Sympy symbolic function
			- 'sage': Sage symbolic function
		"""
		self.Env = env
		self.NumPoints = numpoints
		self.title = ""
		self.color = []
		self.X = []
		self.Y = []
		self.Z = []
		self.xlabel = '$x$'
		self.ylabel = '$y$'
		self.zlabel = '$z$'
		self.xmin = None
		self.xmax = None
		self.ymin = None
		self.ymax = None
		self.zmin = None
		self.zmax = None
		self.legend = []
		self.grid = True
		self.PlotCount = 0
		self.PlotType = ''

	def SetLabelX(self, lbl):
		"""
		Sets the label for `X` axis
		"""
		self.xlabel = lbl

	def SetLabelY(self, lbl):
		"""
		Sets the label for `Y` axis
		"""
		self.ylabel = lbl

	def SetLabelZ(self, lbl):
		"""
		Sets the label for `Z` axis
		"""
		self.zlabel = lbl

	def SetTitle(self, ttl):
		"""
		Sets the title of the graph.
		"""
		self.title = ttl

	def Plot2D(self, func, xrng, color='blue', legend=""):
		"""
		Appends a curve to the Graphics object. The parameters are as follows:
			- `func`: the function to be plotted,
			- `xrng`: a triple of the form `(x, a, b)`, where `x` is the `func`s
			independents variable, over the range `[a, b]`,
			- `color`: the color of the current curve,
			- `legend`: the text for the legend of the current crve.
		"""
		if self.PlotType == '':
			self.PlotType = '2D'
		elif self.PlotType != '2D':
			raise Exception("Cannot combine 2d and 3d plots")

		if self.Env == 'numeric':
			f_ = func
		elif self.Env == 'sympy':
			from sympy import lambdify
			f_ = lambdify(xrng[0], func, "numpy")
		elif self.Env == 'sage':
			from sage.all import fast_callable
			f_ = fast_callable(func, vars=[xrng[0]])
		else:
			raise Exception("The function type is not recognized. Only 'numeric', 'sympy' and 'sage' are accepted.")

		if self.X == []:
			stp_lngt = float(xrng[2] - xrng[1])/self.NumPoints
			self.X = [xrng[1]+i*stp_lngt for i in range(self.NumPoints+1)]
			self.xmin = xrng[1]
			self.xmax = xrng[2]
		TempY = [f_(x) for x in self.X]
		if self.ymin is None:
			self.ymin = min(TempY)
			self.ymax = max(TempY)
		else:
			self.ymin = min(min(TempY), self.ymin)
			self.ymax = max(max(TempY), self.ymax)
		self.Y.append(TempY)
		self.color.append(color)
		self.legend.append(legend)
		self.PlotCount += 1

	def Plot3D(self, func, xrng, yrng):
		"""
		Sets a surface to the Graphics object. The parameters are as follows:
			- `func`: the function to be plotted,
			- `xrng`: a triple of the form `(x, a, b)`, where `x` is the first 
			`func`s independents variable, over the range `[a, b]`,
			- `yrng`: a triple of the form `(y, c, d)`, where `x` is the second
			`func`s	independents variable, over the range `[c, d]`.
		"""
		import numpy as np
		if self.PlotType == '':
			self.PlotType = '3D'
		elif self.PlotType != '3D':
			raise Exception("Cannot combine 2d and 3d plots")
		
		if self.Env == 'numeric':
			f_ = func
		elif self.Env == 'sympy':
			from sympy import lambdify
			f_ = lambdify((xrng[0], yrng[0]), func, "numpy")
		elif self.Env == 'sage':
			from sage.all import fast_callable
			f_ = fast_callable(func, vars=[xrng[0], yrng[0]])
		else:
			raise Exception("The function type is not recognized. Only 'numeric', 'sympy' and 'sage' are accepted.")
		#self.f_ = f_
		if self.X == []:
			x_stp_lngt = float(xrng[2] - xrng[1])/self.NumPoints
			self.X = [xrng[1]+i*x_stp_lngt for i in range(self.NumPoints+1)]
			self.xmin = xrng[1]
			self.xmax = xrng[2]
			y_stp_lngt = float(yrng[2] - yrng[1])/self.NumPoints
			self.Y = [yrng[1]+i*y_stp_lngt for i in range(self.NumPoints+1)]
			self.ymin = yrng[1]
			self.ymax = yrng[2]
		X = np.array(self.X)
		Y = np.array(self.Y)
		self.X, self.Y = np.meshgrid(X, Y)
		self.Z = f_(self.X, self.Y)
		self.intX, self.intY = np.mgrid[xrng[1]:xrng[2]:x_stp_lngt, yrng[1]:yrng[2]:y_stp_lngt]
		self.intZ = f_(self.intX, self.intY)


	def save(self, fname="fig.png"):
		"""
		Saves the outpu of the `Graphics` object to the file `fname`.
		"""
		import numpy as np
		if self.PlotType == '2D':
			import matplotlib.pyplot as plt
			
			plt.axis([self.xmin, self.xmax, self.ymin, self.ymax])
			if self.title != "":
				plt.title(self.title)
			if self.xlabel != '':
				plt.xlabel(self.xlabel, color='gray')
			if self.ylabel != '':
				plt.ylabel(self.ylabel, color='gray')
			for idx in range(self.PlotCount):
				plt.plot(self.X, self.Y[idx], color=self.color[idx])
			plt.legend(self.legend, loc=2)
			plt.grid(self.grid)

		elif self.PlotType == '3D':
			from mpl_toolkits.mplot3d import Axes3D
			import matplotlib.pyplot as plt
			from matplotlib import cm
			
			fig = plt.figure()
			ax = Axes3D(fig)
			ax.set_xlabel(self.xlabel)
			ax.set_ylabel(self.ylabel)
			ax.set_zlabel(self.zlabel)
			#ax.set_xlim3d(self.xmin, self.xmin)
			#ax.set_ylim3d(self.ymin, self.ymax)
			#ax.set_zlim3d(0, 1)
			ax = fig.gca(projection='3d')
			ax.plot_surface(self.X, self.Y, self.Z, rstride=1, cstride=1, linewidth=.1, cmap=cm.jet)

		plt.savefig(fname)

	def interact(self):
		"""
		Shows an interavtive demo of the 3D surface, using `mayavi`,
		so it requires `mayavi` for python to be installed.
		"""
		if self.PlotType != '3D':
			raise Exception("No 3D surface is available.")
		
		import numpy as np
		from mayavi import mlab
		s = mlab.surf(self.intX, self.intY, self.intZ)
		mlab.show()