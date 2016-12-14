# Solving a pde via collocation method
Symbolic = 'sympy'

if Symbolic == 'sympy':
	from sympy import *
	x = Symbol('x')
	t = Symbol('t')
	y = Function('y')(t, x)
elif Symbolic == 'sage':
	from sage.all import *
	t = var('t')
	x = var('x')
	y = function('y')(t, x)
elif Symbolic == 'symengine':
	from symengine import *
	x = Symbol('x')
	t = Symbol('t')
	y = function_symbol('y', t, x)

from pyProximation import *

n = 3

S = OrthSystem([t, x], [(0, 1), (0, 1)])
B = S.PolyBasis(n)
S.Basis(B)
S.FormBasis()

m = int(sqrt(len(S.OrthBase)))
pnts = [(.2, .4), (.4, .6), (.4, .2), (.6, .2), (.2, .8)]

Z = t*sin(pi*x)
R = diff(Z, t) - diff(Z, x)
if Symbolic == 'sympy':
	EQ1 = Eq(diff(y, t) - diff(y, x), R)
if Symbolic == 'symengine':
	EQ1 = (diff(y, t) - diff(y, x) - R)
elif Symbolic == 'sage':
	EQ1 = diff(y, t) - diff(y, x) == R

C = Collocation([t, x], [y])
C.SetOrthSys(S, y)
C.Equation([EQ1])

if Symbolic == 'sympy':
	C.Condition(Eq(y, 0), [0, 0])
	C.Condition(Eq(y, 0), [0, .3])
	C.Condition(Eq(y, 0), [1, 1])
	C.Condition(Eq(y, 1), [1, .5])
	C.Condition(Eq(y, 0), [0, .7])
elif Symbolic == 'symengine':
	C.Condition(y, [0, 0])
	C.Condition(y, [0, .3])
	C.Condition(y, [1, 1])
	C.Condition(y - 1, [1, .5])
	C.Condition(y, [0, .7])
elif Symbolic == 'sage':
	C.Condition(y == 0, [0, 0])
	C.Condition(y == 0, [0, .3])
	C.Condition(y == 0, [1, 1])
	C.Condition(y == 1, [1, .5])
	C.Condition(y == 0, [0, .7])

#C.CollPoints(pnts)
C.setSolver('scipy')

Apprx = C.Solve()
print Apprx[0]

G = Graphics(Symbolic)
G.SetLabelX("$t$")
G.SetLabelY("$x$")
G.SetLabelZ("$y(t, x)$")
G.Plot3D(Apprx[0], (t, 0, 1), (x, 0, 1))

G.save('Exm3-%d.png'%(n))
G.interact()