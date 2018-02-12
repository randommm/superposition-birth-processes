import sympy as sym

# Here we find the derivative of the survival function used in the
# article in order to obtain the probability density function

t, a, beta, x, lamb, mu = sym.symbols('t a beta x lambda mu')

p1 = sym.exp(a+beta*x)

p2 = sym.exp(-lamb * t) * sym.exp(-mu * t)

p3 = (1 - sym.exp(-lamb * t)) * sym.exp(-mu * t)

survival = sym.exp(-p1 * (1 - p2/(1-p3)))

density = sym.diff(1 - survival, t)

sym.simplify(survival)
sym.simplify(density)
sym.simplify(sym.log(survival))
sym.simplify(sym.log(density))
