import sympy as sym

# Here we find the derivative of the survival function used in the
# article in order to obtain the probability density function

t, rho, p, x, lamb, mu = sym.symbols('t rho p x lambda mu')

p1 = rho * p

p2 = sym.exp(-lamb * t) * sym.exp(-mu * t)

p3 = (1 - sym.exp(-lamb * t)) * sym.exp(-mu * t)

survival = sym.exp(-p1 * (1 - p2/(1-p3)))

density = sym.diff(1 - survival, t)

sym.simplify(survival)
sym.simplify(density)
sym.simplify(sym.log(survival))
sym.simplify(sym.log(density))
