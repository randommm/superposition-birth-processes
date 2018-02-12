import numpy as np
import scipy.stats as stats

# Here we demonstre some properties of a Yule Simon process

d = 3
lambda_ = .4

a = [stats.geom.rvs(np.exp(-10.6 * lambda_), size=d).sum() for _ in range(10000)]
a = np.array(a)
print(a.mean())
print(a.var())

b = [stats.geom.rvs(np.exp(-0.4 * lambda_), size=d).sum() for _ in range(10000)]
b = [stats.geom.rvs(np.exp(-0.2 * lambda_), size=x).sum() for x in b]
b = np.array(b)
print(b.mean())
print(b.var())

c = [stats.geom.rvs(np.exp(-0.35 * lambda_), size=d).sum() for _ in range(10000)]
c = [stats.geom.rvs(np.exp(-0.15 * lambda_), size=x).sum() for x in c]
c = [stats.geom.rvs(np.exp(-0.1 * lambda_), size=x).sum() for x in c]
c = [stats.geom.rvs(np.exp(-1.0 * lambda_), size=x).sum() for x in c]
c = np.array(c)
print(c.mean())
print(c.var())

e = np.empty(1000)

for i in range(1000):
    t = np.array([0.])
    v = np.array([d])

    while t[-1] <= 1.6:
        new_t = np.min(stats.expon.rvs(scale=1/lambda_, size=v[-1]))
        t = np.append(t, t[-1]+new_t)
        v = np.append(v, v[-1]+1)
    #print(v)

    e[i] = v[-1]-1

print(e.mean())
print(e.var())


class YuleProcess():
    def __init__(self, lambda_, d):
        self.lambda_ = lambda_
        self.upsilon = stats.expon.rvs(scale=1/lambda_, size=d)

    def get_next(self):
        current_min_index = np.argmin(self.upsilon)
        current_min = self.upsilon[current_min_index]

        self.upsilon = np.append(self.upsilon, np.nan)
        self.upsilon[[current_min_index, -1]] = (current_min +
            stats.expon.rvs(scale=1/lambda_, size=2))

        return current_min

f = np.empty(1000)

for i in range(1000):
    process = YuleProcess(lambda_, d)
    t = process.get_next()
    v = d

    while t <= 10.6:
        t = process.get_next()
        v += 1
    #print(v)

    f[i] = v

print(f.mean())
print(f.var())
