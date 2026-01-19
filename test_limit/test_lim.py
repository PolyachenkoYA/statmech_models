import numpy as np
import matplotlib.pyplot as plt
import mylib as my

def gauss(x0, s, x):
	return np.exp(-np.power((x-x0) / s, 2) / 2) / np.sqrt(2 * np.pi * s**2)

def plot_dist(N, M):
	data = np.random.rand(N, M)
	h = np.mean(data, axis=0)

	plt.hist(h, bins=30, density=True, alpha=0.2, label = 'N = ' + str(N), color=my.get_my_color(n))
	
	s = np.sqrt(1 / 12 / N)
	x_draw = np.linspace(x0 - 5 * s, x0 + 5 * s, 100)
	plt.plot(x_draw, gauss(x0, s, x_draw), color=my.get_my_color(n))


M = 100000
x0 = 0.5

plt.figure

for n in range(2, 5):
	plot_dist(n, M)


plt.legend()
plt.show()

