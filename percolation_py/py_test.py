#import os
import numpy as np
import matplotlib.pyplot as plt
#import copy

from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline

from numba import jit

import mylib as my

@jit
def prob_for_p(p, box_size, N_iter, my_seed=None, verbose=0):
	if(my_seed is None):
		my_seed = np.random.randint(N_iter * 100, dtype=np.uint32)

	np.random.seed(my_seed)

	box_size2 = box_size * box_size

	grid_v = np.empty(box_size2, dtype=np.intc)
	is_checked_v = np.empty(box_size2, dtype=np.intc)
	cluster_v = np.empty(box_size2, dtype=np.intc)
	results = np.empty(N_iter)
	for i in range(N_iter):
		fill_grid_v(grid_v, box_size, p)
		results[i] = is_infinite_grid_v(grid_v, box_size, is_checked_v, cluster_v)

	return np.mean(results)

@jit
def is_infinite_grid_v(grid, N, is_checked, cluster):
	N2 = N * N

	is_checked[:] = 0

	i = 0
	cluster_size = N2   # initial N2 so it's all initialized to -1
	while(i < N):   # we need to check only the top row to find infinite clusters
		if(grid[i]):
			clear_cluster_v(cluster, cluster_size)
			cluster_size = add_to_cluster_v(grid, N, is_checked, cluster, 0, i)
			if(is_infinite_cluster_v(cluster, cluster_size, N)):
				return 1
		else:
			is_checked[i] = 1

		i = i + 1

	return 0

@jit
def add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos):
	if(not is_checked[pos]):
		#print('hi', cluster_size)
		is_checked[pos] = 1
		if(grid[pos]):
			N2 = N * N
			cluster[cluster_size] = pos
			cluster_size = cluster_size + 1
			#for i in range(cluster_size):
				#print(cluster[i], end=' ')

			if(pos - N >= 0):
				cluster_size = add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos - N)
			if((pos - 1) / N == (pos / N)):
				cluster_size = add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos - 1)
			if((pos + 1) / N == (pos / N)):
				cluster_size = add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos + 1)
			if(pos + N < N2):
				cluster_size = add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos + N)
	return cluster_size

@jit
def is_infinite_cluster_v(cluster, cluster_size, N):
	if(cluster[0] >= N):
		return 0

	i = 0;
	next_row_to_add = 0;
	while(i < cluster_size):
		if(np.floor(cluster[i] / N) == next_row_to_add):
			next_row_to_add = next_row_to_add + 1
			if(next_row_to_add == N):
				return 1
		i = i + 1

	return 0

@jit
def fill_grid_v(grid, N, p):
	N2 = N * N;

	for i in range(N2):
		grid[i] = (np.random.rand() < p)

	return 0

@jit
def clear_cluster_v(cluster, cluster_size):
	for i in range(cluster_size):
		cluster[i] = -1;

	return 0;

def interpolate(xs, ys, a, b):
    #xs = np.linspace(a,b, numPoints, endpoint=True)
    spline = CubicSpline(xs, ys, extrapolate=False)

    def my_spline(x):
        if x<a:
            return spline(a)
        elif x>b:
            return spline(b)
        else:
            return spline(x)

    return np.vectorize(my_spline)

def main():
	#my.run_it('./rebuild_py.sh ')
	#import percol

	box_size = 20
#	my_seed = 0
	f0 = 0.5

	N_iter = 1000
	a = 0.6
	p_init = 0.8
	b = 0.9
	Np = 30

	N_iter = 100000
	a = 0.7
	p_init = 0.75
	b = 0.8
	Np = 30

	N_iter = 1000
	a = 0.4
	p_init = 0.8
	b = 0.95
	Np = 100

	#p_min_log = -3
	#probs = np.exp(np.linspace(p_min_log, 0, Np, endpoint=False))
	probs = np.linspace(a, b, Np)

	f = np.empty(Np)
	for i, p in enumerate(probs):
		#f[i] = percol.prob_for_p(p, box_size, N_iter, my_seed)
		f[i] = prob_for_p(p, box_size, N_iter)

	interp_spline = interpolate(probs, f - f0, a, b)
	fit3 = np.polyfit(probs, f - f0, 3);
	p0 = fsolve(lambda p : np.polyval(fit3, p), p_init)
	#p0 = fsolve(lambda p : (prob_for_p(p, box_size, N_iter) - f0), p_init)
	n_dgt = 5

	#fig, ax = my.get_fig('p', 'f', xscl='log', yscl='log')
	fig, ax = my.get_fig('p', 'f')
	p_draw = np.linspace(a, b, 1000)
	ax.plot(probs, f, '.', label='data')
	#ax.plot(p_draw, np.polyval(fit3, p_draw) + f0, label='cube interp')
	ax.plot([min(probs), max(probs)], [f0] * 2, label='$f_0 = ' + my.f2s(f0) + '$')
	#ax.plot([p0] * 2, [min(f), max(f)], label='$p(f = ' + my.f2s(prob_for_p(p0, box_size, N_iter * Np), n_dgt) + ') = ' + my.f2s(p0, n_dgt) + '$')
	ax.legend()
	plt.show()

if __name__ == "__main__":
	main()

