import numpy as np
import os
import sys
import matplotlib.pyplot as plt

import mylib as my

my.run_it('make polymer.so')

import polymer

# pybind11::array_t<double> compute_evolution(int Nx, int Nt, int my_seed, int verbose)

[Nx, Nt, my_seed, verbose, test_mode], _ = \
    my.parse_args(sys.argv, ['-Nx', '-Nt', '-seed', '-verbose', '-test_mode'], \
                  possible_values=[None, None, None, None, None], \
                  possible_arg_numbers=[[1], [1], [0, 1], [0, 1], [0, 1]], \
                  default_values=[None, None, ['0'], ['0'], ['2']])

Nx = int(Nx)
Nt = int(Nt)
my_seed = int(my_seed[0])
verbose = int(verbose[0])
test_mode = int(test_mode[0])

if(test_mode == 1):
	R_data = polymer.compute_evolution(Nx, Nt, my_seed, verbose)

	plt.hist(R_data, bins=50, density=True)

elif(test_mode == 2):
	Ls = np.array([3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20])
	#Ls = np.array([5, 10])
	
	N_L = len(Ls)
	R_means = np.empty(N_L)
	d_R_means = np.empty(N_L)
	for i in range(N_L):
		R_data = polymer.compute_evolution(Ls[i], Nt, my_seed, 0)
		R_means[i] = np.mean(R_data)
		d_R_means[i] = np.std(R_data)
	
	fig, ax = my.get_fig('L', '<R>')
	
	ax.errorbar(Ls, R_means, yerr=d_R_means, fmt='.')

plt.show()
