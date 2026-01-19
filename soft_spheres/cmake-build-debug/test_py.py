import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import mylib as my

my.run_it('make soft_spheres.so')

import soft_spheres

def get_n(ksi, T, n=6):
	return ksi * (Temp**(3/n) * (6 / np.pi))

def proc_ksi(s, ksi, Nt, d, Temp, seed, verbose=0, t_stab=1000):
	n = get_n(ksi, Temp)
	print('n = ', n)

	E_data = soft_spheres.compute_evolution(s, n, Nt, d, Temp, seed, verbose)
	assert E_data.shape[0] == 2 * Nt, ('(len(data) = ' + str(E_data.shape[0]) + ') != (2 * Nt = ' + str(2 * Nt) + ')')
	E = E_data[:Nt]
	P = E_data[Nt : 2 * Nt]
	
	t_inds = np.arange(Nt)
	t_mem = 10
	N_uncor_steps = int(Nt / t_mem)
	stab_inds = t_inds > t_stab
	E0 = np.mean(E[stab_inds])
	P0 = np.mean(P[stab_inds])
	d_E0 = np.std(E[stab_inds]) / np.sqrt(N_uncor_steps)
	d_P0 = np.std(P[stab_inds]) / np.sqrt(N_uncor_steps)

	if(verbose):
		it = np.arange(Nt)

		fig_E, ax_E = my.get_fig(r'step', 'E', 'E(i)')
		ax_E.plot(it, E, label='E')
		ax_E.legend()

		fig_P, ax_P = my.get_fig(r'step', 'P', 'P(i)')
		ax_P.plot(it, P, label='P')
		ax_P.legend()

	return E0, P0, d_E0, d_P0

args = sys.argv[1:]
argc = len(args)

if(argc not in [5]):
	print('usage:\n' + sys.argv[0] + '   L/a   Nt   d   Temp   seed')
	exit()

s = int(args[0])
#ksi = float(args[1])
Nt = int(args[1])
d = float(args[2])
Temp = float(args[3])
my_seed = int(args[4])

test_ksi = 1.5
# 1.18
#test_ksi = -1

if(test_ksi > 0):
	proc_ksi(s, test_ksi, Nt, d, Temp, my_seed, verbose=1)
else:
	N_ksi = 10
	ksi_s = np.linspace(0.5, 1.5, N_ksi)


	E0 = np.empty(N_ksi)
	P0 = np.empty(N_ksi)
	d_E0 = np.empty(N_ksi)
	d_P0 = np.empty(N_ksi)
	for i in range(N_ksi):
		E0[i], P0[i], d_E0[i], d_P0[i] = proc_ksi(s, ksi_s[i], Nt, d, Temp, my_seed, verbose=0, t_stab=5000)
	n = get_n(ksi_s, Temp)
	z = P0 / (n * Temp)

	fig_PV, ax_PV = my.get_fig(r'$\xi$', 'z', r'z($\xi$)')
	ax_PV.plot(ksi_s, z, label='z')
	ax_PV.legend()


plt.show()

