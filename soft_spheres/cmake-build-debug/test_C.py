import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import mylib as my

my.run_it('make soft_spheres')

def get_n(ksi, T, n=6):
	return ksi * (Temp**(3/n) * (6 / np.pi))

def proc_ksi(s, ksi, Nt, d, Temp, seed, verbose=0):
	n = get_n(ksi, Temp)
	print('n = ', n)

	my.run_it(r'./soft_spheres %d %lf %d %lf %lf %d' % (s, n, Nt, d, Temp, my_seed), verbose=1)

	filename = r's%d_n%lf_Nt%d_d%lf_Temp%lf_seed%d.dat' % (s, n, Nt, d, Temp, my_seed)

	E_data = np.loadtxt(filename)
	assert E_data.shape[0] == Nt, ('(len(data) = ' + str(E_data.shape[0]) + ') != (Nt = ' + str(Nt) + ')')
	E = E_data[:, 0]
	P = E_data[:, 1]

	if(verbose):
		it = np.arange(Nt)
		
		fig_E, ax_E = my.get_fig(r'step', 'E', 'E(i)')
		ax_E.plot(it, E, label='E')
		ax_E.legend()

		fig_P, ax_P = my.get_fig(r'step', 'P', 'P(i)')
		ax_P.plot(it, P, label='P')
		ax_P.legend()
		
	return E, P

args = sys.argv[1:]
argc = len(args)

if(argc not in [5]):
	#print('usage:\n' + sys.argv[0] + '   L/a   ksi   Nt   d   Temp   seed')
	print('usage:\n' + sys.argv[0] + '   L/a   Nt   d   Temp   seed')
	exit()

s = int(args[0])
#ksi = float(args[1])
Nt = int(args[1])
d = float(args[2])
Temp = float(args[3])
my_seed = int(args[4])

test_ksi = -1
# 1.18

if(test_ksi > 0):
	proc_ksi(s, test_ksi, Nt, d, Temp, my_seed, verbose=1)
else:
	N_ksi = 10
	ksi_s = np.linspace(0.5, 1.5, N_ksi)
	
	E = np.empty((Nt, N_ksi))
	P = np.empty((Nt, N_ksi))
	for i in range(N_ksi):
		E[:, i], P[:, i] = proc_ksi(s, ksi_s[i], Nt, d, Temp, my_seed, verbose=0)
	
	
	t_stab = 5000
	t_inds = np.arange(Nt)
	stab_inds = t_inds > t_stab
	E0 = np.empty(N_ksi)
	P0 = np.empty(N_ksi)
	for i in range(N_ksi):
		E0[i] = np.mean(E[stab_inds, i])
		P0[i] = np.mean(P[stab_inds, i])
	n = get_n(ksi_s, Temp)
	z = P0 / (n * Temp)
	
	fig_PV, ax_PV = my.get_fig(r'$\xi$', 'z', r'z($\xi$)')
	ax_PV.plot(ksi_s, z, label='z')
	ax_PV.legend()


plt.show()
