import numpy as np
import matplotlib.pyplot as plt
import os

import mylib as my

my.run_it('make izing.so')

import izing

def proc_T(Temp, box, Nt, my_seed, verbose=0, equil_time=None):
	E = izing.compute_evolution(box, Temp, Nt, my_seed)
	assert len(E) == 2*Nt, 'returned EM.shape = ' + str(E.shape) + ', but 2*Nt = ' + str(2*Nt)
	M = E[Nt : 2*Nt] / box**2
	E = E[0  :   Nt] / box**2
	
	steps = np.arange(Nt)
	if(equil_time is None):
		equil_time = np.exp(2 / Temp) * 3
	stab_step = int(min(box**2 * equil_time, Nt / 2))
	stab_ind = (steps > stab_step)
	
	E_stab = E[stab_ind]
	E_mean = np.mean(E_stab)
	E_std = np.std(E_stab)

	M_stab = M[stab_ind]
	M_mean = np.mean(M_stab)
	M_std = np.std(M_stab)
	
	C = E_std**2 / Temp**2 * box**2

	dE_step = 2
	memory_time = E_std / dE_step * np.exp(dE_step / Temp)
	d_E_mean = E_std / np.sqrt(Nt / memory_time)
	d_M_mean = M_std / np.sqrt(Nt / memory_time)
	
	if(verbose):
		fig_E, ax_E = my.get_fig('step', 'E / N')
		ax_E.plot(steps, E, label='data')
		ax_E.plot([stab_step] * 2, [min(E), max(E)], label='equilibr')
		ax_E.plot([stab_step, Nt], [E_mean] * 2, label=('$E = ' + my.errorbar_str(E_mean, d_E_mean) + '$'))
		ax_E.legend()
		
		fig_M, ax_M = my.get_fig('step', 'M / N')
		ax_M.plot(steps, M, label='data')
		ax_M.plot([stab_step] * 2, [min(M), max(M)], label='equilibr')
		ax_M.plot([stab_step, Nt], [M_mean] * 2, label=('$M = ' + my.errorbar_str(M_mean, d_M_mean) + '$'))
		ax_M.legend()
		
	return E_mean, d_E_mean, M_mean, d_M_mean, C

def get_E_T(box, Nt, my_seed, T_arr, time_verb=0, E_verb=0):
	N_T = len(T_arr)
	E_arr = np.empty(N_T)
	d_E_arr = np.empty(N_T)
	M_arr = np.empty(N_T)
	d_M_arr = np.empty(N_T)
	C_fluct = np.empty(N_T)
	for i in range(N_T):
		E_arr[i], d_E_arr[i], M_arr[i], d_M_arr[i], C_fluct[i] = proc_T(T_arr[i], box, Nt, my_seed, verbose=time_verb)
		print((i + 1) / N_T)
	
	return E_arr, d_E_arr, M_arr, d_M_arr, C_fluct

def get_deriv(x, y, n_gap=2, degr=1):
	N = len(x)
	
	deriv = []
	x_der = []
	#for i in range(n_gap, N - n_gap - 2):
	for i in range(N):
		fit_ind = np.arange(max(0, i - n_gap), min(N - 1, i + n_gap))
		fit = np.polyfit(x[fit_ind] - x[i], y[fit_ind], degr)
		deriv.append(fit[degr - 1])
		x_der.append(x[i])
	
	return np.array(x_der), np.array(deriv)
	
def proc_N(box, Nt, my_seed, ax_C=None, ax_E=None, ax_M=None, ax_M2=None, recomp=0, T_min_log10=-0.1, T_max_log10=1.5, N_T=100):
	N_particles = box**2
		
	T_arr = np.power(10, np.linspace(T_min_log10, T_max_log10, N_T))
	#T_arr = [1, 10]
	res_filename = r'N%d_Nstep%d_Tmin%lf_Tmax%lf_NT%d_ID%d.npz' % (box, Nt, T_min_log10, T_max_log10, N_T, my_seed)
	
	if(os.path.isfile(res_filename) and (not recomp)):
		print('loading ' + res_filename)
		E_data = np.load(res_filename)
		E_arr = E_data['E_arr']
		d_E_arr = E_data['d_E_arr']
		M_arr = E_data['M_arr']
		d_M_arr = E_data['d_M_arr']
		C_fluct = E_data['C_fluct']
	else:
		print('computing E & M')
		E_arr, d_E_arr, M_arr, d_M_arr, C_fluct = get_E_T(box, Nt, my_seed, T_arr)
		np.savez(res_filename, E_arr=E_arr, d_E_arr=d_E_arr, M_arr=M_arr, d_M_arr=d_M_arr, C_fluct=C_fluct)
		print('saved ' + res_filename)
		
	if(ax_E is not None):
		#fig_E, ax_E = my.get_fig('T', 'E', xscl='log')
		ax_E.errorbar(T_arr, E_arr, yerr=d_E_arr, label=('box = ' + str(box)))
	if(ax_M is not None):
		ax_M.errorbar(T_arr, M_arr, yerr=d_E_arr, label=('box = ' + str(box)))
	
	if(ax_M2 is not None):
		T_small_ind = T_arr < 3.3
		ax_M2.errorbar(T_arr[T_small_ind], np.power(M_arr[T_small_ind], 2), yerr=d_E_arr[T_small_ind], label=('box = ' + str(box)))

	T_C_1, C_1 = get_deriv(T_arr, E_arr, degr=1, n_gap=3)
	#T_C_3, C_3 = get_deriv(T_arr, E_arr, degr=3)
	
	if(ax_C is not None):
		#fig_C, ax_C = my.get_fig('T', 'C', yscl='log')
		ax_C.plot(T_arr, C_fluct * box**2, label=('box = ' + str(box)))
		#ax_C.plot(T_C_1, C_1, label=('box = ' + str(box)))
		#ax_C.plot(T_C_3, C_3, label='d=3')
		#ax_C.legend()
		
	return T_arr

import Izing_data

T_table = 1 / Izing_data.data[:, 0]
E_table = Izing_data.data[:, 1]

Nt = 10000000
my_seed = 0
recomp = 0
test_C = 0
T_min_log10 = -0.1
T_max_log10 = 1.5
N_T = 100

N_s = [2, 3, 4, 5, 10, 100, 1000]
N_s = [3, 4, 5, 10, 100]
#N_s = [100]

if(test_C):
	Temp = 2.0
	proc_T(Temp, N_s[0], Nt, my_seed, verbose=1)
else:
	fig_E, ax_E = my.get_fig('T', '$E/N$', xscl='log')
	fig_M, ax_M = my.get_fig('T', '$M/N$', xscl='log')
	fig_C, ax_C = my.get_fig('T', '$C/N$', xscl='log', yscl='log')
	fig_M2, ax_M2 = my.get_fig('T', '$(M/N)^2$')
	
	T_arr = np.power(10, np.linspace(T_min_log10, T_max_log10, N_T))
	
	for n in N_s:
		T_arr = proc_N(n, Nt, my_seed, ax_C=ax_C, ax_E=ax_E, ax_M=ax_M, ax_M2=ax_M2, recomp=recomp, T_min_log10=T_min_log10, T_max_log10=T_max_log10, N_T=N_T)
	
	ax_E.plot(T_arr, -np.tanh(1/T_arr)*2, '--', label='$-2 th(1/T)$')
	T_simulated_ind = (min(T_arr) < T_table) & (T_table < max(T_arr))
	ax_E.plot(T_table[T_simulated_ind], E_table[T_simulated_ind], '.', label='Onsager')

	ax_C.legend()
	ax_E.legend()
	ax_M.legend()

plt.show()
