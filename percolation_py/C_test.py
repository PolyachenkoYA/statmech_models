import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
#import copy

from scipy.optimize import fsolve
from scipy.interpolate import CubicSpline

from numba import jit

import mylib as my

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
    
def sigm_fnc_inv(f):
	return - np.log(1 / f - 1)

def sigm_fnc(x, x0, sgm):
	return 1 / (1 + np.exp(-(x - x0)/sgm))

def atan_fnc_inv(f):
	return np.tan((f - 1 / 2) * np.pi)

def atan_fnc(x, x0, sgm):
	return 1 / 2 + np.arctan((x - x0) / sgm) / np.pi
	
def tanh_fnc_inv(f):
	return np.arctanh((f - 1 / 2) * 2)

def tanh_fnc(x, x0, sgm):
	return 1 / 2 + np.tanh((x - x0) / sgm) / 2


def proc_box_size(box_size, f0, a, b, Np, N_iter, ax, lin_axes, lin_fncs, lin_lbls, fit_bounds, model_ind):
	my_seed = 0
	run_id = 0
	
	filename = os.path.join('data', r'Ni%d_Box%d_Seed%d_P%1.3lf_%1.3lf_%d_ID%d.dat' % (N_iter, box_size, my_seed, a, b, Np, run_id))
	if(not os.path.isfile(filename)):
		my.run_it('./percolation %d %lf %lf %d %d %d %d' % (N_iter, a, b, Np, box_size, my_seed, run_id))
		# %s   Ni   p1   p2   Np   BoxL   seed   id
		
	data = np.loadtxt(filename)
	p = data[:, 0]
	f = data[:, 1]
	fit_ind = (f > 1e-10) & (1-f > 1e-10)
	p_fit = p[fit_ind]
	N_data_fit = len(p_fit)
	
	N_fncs = len(lin_axes)
	lin_f = np.empty((N_data_fit, N_fncs))
	lin_fit = []
	p3_fit = []
	sgm_arr = np.empty(N_fncs)
	p_center = np.empty(N_fncs)
	p0 = np.empty(N_fncs)
	p3_sols = np.empty(N_fncs)
	n_dgt = 5
	p_draw = p_fit
	for i in range(N_fncs):
		clr = my.get_my_color(box_size)
		lin_f[:, i] = lin_fncs[i][1](f[fit_ind])
		
		fit_ind = (fit_bounds[0] < p_fit) & (p_fit < fit_bounds[1])
		lin_fit.append(np.polyfit(p_fit[fit_ind], lin_f[fit_ind, i], 1))
		p3_fit.append(np.polyfit(p_fit[fit_ind], lin_f[fit_ind, i], 3))
		sgm_arr[i] = 1 / lin_fit[i][0]
		p_center[i] = - lin_fit[i][1] * sgm_arr[i]
		lin_f0 = lin_fncs[i][1](f0)
		p0[i] = p_center[i] + sgm_arr[i] * lin_f0
		p3_sols[i] = fsolve(lambda p: (np.polyval(p3_fit[i], p) - lin_f0), p0[i])
		
		#print('pc: ', p_center[i], '; sgm: ', sgm_arr[i])
	
		lin_axes[i].plot(p_draw, np.polyval(lin_fit[i], p_draw), '--', label=None, color=clr)
		lin_axes[i].plot(p_draw, lin_f[:, i], '.', label='$L_{box} = ' + str(box_size) + '$', color=clr)
		lin_axes[i].plot(p_draw, np.polyval(p3_fit[i], p_draw), label=None, color=clr)
	#	lin_axes[i].plot([p0] * 2, [min(f), max(f)], label='$p(f = ' + my.f2s(f0, n_dgt) + ') = ' + my.f2s(p0, n_dgt) + '$', color=clr)
	
		ax.plot(p_draw, lin_fncs[i][0](p_draw, p_center[i], sgm_arr[i]), label=(lin_lbls[i] if(model_ind == 0) else None), color=clr)

	
	p_draw = p
	ax.plot(p, f, '.', label='$L_{box} = ' + str(box_size) + '$', color=my.get_my_color(box_size))
#	ax.plot([p0] * 2, [min(f), max(f)], label='$p(f = ' + my.f2s(f0, n_dgt) + ') = ' + my.f2s(p0, n_dgt) + '$', color=my.get_my_color(box_size))

	return p_center, sgm_arr, p3_sols

def main():
	my.run_it('./rebuild_py.sh ')
	
	box_sizes = np.array([1, 2, 4, 10, 20, 40, 80, 160])
	#box_sizes = [10]
	f0 = 0.5
	fit_bounds = [0.5, 0.7]
	N_big = 35

	N_iter = 5000
	a = 0.45
	b = 0.75
	Np = 30
	
	N_sizes = len(box_sizes)
	lin_fncs = [[sigm_fnc, sigm_fnc_inv], [tanh_fnc, tanh_fnc_inv], [atan_fnc, atan_fnc_inv]]
	lin_lbls_inv = ['-\ln(1/f-1)', 'arctanh(2f - 1)', 'tan((f - 1/2)\pi)']
	lin_lbls = ['sigmoid', 'tanh', 'atan']
	N_fncs = len(lin_lbls)
	N_fncs = 1
	
	#fig, ax = my.get_fig('p', 'f', xscl='log', yscl='log')
	fig, ax = my.get_fig('p', 'f', 'P(p)')
	
	lin_axes = [[]] * N_fncs
	lin_figs = [[]] * N_fncs
	for i in range(N_fncs):
		lin_figs[i], lin_axes[i] = my.get_fig('p', '$' + lin_lbls_inv[i] + '$', lin_lbls[i] + ' linearization')

	sgm_arr = np.empty((N_fncs, N_sizes))
	pc_arr = np.empty((N_fncs, N_sizes))
	p3c_arr = np.empty((N_fncs, N_sizes))
	p0_0 = np.empty(N_fncs)
	p0_log = np.empty(N_fncs)
	log_fit = np.empty((2, N_fncs))
	log_fit_best = np.empty((2, N_fncs))
	#R_log = np.empty(N_fncs) 
	nu = np.empty(N_fncs)
	nu_best = np.empty(N_fncs)
	log_optimize = []
	for i in range(N_sizes):
		pc_arr[:, i], sgm_arr[:, i], p3c_arr[:, i] = proc_box_size(box_sizes[i], f0, a, b, Np, N_iter, ax, lin_axes, lin_fncs, lin_lbls, fit_bounds, i)
	big_box_inds = (box_sizes > N_big)
	for i in range(N_fncs):
		p0_0[i] = np.mean(p3c_arr[i, big_box_inds])
		
		#small_boxes = box_sizes[~big_box_inds]
		#dp_0 = p0_0[i] - p3c_arr[i, ~big_box_inds]
		log_fit[:, i] = np.polyfit(np.log(box_sizes[~big_box_inds]), np.log(p0_0[i] - p3c_arr[i, ~big_box_inds]), 1)
		#R_log[i], _ = scipy.stats.pearsonr(np.log(small_boxes), np.log(np.abs(dp_0)))
		nu[i] = -1 / log_fit[0, i]
		
		min_R_fnc = lambda c: \
				scipy.stats.pearsonr( \
					np.log(box_sizes[~big_box_inds]), \
					np.log(np.abs(c - p3c_arr[i, ~big_box_inds])))[0]
		
		log_optimize.append(scipy.optimize.minimize(min_R_fnc, p0_0[i]))
		p0_log[i] = log_optimize[i].x[0]
		log_fit_best[:, i] = np.polyfit(np.log(box_sizes[~big_box_inds]), np.log(p0_log[i] - p3c_arr[i, ~big_box_inds]), 1)
		nu_best[i] = -1 / log_fit_best[0, i]
	
	ax.plot([a, b], [f0] * 2, label='$f_0 = ' + my.f2s(f0) + '$')
	ax.legend()
	
	fig_sgm, ax_sgm = my.get_fig('$L$', '$\sigma$', xscl='log', yscl='log')
	sgm_fit = np.empty((N_fncs, 2))
	for i in range(N_fncs):
		lin_axes[i].plot([a, b], [lin_fncs[i][1](f0)] * 2, label='$f_0 = ' + my.f2s(f0) + '$', color=my.get_my_color(i))
		lin_axes[i].legend()
		
		sgm_fit[i, :] = np.polyfit(np.log(box_sizes), np.log(sgm_arr[i, :]), 1)
		ax_sgm.plot(box_sizes, sgm_arr[i, :], 'o', label=lin_lbls[i], color=my.get_my_color(i))
		ax_sgm.plot(box_sizes, np.exp(np.polyval(sgm_fit[i, :], np.log(box_sizes))), label='k = ' + my.f2s(sgm_fit[i, 0]), color=my.get_my_color(i))
	ax_sgm.legend()
	
	fig_p0, ax_p0 = my.get_fig('$L$', '$p_0$', title='$p_0(f_0 = ' + str(f0) + ', L)$', xscl='log')
	for i in range(N_fncs):
		ax_p0.plot(box_sizes, p3c_arr[i, :], 'o', label=lin_lbls[i], color=my.get_my_color(i))
		ax_p0.plot([N_big, max(box_sizes)], [p0_0[i]] * 2, label='$p_{lim} = ' + my.f2s(p0_0[i], 4) + '$', color=my.get_my_color(i))
	ax_p0.legend()

	fig_p0_N, ax_p0_N = my.get_fig('$L$', '$p_{lim} - p_0(L)$', title='$p_0(f_0 = ' + str(f0) + ', L)$', xscl='log', yscl='log')
	for i in range(N_fncs):
		ax_p0_N.plot(box_sizes[~big_box_inds], p0_0[i] - p3c_arr[i, ~big_box_inds], 'o', label=lin_lbls[i], color=my.get_my_color(i))
		ax_p0_N.plot(box_sizes[~big_box_inds], np.exp(np.polyval(log_fit[:, i], np.log(box_sizes[~big_box_inds]))), \
		             label=r'$\nu = ' + my.f2s(nu[i]) + '$', color=my.get_my_color(i))
		
		ax_p0_N.plot(box_sizes[~big_box_inds], p0_log[i] - p3c_arr[i, ~big_box_inds], '+', label='$p_{optimize} = ' + my.f2s(p0_log[i]) + '$', color=my.get_my_color(i + N_fncs))
		ax_p0_N.plot(box_sizes[~big_box_inds], np.exp(np.polyval(log_fit_best[:, i], np.log(box_sizes[~big_box_inds]))), \
		             label=r'$\nu_{optimize} = ' + my.f2s(nu_best[i]) + '$', color=my.get_my_color(i + N_fncs))
		
	ax_p0_N.legend()

	fig_p0_inv, ax_p0_inv = my.get_fig('$1/L$', '$p_0$', title='$p_0(f_0 = ' + str(f0) + ', L)$')
	for i in range(N_fncs):
		ax_p0_inv.plot(1/box_sizes, p3c_arr[i, :], 'o', label=lin_lbls[i], color=my.get_my_color(i))
	ax_p0_inv.legend()

	plt.show()

if __name__ == "__main__":
	main()

