import numpy as np
import os
import sys
import matplotlib.pyplot as plt

import mylib as my

my.run_it('make harmonic.so')

import harmonic

# pybind11::array_t<double> compute_evolution(int Nx, int Nt, double L, double d, int my_seed, int verbose)

[Nx, Nt, L, d, my_seed, verbose], _ = \
    my.parse_args(sys.argv, ['-Nx', '-Nt', '-L', '-d', '-seed', '-verbose'], \
                  possible_values=[None, None, None, None, None, None], \
                  possible_arg_numbers=[[1], [1], [1], [1], [0, 1], [0, 1]], \
                  default_values=[None, None, None, None, ['0'], ['1']])

Nx = int(Nx)
Nt = int(Nt)
L = float(L)
d = float(d)
my_seed = int(my_seed[0])
verbose = int(verbose[0])

test_mode = 1

if(test_mode == 2):
	L = np.pi/2
	x = np.linspace(-L, L, Nx + 1)
	f = np.sin((x / L + 1) / 2 * np.pi)
	dx = 2 * L / Nx;

	d2 = harmonic.comp_d2(f, dx)
	d4 = harmonic.comp_d2(d2, dx)
	d6 = harmonic.comp_d2(d4, dx)
	d8 = harmonic.comp_d2(d6, dx)

	fig, ax = my.get_fig('x', 'f')
	ax.plot(x, d4, label='d4')
	ax.plot(x, d6, label='d6')
	ax.plot(x, f, label='f')
	ax.plot(x, d2, label='d2')
	ax.legend()

	fig2, ax2 = my.get_fig('x', 'f', yscl='log')
	ax2.plot(x, np.abs(f - d4), '.', label='|f - d4|')
	ax2.plot(x, np.abs(f + d4), '.', label='|f + d6|')
	ax2.plot(x, np.abs(f - d8), '.', label='|f - d8|')
	ax2.legend()

elif(test_mode == 1):
	Ef_data = harmonic.compute_evolution(Nx, Nt, L, d, my_seed, verbose)
	E = Ef_data[:Nt]
	f = Ef_data[Nt:]
	assert len(f) == (Nx + 1), ('len(x) = ' + str(Nx + 1) + ', but len(f) = ' + str(len(f)))
	t = np.arange(Nt)
	x = np.linspace(-L, L, Nx + 1)
	
	log_f = np.log(np.abs(f) + 1e-10)
	#log_f = -np.log(f)
	a = 0.8
	x_ind = (x > -a) & (x < a)
	fit2 = np.polyfit(x[x_ind], log_f[x_ind], 2)
	sgm = 1/np.sqrt(-fit2[0])

	fig_f, ax_f = my.get_fig('x', r'$\Psi$')
	ax_f.plot(x, f, label = 'result')
	ax_f.legend()

	fig_logf, ax_logf = my.get_fig('x', r'$\Psi$', yscl='log')
	ax_logf.plot(x, f, label = 'data')
	ax_logf.plot(x[x_ind], np.exp(np.polyval(fit2, x[x_ind])), label = r'$sigma = %lf$' % (sgm))
	ax_logf.legend()

	fig_E, ax_E = my.get_fig('t (step)', r'$E$')
	ax_E.plot(t, E, label = 'evol')
	ax_E.legend()

plt.show()
