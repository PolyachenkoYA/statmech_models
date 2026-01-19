import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import mylib as my

args = sys.argv[1:]
argc = len(args)
if(argc not in [4]):
	print('usage:\n' + sys.argv[0] + '   N   T   Nt   seed')
	exit()

N = int(args[0])
T = float(args[1])
Nt = int(args[2])
my_seed = int(args[3])

my.run_it('make Izing')

my.run_it(r'./Izing %d %lf %d %d' % (N, T, Nt, my_seed), verbose=1)

filename = r'N%d_T%lf_Nt%d_seed%d.dat' % (N, T, Nt, my_seed)

E_data = np.loadtxt(filename)
assert len(E_data) == Nt, ('(len(data) = ' + str(len(E_data)) + ') != (Nt = ' + str(Nt) + ')')
it = np.arange(Nt)

#print(data)

fig, ax = my.get_fig(r'step', 'E', 'E(i)')
ax.plot(it, E_data)

plt.show()
