import numpy as np
import matplotlib.pyplot as plt
import sys

import mylib as my

my.run_it('cmake ..')
my.run_it('make ')

import test as izing

args = sys.argv[1:]
argc = len(args)
if(argc not in [4]):
	print('usage:\n' + sys.argv[0] + '   N   T   Nt   seed')
	exit()

N = int(args[0])
Temp = float(args[1])
Nt = int(args[2])
my_seed = int(args[3])

E = izing.compute_evolution(N, Temp, Nt, my_seed)
steps = np.arange(Nt)

fig, ax = my.get_fig('step', 'E')
ax.plot(steps, E)

Nn = 10
a = np.arange(Nn)
b = np.arange(Nn) * 2

c = izing.add_arrays(a,b)

print(a)
print(b)
print(c)

my.run_it(r'./Izing %d %lf %d %d' % (N, T, Nt, my_seed), verbose=1)

plt.show()
