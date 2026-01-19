import numpy as np
import matplotlib.pyplot as plt

N = 40
E0 = 10
d0 = 0.5
Ed0 = 0.5
Nt = 100000
Nt_stab = 1000

v0 = np.sqrt(2 * E0) / N

E = np.empty((Nt + 1, 1))
v = np.empty((N, Nt + 1))

#rnd = np.random.rand(Nt, 1)
E[0] = np.sum(np.power(v, 2)) / 2
v[:, 0] = np.squeeze(np.ones((N, 1))) * v0
Ed = Ed0
for ti in range(Nt):
	d = (np.random.rand() * 2 - 1) * d0
	pi = np.random.randint(N)
	v[:, ti + 1] = v[:, ti]
	
	v[pi, ti + 1] = v[pi, ti + 1] + d
	Enew = np.sum(np.power(v[:, ti + 1], 2)) / 2
	Ed_new = Ed - (Enew - E[ti])
	
	if(Ed_new > 0):
		E[ti + 1] = Enew
		Ed = Ed_new
	else:
		E[ti + 1] = E[ti]
		v[pi, ti + 1] = v[pi, ti + 1] - d
		
	print((ti + 1) / (Nt))

plt.figure()
plt.plot(range(Nt + 1), E)
plt.xlabel('t')
plt.ylabel('E')

plt.figure()
plt.hist(v[:, Nt_stab : ].flatten(), bins=50, density=True)
plt.xlabel('v')
plt.ylabel('hist')

plt.figure()
plt.hist(np.power(v[:, Nt_stab : ].flatten(), 2), bins=50, density=True)
plt.xlabel('v^2')
plt.ylabel('hist')

plt.show()
