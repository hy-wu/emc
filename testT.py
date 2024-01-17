from mc import *
from matplotlib import pyplot as plt
import numpy as np
import time

n_particles=1e3
cell_size=2
cubic_size=10
r=1
m=1e-26
T=0.01
n_steps=1e3
dt=1e-7

time_start = time.time()
cubic = Cubic(n_particles=n_particles, cell_len=cell_size, cubic_size=cubic_size, r=r, m=m, T=T)
cubic.run(n_steps=n_steps, dt=dt)
# print(max(cubic.ωs))

times = np.arange(0, n_steps * dt - dt/2, dt)
print('T_start =', cubic.Ts[0], 'T_end =', cubic.Ts[-1], 'T_mean =', np.mean(cubic.Ts), 'T_max =', max(cubic.Ts), 'T_min =', min(cubic.Ts))

plt.plot(times, cubic.Ts)
plt.title('Temperature: n_particles = {}, cell_size = {}, cubic_size = {},\n r = {}, m = {}, T = {}, n_steps = {}, dt = {}, time = {}s'.format(n_particles, cell_size, cubic_size, r, m, T, n_steps, dt, round(time.time() - time_start)))
plt.savefig('T.pdf')
