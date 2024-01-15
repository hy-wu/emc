from mc import *
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time

n_particles=1e4
cell_size=2
cubic_size=25
box_size=5
r=1
m=1e-26
T=300
# n_steps=4
dt=1e-9

time_start = time.time()
cubic = Cubic(n_particles=n_particles, cell_len=cell_size, cubic_size=cubic_size, box_size=box_size, r=r, m=m, T=T)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.zeros((cubic.cubic_box_size, cubic.cubic_box_size, cubic.cubic_box_size))
y = np.zeros((cubic.cubic_box_size, cubic.cubic_box_size, cubic.cubic_box_size))
z = np.zeros((cubic.cubic_box_size, cubic.cubic_box_size, cubic.cubic_box_size))

for i in range(cubic.cubic_box_size):
    for j in range(cubic.cubic_box_size):
        for k in range(cubic.cubic_box_size):
            x[i, j, k] = cubic.box_len * i
            y[i, j, k] = cubic.box_len * j
            z[i, j, k] = cubic.box_len * k

x = x.flatten()
y = y.flatten()
z = z.flatten()

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

def update(num):
    ax.clear()
    cubic.run(n_steps=1, dt=1e-9, quiet=True)
    qk_normed = cubic.qk / np.linalg.norm(cubic.qk, axis=3, keepdims=True)

    u = qk_normed[:, :, :, 0].flatten()
    v = qk_normed[:, :, :, 1].flatten()
    w = qk_normed[:, :, :, 2].flatten()
    # u = qk[:, :, :, 0].flatten()
    # v = qk[:, :, :, 1].flatten()
    # w = qk[:, :, :, 2].flatten()

    ax.quiver(x, y, z, u, v, w, length=5, color='red')
    ax.set_xlim([0, cubic.cubic_len])
    ax.set_ylim([0, cubic.cubic_len])
    ax.set_zlim([0, cubic.cubic_len])

ani = animation.FuncAnimation(fig, update, frames=100, repeat=False)

plt.savefig('quiver.pdf')
plt.show()