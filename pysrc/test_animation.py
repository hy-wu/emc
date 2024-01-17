from mc import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# from mpl_toolkits.mplot3d import Axes3D

cubic = Cubic(n_particles=1e4, cell_len=2, cubic_size=25, box_size=5, r=1, m=1e-26, T=300)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def update(num):
    ax.clear()
    cubic.run(n_steps=1, dt=1e-8, quiet=True)
    ax.scatter(cubic.particles()[:, 0], cubic.particles()[:, 1], cubic.particles()[:, 2], s=1)
    ax.set_xlim([0, cubic.cubic_len])
    ax.set_ylim([0, cubic.cubic_len])
    ax.set_zlim([0, cubic.cubic_len])

ani = animation.FuncAnimation(fig, update, frames=100, repeat=False)

plt.show()
