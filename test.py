# from mc import *
# from matplotlib import pyplot as plt

# cubic = Cubic(n_particles=10000, cell_size=0.5, cubic_size=10, r=0.1, m=1e-21, T=300)

# cubic.run(n_steps=10, dt=1e-14)

# plt.scatter(cubic.particles()[:, 0], cubic.particles()[:, 1], s=1)
# plt.show()

from mc import *
from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

cubic = Cubic(n_particles=10000, cell_len=0.2, cubic_size=25, r=0.1, m=1e-21, T=300)

cubic.run(n_steps=100, dt=1e-3)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(cubic.particles()[:, 0], cubic.particles()[:, 1], cubic.particles()[:, 2], s=1)

plt.show()

# from mc import *
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# # from mpl_toolkits.mplot3d import Axes3D

# cubic = Cubic(n_particles=10000, cell_size=0.5, cubic_size=10, r=0.1, m=1e-21, T=300)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# def update(num):
#     ax.clear()
#     cubic.run(n_steps=1, dt=1e-3, quiet=True)
#     ax.scatter(cubic.particles()[:, 0], cubic.particles()[:, 1], cubic.particles()[:, 2], s=1)

# ani = animation.FuncAnimation(fig, update, frames=100, repeat=False)

# plt.show()