import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.animation as animation
n_steps = 1000
n_particles = 100000
fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot(111)

def update(num):
    ax.clear()
    data = pd.read_csv('particle_records.csv', header=None, names=['step','i','x', 'y', 'z', 'vx', 'vy', 'vz'], skiprows=num*n_particles*10, nrows=n_particles)
    ax.hist(data['x'][:10000], bins=100)

ani = animation.FuncAnimation(fig, update, frames=n_steps // 10, repeat=False)

plt.show()