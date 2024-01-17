from mc import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

cubic = Cubic(n_particles=1e1, cell_len=0.2, cubic_size=100, r=0.1, m=1e-21, T=3000)

# Initialize a list to hold the particle positions at each time step
particle_positions = []

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def update(num):
    ax.clear()
    cubic.run(n_steps=1, dt=1e-2, quiet=True)
    
    # Save the current particle positions
    particle_positions.append(cubic.particles().copy())
    
    # Draw the particle trajectories
    for i in range(len(cubic.particles())):
        positions = np.array([pos[i] for pos in particle_positions])
        ax.plot(positions[:, 0], positions[:, 1], positions[:, 2])
    
    ax.set_xlim([0, cubic.cubic_len])
    ax.set_ylim([0, cubic.cubic_len])
    ax.set_zlim([0, cubic.cubic_len])

ani = animation.FuncAnimation(fig, update, frames=1000, repeat=False)

plt.show()