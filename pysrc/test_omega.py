from mc import *
from matplotlib import pyplot as plt
import numpy as np
import time

n_particles=1e5
cell_size=2
cubic_size=50
box_size=5
r=1
m=1e-26
T=300
n_steps=1e2
dt=1e-9

time_start = time.time()
cubic = Cubic(n_particles=n_particles, cell_len=cell_size, cubic_size=cubic_size, r=r, m=m, T=T)

cubic.run(n_steps=n_steps, dt=dt)
print(max(cubic.ωs))

fig, axs = plt.subplots(2)
counts, bins, patches = axs[0].hist(cubic.ωs, bins=1000)
axs[0].set_xlabel('ω')
axs[0].set_ylabel('N')
axs[0].set_title('ω distribution: n_particles = {}, cell_size = {}, cubic_size = {},\n r = {}, m = {}, T = {}, n_steps = {}, dt = {}, time = {}s'.format(n_particles, cell_size, cubic_size, r, m, T, n_steps, dt, round(time.time() - time_start)))

# Store the mean probabilities at each step
mean_probs = []

total_trials = 0
accepted_trials = 0

for ω in cubic.ωs:
    # Generate a random number
    rand_num = np.random.random()

    # If the random number is less than ω, increment the number of accepted trials
    if rand_num < ω:
        accepted_trials += 1

    # Increment the total number of trials
    total_trials += 1

    # Calculate the mean probability of acceptance and store it
    mean_prob = accepted_trials / total_trials
    mean_probs.append(mean_prob)

# Plot the mean probabilities
axs[1].plot(mean_probs)
axs[1].set_xlabel('Step number')
axs[1].set_ylabel('Mean probability of acceptance')

# Adjust the space between the subplots
plt.tight_layout()
plt.savefig('omega.pdf')
plt.show()
