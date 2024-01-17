
# tmp

$\sigma\hat{\boldsymbol{\sigma}}$

```python
            # for i in range(self.cubic_size):
            #     for j in range(self.cubic_size):
            #         for k in range(self.cubic_size):
            #             n_particles = len(self.cells[i][j][k])
            #             for _ in range(n_particles):
            #                 particle = self.cells[i][j][k].pop()
            #                 particle[:3] += particle[3:] * dt
            #                 x, y, z = particle[:3]
            #                 # elastic reflection
            #                 if x < 0:
            #                     x, particle[3] = -x, -particle[3]
            #                 elif x > self.cubic_length:
            #                     x, particle[3] = 2 * self.cubic_length - x, -particle[3]
            #                 if y < 0:
            #                     y, particle[4] = -y, -particle[4]
            #                 elif y > self.cubic_length:
            #                     y, particle[4] = 2 * self.cubic_length - y, -particle[4]
            #                 if z < 0:
            #                     z, particle[5] = -z, -particle[5]
            #                 elif z > self.cubic_length:
            #                     z, particle[5] = 2 * self.cubic_length - z, -particle[5]
            #                 self.cells[int(x // self.box_size)][int(y // self.box_size)][int(z // self.box_size)].append(particle)
```

Create a second y-axis

```python
# # Create a second y-axis
# ax2 = plt.twinx()

# # Calculate the normalized values (percentage)
# counts_norm = counts/counts.sum()*100

# n_ticks = 10
# ax2.yaxis.set_major_locator(plt.MaxNLocator(n_ticks))

# # Calculate the normalized values for these y-ticks
# ticks = np.linspace(0, counts.max(), n_ticks)
# ticks_norm = ticks/ticks.sum()*100

# # Set the y-tick labels on the second y-axis to the normalized values
# ax2.set_yticklabels(ticks_norm)

# # Label the second y-axis
# ax2.set_ylabel('Normalized N (%)')

# set X range
# plt.xlim(0, 0.01)
```
