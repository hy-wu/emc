
# Enskog Monte Carlo simulation

Monte Carlo simulation for Enskog dense gas

## Python version

```python
from mc import *
n_particles=int(1e4)
cell_len=2
cubic_size=25
box_size=5
r=1
m=1e-26
T=300
n_steps=1000
dt=1e-8
k_B = Boltzmann * 1e9
cubic_len = cubic_size * cell_len

vs = np.concatenate((np.random.normal(loc=0, scale=np.sqrt(k_B * T / m), size=(n_particles, 2)), np.random.uniform(low=0, high=10 * np.sqrt(k_B * T / m), size=(n_particles, 1))), axis=1)
cubic = Cubic(n_particles=n_particles, cell_len=cell_len, cubic_size=cubic_size, box_size=box_size, r=r, m=m, T=T, equ=False, vs=vs)
cubic.run(n_steps=n_steps, dt=dt)
```

## using ROOT

```shell
$ root
root [0] .L nt.cpp
root [1] cubic = nt()
num_steps: 1000
sampling_step: 1
nσ^3 = 0.8
simulation[==>              ] 2% [00m:08s<06m:41s] 22/1000
```

## libTorch version

debug

```shell
$ mkdir build && cd build
$ cmake .. && make && gdb emc
(gdb) run
```
