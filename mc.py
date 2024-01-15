import numpy as np
from tqdm import tqdm
from scipy.constants import Boltzmann

class Cubic:
    def __init__(self, n_particles, cell_len, cubic_size, box_size, r, m, T, equ=True, vs=None):
        # unit: nm, kg, K
        self.k_B = Boltzmann * 1e9
        self.n_particles = n_particles = int(n_particles)
        self.cell_len = cell_len
        self.box_size = box_size
        self.cell_volume = cell_len ** 3
        self.cubic_size = cubic_size = int(cubic_size)
        assert cubic_size % box_size == 0, 'cubic_size % box_size != 0'
        self.cubic_box_size = cubic_size // box_size
        self.cubic_len = cubic_size * cell_len
        self.box_len = box_size * cell_len
        self.r = r
        self.m = m
        self.T0 = T
        self.t = 0
        self.ts = []
        self.Ts = []
        
        # [x, y, z, vx, vy, vz]
        particles = np.concatenate((np.random.uniform(low=0, high=self.cubic_len, size=(n_particles, 3)), np.random.normal(loc=0, scale=np.sqrt(self.k_B * T / m), size=(n_particles, 3))), axis=1)
        if not equ:
            particles[:, 3:] = vs
        self.cells = [[[[] for _ in range(cubic_size)] for _ in range(cubic_size)] for _ in range(cubic_size)]
        for particle in particles:
            x, y, z = particle[:3]
            self.cells[int(x // self.cell_len)][int(y // self.cell_len)][int(z // self.cell_len)].append(particle)
        print('nσ^3 =', n_particles / cubic_size ** 3 * (2 * r / cell_len) ** 3)
        self.ωs = []
        self.us = []
        self.pks = []
        self.qks = []
    
    def run(self, n_steps, dt, quiet=False):
        # dt to big (dt > 1e-1) may cause error
        n_steps = int(n_steps)
        particles = []
        def refresh_particles():
            particles.clear()
            for i in range(self.cubic_size):
                    for j in range(self.cubic_size):
                        for k in range(self.cubic_size):
                            particles.extend(self.cells[i][j][k])
        def bound(x):
            if x < 0:
                return -x
            elif x > self.cubic_len:
                return 2 * self.cubic_len - x
            else:
                return x
        def cycle_bound(x):
            if x < 0:
                return self.cubic_len + x
            elif x > self.cubic_len:
                return x - self.cubic_len
            else:
                return x
        refresh_particles()
        # print(len(particles))
        pbar = tqdm(total=n_steps, disable=quiet)
        pbar.set_description(f"{len(particles)} particles")
        for _ in range(n_steps):
            self.t += dt
            self.ts.append(self.t)
            # free movement
            refresh_particles()
            # pbar.set_description(f"{len(particles)} particles")
            if not quiet:
                pbar.update(1)
            self.cells = [[[[] for _ in range(self.cubic_size)] for _ in range(self.cubic_size)] for _ in range(self.cubic_size)]
            for particle in particles:
                particle[:3] += particle[3:] * dt
                x, y, z = particle[:3]
                # elastic reflection
                if x < 0:
                    x, particle[3] = -x, -particle[3]
                elif x > self.cubic_len:
                    x, particle[3] = 2 * self.cubic_len - x, -particle[3]
                if y < 0:
                    y, particle[4] = -y, -particle[4]
                elif y > self.cubic_len:
                    y, particle[4] = 2 * self.cubic_len - y, -particle[4]
                # z loop
                if z < 0:
                    z, particle[5] = self.cubic_len + z, particle[5]
                elif z > self.cubic_len:
                    z, particle[5] = z - self.cubic_len, particle[5]
                # if z < 0:
                #     z, particle[5] = -z, -particle[5]
                # elif z > self.cubic_len:
                #     z, particle[5] = 2 * self.cubic_len - z, -particle[5]
                try:
                    self.cells[int(x // self.cell_len)][int(y // self.cell_len)][int(z // self.cell_len)].append(particle)
                except:
                    print(x, y, z)
                    print(int(x // self.cell_len), int(y // self.cell_len), int(z // self.cell_len))
                    raise
            
            # collision
            refresh_particles()            
            # random generate $\sigma\hat{\boldsymbol{sigma}}$
            θs = 2 * np.pi * np.random.random(self.n_particles)
            φs = np.arccos(2 * np.random.random(self.n_particles) - 1)
            σhats = np.array([np.sin(φs) * np.cos(θs), np.sin(φs) * np.sin(θs), np.cos(φs)]).T
            σs = 2 * self.r * np.array([np.sin(φs) * np.cos(θs), np.sin(φs) * np.sin(θs), np.cos(φs)]).T
            accepts = np.random.random(self.n_particles)
            self.cells = [[[[] for _ in range(self.cubic_size)] for _ in range(self.cubic_size)] for _ in range(self.cubic_size)]
            for i, particle in enumerate(particles):
                x, y, z = particle[:3] = bound(particle[0]), bound(particle[1]), cycle_bound(particle[2])                
                σ = σs[i]
                xJ, yJ, zJ = x + σ[0], y + σ[1], z + σ[2]
                xJ, yJ, zJ = bound(xJ), bound(yJ), cycle_bound(zJ)
                cell = self.cells[int(xJ // self.cell_len)][int(yJ // self.cell_len)][int(zJ // self.cell_len)]
                if len(cell) == 0:
                    pass
                else:  # random choose a particle j in the cell
                    n_J = len(cell)
                    j = np.random.randint(len(cell))
                    particle_j = cell[j]
                    g_ij = particle[3:] - particle_j[3:]
                    σhat_i_dot_g_ij = np.dot(σhats[i], g_ij)
                    if σhat_i_dot_g_ij > 0:  # Θ(σhat_i_dot_g_ij) = 1
                        ω_ij = 16 * np.pi * self.r**2 * σhat_i_dot_g_ij * (n_J / self.cell_volume) * dt
                        # assert 0 <= ω_ij <= 1, 'ω_ij = {}'.format(ω_ij)
                        self.ωs.append(ω_ij)
                        if accepts[i] < ω_ij:
                            # print('ω_ij = {}'.format(ω_ij))
                            particle[3:] -= σ                        
                self.cells[int(x // self.cell_len)][int(y // self.cell_len)][int(z // self.cell_len)].append(particle)
                
            # temperature
            refresh_particles()
            particles_array = np.array(particles)
            E_k = np.sum(particles_array[:, 3:] ** 2) * self.m / 2
            T = 2 / (3 * self.k_B * self.n_particles) * E_k
            self.Ts.append(T)
            if not quiet:
                pbar.set_description(f"T = {T:.6f} K")
            
            self.boxs = [[[[] for _ in range(self.cubic_box_size)] for _ in range(self.cubic_box_size)] for _ in range(self.cubic_box_size)]
            u = np.zeros((self.cubic_box_size, self.cubic_box_size, self.cubic_box_size, 3))
            pk = np.zeros((self.cubic_box_size, self.cubic_box_size, self.cubic_box_size, 3, 3))
            qk = np.zeros((self.cubic_box_size, self.cubic_box_size, self.cubic_box_size, 3))
            # self.pc = np.zeros((self.cubic_box_size, self.cubic_box_size, self.cubic_box_size, 3, 3))
            # self.qc = np.zeros((self.cubic_box_size, self.cubic_box_size, self.cubic_box_size, 3))
            for particle in particles:
                x, y, z = particle[:3]
                self.boxs[int(x // self.box_len)][int(y // self.box_len)][int(z // self.box_len)].append(particle)
            for i in range(self.cubic_box_size):
                for j in range(self.cubic_box_size):
                    for k in range(self.cubic_box_size):
                        box = self.boxs[i][j][k]
                        for particle in box:
                            u[i][j][k] += particle[3:]
                        u[i][j][k] = u[i][j][k] / len(box) if len(box) != 0 else 0
                        for particle in box:
                            dv = particle[3:] - u[i][j][k]
                            pk[i][j][k] += np.outer(dv, dv)
                            qk[i][j][k] += np.dot(dv, dv) * dv / 2
                        pk[i][j][k] *= self.m
                        qk[i][j][k] *= self.m
            self.us.append(u)
            self.pks.append(pk)
            self.qks.append(qk)
        # pbar.close()
                            
    def particles(self):
        particles = []
        for i in range(self.cubic_size):
                for j in range(self.cubic_size):
                    for k in range(self.cubic_size):
                        particles.extend(self.cells[i][j][k])
        # print(len(particles))
        return np.array(particles)
                            

