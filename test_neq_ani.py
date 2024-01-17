from mc import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# from mpl_toolkits.mplot3d import Axes3D

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

fig = plt.figure()
ax = fig.add_subplot(231, projection='3d')
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

u_x, u_y, u_z, pk_xx, pk_yy, pk_zz, pk_xy, pk_yz, pk_zx, pk_yx, pk_zy, pk_xz, qk_x, qk_y, qk_z = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

def update(num):
    ax.clear()
    cubic.run(n_steps=1, dt=dt, quiet=True)
    ax.scatter(cubic.particles()[:, 0], cubic.particles()[:, 1], cubic.particles()[:, 2], s=1)
    ax.set_xlim([0, cubic.cubic_len])
    ax.set_ylim([0, cubic.cubic_len])
    ax.set_zlim([0, cubic.cubic_len])
    ax2.clear()
    pk_xx.append(np.mean(cubic.pks[-1][:,0,0]))
    pk_yy.append(np.mean(cubic.pks[-1][:,1,1]))
    pk_zz.append(np.mean(cubic.pks[-1][:,2,2]))
    ax2.plot(cubic.ts, pk_xx, label=r'$p^k_{xx}$')
    ax2.plot(cubic.ts, pk_yy, label=r'$p^k_{yy}$')
    ax2.plot(cubic.ts, pk_zz, label=r'$p^k_{zz}$')
    ax2.legend()
    ax3.clear()
    u_x.append(np.mean(cubic.us[-1][:,0]))
    u_y.append(np.mean(cubic.us[-1][:,1]))
    u_z.append(np.mean(cubic.us[-1][:,2]))
    ax3.plot(cubic.ts, u_x, label=r'$u_x$')
    ax3.plot(cubic.ts, u_y, label=r'$u_y$')
    ax3.plot(cubic.ts, u_z, label=r'$u_z$')
    ax3.legend()
    ax4.clear()
    qk_x.append(np.mean(cubic.qks[-1][:,0]))
    qk_y.append(np.mean(cubic.qks[-1][:,1]))
    qk_z.append(np.mean(cubic.qks[-1][:,2])) 
    ax4.plot(cubic.ts, qk_x, label=r'$q_x$')
    ax4.plot(cubic.ts, qk_y, label=r'$q_y$')
    ax4.plot(cubic.ts, qk_z, label=r'$q_z$')  
    ax4.legend()
    ax5.clear()
    pk_xy.append(np.mean(cubic.pks[-1][:,0,1]))
    pk_yz.append(np.mean(cubic.pks[-1][:,1,2]))
    pk_zx.append(np.mean(cubic.pks[-1][:,2,0]))
    pk_yx.append(np.mean(cubic.pks[-1][:,1,0]))
    pk_zy.append(np.mean(cubic.pks[-1][:,2,1]))
    pk_xz.append(np.mean(cubic.pks[-1][:,0,2]))
    ax5.plot(cubic.ts, pk_xy, label=r'$p^k_{xy}$')
    ax5.plot(cubic.ts, pk_yz, label=r'$p^k_{yz}$')
    ax5.plot(cubic.ts, pk_zx, label=r'$p^k_{zx}$') 
    ax5.plot(cubic.ts, pk_yx, label=r'$p^k_{yx}$')
    ax5.plot(cubic.ts, pk_zy, label=r'$p^k_{zy}$')
    ax5.plot(cubic.ts, pk_xz, label=r'$p^k_{xz}$')
    ax5.legend()
    ax6.clear()
    ax6.plot(cubic.ts, cubic.Ts)
    ax6.set_ylabel('Temperature')

ani = animation.FuncAnimation(fig, update, frames=100, repeat=False)
# set figure title
fig.suptitle(f'{n_particles} particles, cell_len={cell_len}, cubic_size={cubic_size}, box_size={box_size}, r={r}, m={m}, T0={T}, n_steps={n_steps}, dt={dt}')
plt.show()
# ani.save(f'{n_particles}_{cell_len}_{cubic_size}_{box_size}_{r}_{m}_{T}_{n_steps}_{dt}.mp4')
