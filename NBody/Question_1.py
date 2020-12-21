import numpy as np
import matplotlib.pyplot as plt
import particle_properties as p
import nbody as nb
import matplotlib.animation as animation

nparticles = 1
gridsize = 2**5
ndim = 3
mass = 10
init_mass = [mass for t in range(nparticles)]
pos = [(gridsize/2 ,gridsize/2, gridsize/2)]
vel = [(0, 0, 0)]

system = p.initial_system(nparticles, gridsize, ndim,  init_mass, n_part_specific = pos, n_part_specficVel = vel)

dt = 1
soft = 0.1
niter = 100

model = nb.NBody(gridsize, ndim, system, dt, soft=soft)

def animate(i):
    global system,ax,fig
    model.evolve_system()
    ptcl.set_data(model.positionP[:,0],model.positionP[:,1])
    return ptcl,

title = 'Question_1.gif'
T = f'Stationary object with $v_o$=0, $m$={mass}, $dt$={dt} with {niter} frames'

labelsize = 15
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,model.ngrid),ylim=(0,model.ngrid))
ax.tick_params(labelsize=labelsize)
ax.set_xlabel("X Position",fontsize=labelsize)
ax.set_ylabel("Y Position",fontsize = labelsize)
ax.set_title(T,fontsize=labelsize)

ptcl, = ax.plot([],[],'*',markersize=10,color='black')

an = animation.FuncAnimation(fig,animate,frames=niter,interval=dt,repeat=True)
an.save(title, writer='imagemagick')