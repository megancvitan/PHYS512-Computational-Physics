import numpy as np
import matplotlib.pyplot as plt
import particle_properties as p
import nbody as nb
import matplotlib.animation as animation

nparticles = 2
gridsize = 2**8
ndim = 3
mass = 5
init_mass = [mass for t in range(nparticles)]
pos = np.array([[gridsize/2 ,gridsize/2+10, 0], [gridsize/2, gridsize/2-10, 0]])
vel = np.array([[0.1, 0, 0], [-0.1,0,0]])

system = p.initial_system(nparticles, gridsize, ndim, init_mass, n_part_specific = pos, n_part_specficVel = vel)

dt = 1
soft = 0.1
niter = 200

x = []
y = []

model = nb.NBody(gridsize, ndim, system, dt, soft=soft)
print ("Finished calculating the densities and the green function")

title = 'Question_2.gif'
T = f'Orbiting objects with $v_o$=0.1, $m$={mass}, $dt$={dt} with {niter} frames'

def animate(i):
    global g,ax,fig
    model.evolve_system()
    x.append(model.positionP[:,0][0])s
    y.append(model.positionP[:,1][0])
    ptcl.set_data(model.positionP[:,0],model.positionP[:,1])
    return ptcl,

labelsize = 15
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,model.ngrid),ylim=(0,model.ngrid))
ax.tick_params(labelsize=labelsize)
ax.set_xlabel("X Position",fontsize=labelsize)
ax.set_ylabel("Y Position",fontsize = labelsize)
ax.set_title(T,fontsize=labelsize)

ptcl, = ax.plot([],[],'o',markersize=10,color='blue')

an = animation.FuncAnimation(fig,animate,frames=niter,interval=dt,repeat=False)
#Writer = animation.writers['ffmpeg']
#Writer = animation.writers['imagemagick']

#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
an.save(title, writer='imagemagick',fps=30)
#an.save(title, writer=writer)

print (y)
