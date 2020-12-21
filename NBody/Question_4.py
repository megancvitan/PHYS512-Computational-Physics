import numpy as np
import matplotlib.pyplot as plt
import particle_properties as P
import nbody as nb
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm, Normalize
    
def periodic_cosmos():
    def animate(i):
        global g,ax,fig,file_energy
        model.evolve_system(nsteps=10,file_save=energy_file)
        ptcl.set_data(model.densities)
        ptcl.set_cmap(plt.get_cmap('magma'))
        print ("Step size {}".format(i))
        return ptcl,
    
    nparticles = 2**16
    gridsize = 2**8
    size = (gridsize,gridsize,gridsize)
    mass = 40
    init_mass = [mass for t in range(nparticles)]
    init_velocity = 0
    dt = 330
    soft = 10
    niter = 450
    ndim = 2
    
    s = P.initial_system(nparticles,gridsize,ndim,init_mass,n_part_specific=None,soft=soft,n_part_specficVel=init_velocity, cosmos=True)
    energy_file = open('Question4_Energy_2D.txt','w')
    model = nb.NBody(gridsize, ndim, s, dt, soft=soft)

    Title = 'Question_4.gif'
    T = f'Orbiting objects with $v_o$=0.1, $m$={mass}, $dt$={dt} with {niter} frames'
    x = []
    y = []

    labelsize = 15
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,gridsize), ylim=(0,gridsize))
    ax.tick_params(labelsize=labelsize)
    ax.set_xlabel("X Position", fontsize=labelsize)
    ax.set_ylabel("Y Position", fontsize = labelsize)
    ax.set_title(T, fontsize=labelsize)
    
    modelCopy = model.densities.copy()
    modelCopy[modelCopy==0] = modelCopy[modelCopy!=0].min()*1e-3
    
    #ptcl, = ax.plot(bod.positionP[:,0], bod.positionP[:,1], bod.positionP[:,2], 'o', markersize = 2 ,color= 'blue')
    #ptcl = ax.imshow(model.densities,origin='lower',norm=LogNorm(vmin=modelCopy.min(),vmax=modelCopy.max()))
    ptcl = ax.imshow(model.densities,origin='lower',vmin=model.densities.min(),vmax=model.densities.max())    
    plt.colorbar(ptcl)

    an = animation.FuncAnimation(fig,animate,frames=niter,interval=10,repeat=True)
    an.save(Title, writer='imagemagick')

periodic_cosmos()