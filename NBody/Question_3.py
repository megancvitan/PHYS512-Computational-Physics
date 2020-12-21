import numpy as np
import matplotlib.pyplot as plt
import particle_properties as P
import nbody as nb
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  

def periodicInitial():
    def animate(i):
        global g,ax,fig,file_energy
        model.evolve_system(nsteps=10,file_save=energy_file)
        ptcl.set_data(model.densities)
        ptcl.set_cmap(plt.get_cmap('magma'))
        print ("Step size {}".format(i))
        return ptcl,
 
    def animate3D(i): 
        global g,ax,fig,file_energy 
        model.evolve_system(nsteps=10,file_save = energy_file)
        ptcl.set_data(model.positionP[:,0], model.positionP[:,1])
        ptcl.set_3d_properties(model.positionP[:,2])
        print ("Step size{}".format(i))
        return ptcl, 

    npart = 2**13
    gridsize = 2**7
    ndim = 3
    velInit = 0 
    mass = 1/npart
    init_mas = [mass for t in range(npart)]
    soft = 0.8
    s = P.initial_system(npart,gridsize,ndim,init_mas,n_part_specific=None,n_part_specficVel=velInit)
    dt = 5
    model = nb.NBody(gridsize, ndim, s,dt,soft=soft)

    niter = 50

    if ndim == 2: 
        Title = 'Question3_Periodic_2D.gif'
        energy_file = open('Question3_Periodic_Energy_2D.txt','w')
    elif ndim == 3: 
        Title = 'Question3_Periodic_3D.gif'
        energy_file = open('Question3_Periodic_Energy_3D.txt','w')

    T = r'Simulation with $2^{13}$ Particles with $m$=1.2e-4, $dt$=5 with 50 frames'
    x = []
    y = []

    labelsize = 15
    fig = plt.figure(figsize=(10,6))
    if ndim ==2:
        ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,gridsize),ylim=(0,gridsize))
    elif ndim == 3: 
        ax = fig.add_subplot(111,projection='3d')
    ax.tick_params(labelsize=labelsize)
    ax.set_xlabel("X Position",fontsize=labelsize)
    ax.set_ylabel("Y Position",fontsize = labelsize)
    if ndim == 3:
        ax.set_zlabel("Z Position", fontsize= labelsize)
        ax.xaxis.labelpad = 40
        ax.yaxis.labelpad = 40
        ax.zaxis.labelpad = 40
        ax.set_xlim(0, gridsize)
        ax.set_ylim(0, gridsize)
        ax.set_zlim(0, gridsize)
    ax.set_title(T,fontsize=labelsize)

    if ndim ==2 :
        ptcl = ax.imshow(model.densities,origin='lower',vmin=model.densities.min(),vmax=model.densities.max())
        plt.colorbar(ptcl)
        an = animation.FuncAnimation(fig,animate,frames=niter,interval=10,repeat=True)
    elif ndim ==3: 
        ptcl, = ax.plot(model.positionP[:,0], model.positionP[:,1], model.positionP[:,2], ".", markersize = 1)
        an = animation.FuncAnimation(fig,animate3D,frames=niter,interval=10,repeat=True)

    an.save(Title, writer='imagemagick')


def nonPeriodicInitial(): 
    def animate(i):
        global g,ax,fig,file_energy
        model.evolve_system(nsteps=10,file_save=energy_file)
        ptcl.set_data(model.densities)
        ptcl.set_cmap(plt.get_cmap('magma'))
        print ("Step size {}".format(i))
        return ptcl,
 
    def animate3D(i): 
        global g,ax,fig,file_energy 
        model.evolve_system(nsteps=10,file_save = energy_file)
        ptcl.set_data(model.positionP[:,0], model.positionP[:,1])
        ptcl.set_3d_properties(model.positionP[:,2])
        print ("Step size{}".format(i))
        return ptcl, 

    npart = 2**12
    gridsize = 2**7
    ndim = 3
    velInit = 0 
    mass = 1/npart
    init_mas = [mass for t in range(npart)]
    soft = 0.8
    s = P.initial_system(npart,gridsize,ndim,init_mas,n_part_specific=None,n_part_specficVel=velInit,boundaryType='Non-Periodic')
    dt = 10
    model = nb.NBody(gridsize, ndim, s,dt,soft=soft,boundaryType='Non-Periodic')

    niter = 75

    if ndim == 2: 
        Title = 'Question3_NonPeriodic_2D.gif'
        energy_file = open('Question3_Non_Periodic_Energy_2D.txt','w')
    elif ndim == 3: 
        Title = 'Question3_NonPeriodic_3D.gif'
        energy_file = open('Question3_Non_Periodic_Energy_3D.txt','w')

    T = r'Simulation with $2^{12}$ Particles with $m$=1.2e-4, $dt$=10 with 75 frames'
    x = []
    y = []

    labelsize = 15
    fig = plt.figure(figsize=(10,6))
    if ndim ==2:
        ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,gridsize),ylim=(0,gridsize))
    elif ndim == 3: 
        ax = fig.add_subplot(111,projection='3d')
    ax.tick_params(labelsize=labelsize)
    ax.set_xlabel("X Position",fontsize=labelsize)
    ax.set_ylabel("Y Position",fontsize = labelsize)
    if ndim == 3:
        ax.set_zlabel("Z Position", fontsize= labelsize)
        ax.xaxis.labelpad = 40
        ax.yaxis.labelpad = 40
        ax.zaxis.labelpad = 40
        ax.set_xlim(0, gridsize)
        ax.set_ylim(0, gridsize)
        ax.set_zlim(0, gridsize)
    ax.set_title(T,fontsize=labelsize)

    if ndim ==2 :
        ptcl = ax.imshow(model.densities,origin='lower',vmin=model.densities.min(),vmax=model.densities.max())
        plt.colorbar(ptcl)
        an = animation.FuncAnimation(fig,animate,frames=niter,interval=10,repeat=True)
    elif ndim ==3: 
        ptcl, = ax.plot(model.positionP[:,0], model.positionP[:,1], model.positionP[:,2], ".", markersize = 1)
        an = animation.FuncAnimation(fig,animate3D,frames=niter,interval=10,repeat=True)


    an.save(Title, writer='imagemagick')
  
nonPeriodicInitial()
#periodicInitial()
