import random as rn 
import numpy as np
import matplotlib.pyplot as plt 

class particle:
    def __init__(self, mass, posX, posY, posZ=None, velX = 0, velY = 0, velZ = 0):

        """
        Class to define a particle. Every particle has a position, a velocity, and a mass
        Input parameters:
            - mass (float): mass of the particle
            - posX, posY, posZ (float, float, float): x,y,z components of a particle's position
            - velX, velY, velZ (float, float, float): x,y,z components of a particle's velocity
        """
        
        #Initializing the mass, position and velocity for a given particle in 3D
        self.mass = mass
        if posZ is not None:
            self.position = (posX, posY, posZ)
            self.velocity = (velX, velY, velZ)
        else: 
            self.position = (posX,posY)
            self.velocity = (velX, velY)

class initial_system: 
    def __init__(self, nparticles, size, ndim,  init_mass, n_part_specific = None, n_part_specficVel = None, boundaryType = 'Periodic', soft=None, cosmos=False):
        """
        Class to define a system. 
        Input parameters:
            - nparticles (int): number of particles in the system
            - size (float): size of the grid in (size * size * size)
            - ndim (float): dimension of the problem
            - init_mass (float): given mass of the particle
            - n_part_specific ([][]) : to initialize specific particles positions, default is None
            - n_part_specificVel ([][]) : to initialize specific particles velocities, default is None
            - boundaryType (string): either 'Periodic' or 'Non-Periodic' boundary conditions
            - soft (float): for smoothing to avoid infinite behaviour, default None
            - cosmos (boolean): True if simulating the entire universe, default False if other simulation
        """
        #Specify the cosmos and boundary type
        self.cosmos = cosmos
        self.boundaryType = boundaryType
        self.ngrid = size
        self.ndim = ndim

        if self.ndim != 3 and self.ndim !=2: 
            raise ValueError("This simulation is only implemented in either 2D or 3D")

        def initial_pos():
            #Subroutine to generate particle's position in a grid
            #Generate random particles and also particles that are specified as inputs
            if self.boundaryType == 'Periodic':
                init_cond = []
                #Defines specific particles. If None are given, generate random positions
                if n_part_specific is None: 
                    l = 0
                else:
                    l = len(n_part_specific)
                    for k in range(l):
                        #For a 3D system
                        if self.ndim == 3:
                            pos = (n_part_specific[k][0], n_part_specific[k][1], n_part_specific[k][2])
                        elif self.ndim == 2: 
                            pos = (n_part_specific[k][0], n_part_specific[k][1])
                        #Save the initial positions in all three directions
                        init_cond.append(pos)
                for p in range(nparticles-l):
                    #Generate a random position for all three axes
                    if self.ndim ==3: 
                        pos = (rn.random())*(self.ngrid-1), (rn.random())*(self.ngrid-1), (rn.random())*(self.ngrid-1)
                    elif self.ndim ==2:
                        pos = ((rn.random())*(self.ngrid-1), (rn.random())*(self.ngrid-1))

                    init_cond.append(pos)
            #For non-periodic, we don't want particles to be generated on the boundaries of the grid
            elif self.boundaryType == 'Non-Periodic':
                xmin, xmax = 1, self.ngrid-1
                ymin, ymax = 1, self.ngrid-1
                if self.ndim ==3 :
                    zmin, zmax = 1, self.ngrid-1
                init_cond = []
                if n_part_specific is None: 
                    l = 0
                else: 
                    l = len(n_part_specific)
                    for k in range(l):
                        #Print out error message if the position is outside the grid
                        if self.ndim ==3:
                            if n_part_specific[k][0].max() > xmax or n_part_specific[k][1].max() > ymax or n_part_specific[k][2].max() > zmax or n_part_specific[k][0].min() < xmin or n_part_specific[k][1].min() < ymin or n_part_specific[k][2].min() < zmin: 
                                raise ValueError('The position of the given particle must be within the boundaries of the grid')
                            pos = (n_part_specific[k][0], n_part_specific[k][1], n_part_specific[k][2])
                        if self.ndim ==2: 
                            if n_part_specific[k][0].max() > xmax or n_part_specific[k][1].max() > ymax or n_part_specific[k][0].min() < xmin or n_part_specific[k][1].min() < ymin: 
                                raise ValueError('The position of the given particle must be within the boundaries of the grid')
                            pos = (n_part_specific[k][0], n_part_specific[k][1])
                      
                        init_cond.append(pos)
                for p in range(nparticles-l): 
                    #Generate a random number between the given bounds for each direction
                    if self.ndim ==3:
                        pos = (rn.uniform(1.0001, self.ngrid- 1.0001), rn.uniform(1.0001, self.ngrid-1), rn.uniform(1.0001, self.ngrid-1))
                    if self.ndim ==2:
                        pos = (rn.uniform(1.0001, self.ngrid- 1.0001), rn.uniform(1.0001, self.ngrid-1))

                    init_cond.append(pos)
            return init_cond

        def initial_vel(maxSpeed): 
            #Subroutine to generate particle's velocity in a grid
            init_cond = []
            if np.isscalar(n_part_specficVel) == False: 
                if n_part_specific is None: 
                    l = 0 
                else: 
                    l = len(n_part_specficVel)
                    for k in range(l):
                        if self.ndim == 3:
                            vel = (n_part_specficVel[k][0], n_part_specficVel[k][1], n_part_specficVel[k][2])
                        elif self.ndim ==2:
                            vel = (n_part_specficVel[k][0], n_part_specficVel[k][1])

                        init_cond.append(vel)
                for p in range(nparticles-l):
                    #Generate a random value for each direction
                    if self.ndim == 3: 
                        vel = (np.random.normal()*maxSpeed, np.random.normal()*maxSpeed, np.random.normal()*maxSpeed)
                    elif self.ndim == 2: 
                        vel = (np.random.normal()*maxSpeed, np.random.normal()*maxSpeed)

                    init_cond.append(vel)
            else: 
                for p in range(nparticles):
                    #Initialize to zero and save the velocity
                    if self.ndim == 2: 
                        vel = (0.0,0.0)
                    elif self.ndim == 3: 
                        vel = (0.0,0.0,0.0)
                    init_cond.append(vel)
            return init_cond

        #Subroutine to initialize the mass of a particle
        def initial_mass(posP, soft):
            #If we are not simulating the universe
            if self.cosmos == False:
                #Transpose the matrix and return
                return np.array([init_mass.copy()]).T
            else: 
                #Need to calculate the mass; first we need to find the positions
                positionPx = np.rint(posP[:,0]).astype('int') % self.ngrid
                positionPy = np.rint(posP[:,1]).astype('int') % self.ngrid
                if self.ndim ==3:
                    positionPz = np.rint(posP[:,2]).astype('int') % self.ngrid

                #Then find the power spectrum
                #Need to do fourier transforms of positions
                k_x = np.real(np.fft.fft(positionPx))
                k_y = np.real(np.fft.fft(positionPy))
                if self.ndim ==3:
                    k_z = np.real(np.fft.fft(positionPz))
                
                #Get the total k 
                if self.ndim == 3:
                    k_total = np.sqrt(k_x**2+k_y**2+k_z**2)
                elif self.ndim == 2:
                    k_total = np.sqrt(k_x**2+k_y**2)
                    
                #Consider the softening parameter to prevent any funky behaviour
                k_total[k_total<soft] = soft
                #Finally calculate the mass value
                m = init_mass/k_total**3
                
                return np.array([m.copy()]).T

        #Initialize the rest of the properties given the constructors above
        init_cond = initial_pos()
        init_vel = initial_vel(1)
        self.nparticles = nparticles
        
        if self.ndim == 3:
            self.particles = np.asarray([particle(m,s[0],s[1],s[2],velX=v[0],velY=v[1],velZ=v[2]) for m,s,v in zip(init_mass,init_cond,init_vel)])   
        elif self.ndim ==2:
            self.particles = np.asarray([particle(m,s[0],s[1],velX=v[0],velY=v[1]) for m,s,v in zip(init_mass,init_cond,init_vel)])   

        self.velocity = np.asarray([self.particles[i].velocity for i in range(nparticles)])
        self.position = np.asarray([self.particles[i].position for i in range(nparticles)])
        self.mass = initial_mass(posP=self.position,soft=soft)
