import numpy as np
import matplotlib.pyplot as plt 
import warnings as w

BC_OPTIONS = ['Periodic', 'Non-Periodic']
 
class NBody: 
    def __init__(self, size, ndim, particleList, dt, soft=0.1, G=1, boundaryType = 'Periodic'):

        """
        Class to carry out the simulation 
        Inputs: 
            - size (float): size of grid (should be a SQUARE)
            - ndim (int) : dimensions (2D or 3D)
            - particleList: initial list of particles along with their respective positions and velocities
            - dt (float): time step in the simulation
            - soft (float): softener used for the green function
            - G (float): gravitational constant (default 1.0, natural units)
            - boundaryType: either 'Periodic' or 'Non-Periodic' 
        """
        if boundaryType not in BC_OPTIONS: 
            raise ValueError("Boundary Condition must be one of {}".format(BC_OPTIONS))
        self.boundaryType = boundaryType
        #If Non-Periodic, we extend the grid by doubling it so particles dont feel the force from each other
        if size%2 == 0:
            self.ngrid = int(size)
        else: 
            w.warn("The size of the grid must be even, adding one to the input", RuntimeWarning)
            self.ngrid = int(size) + 1
        
        if self.boundaryType == 'Non-Periodic':
            self.ngrid = 2*self.ngrid

        self.ndim = int(ndim)
        self.particleList = particleList
        self.dt = dt
        self.positionP = particleList.position
        self.velocityP = particleList.velocity
        self.mass = particleList.mass
        
        x,y,z = np.arange(self.ngrid, dtype=float), np.arange(self.ngrid, dtype=float), np.arange(self.ngrid, dtype=float)
        
        if self.ndim == 2:
            self.mesh = np.array(np.meshgrid(x,y))
        elif self.ndim == 3:
            self.mesh = np.array(np.meshgrid(x,y,z))
            
        self.soft = soft
        self.G = G
        self.ngp()
        self.green()
        #xydensities = np.delete(self.densities, 1, axis=2)
        #print (xydensities.shape)
       # print (self.densities.shape)
        #fig,ax = plt.subplots(1,1)

        #ax.imshow(self.densities[:,:,0],origin='lower')
        #ax.plot(self.positionP[:,1], self.positionP[:,0],'o',color='red')
        #fig.savefig("testing_densities")
    
    def ngp(self):
        """
        Defining a function that assigns the density of the grid according to the 
        Nearest Grid Points (NGP) scheme. 
        We assign the mass of a given particle to its nearest gridpoint 
        (cells have unit length)
        """
        #Assign to closest grid point
        self.intPos = np.rint(self.positionP).astype('int') % self.ngrid
        
        #To access grid points
        self.mesh_modified = tuple(self.intPos[:, i] for i in range(self.ndim))

        #Bin the particles to their nearest grid point without looping
        edges = np.linspace(0, self.ngrid-1, num=self.ngrid+1)
        edges = np.repeat([edges], self.ndim, axis=0)
       
        #Generate histogram
        hist = np.histogramdd(self.intPos, bins=edges, weights=self.mass.flatten())
        self.densities = hist[0]

    def green(self):
        """
        Function that defines the greenfunction. 
        Note: this will only work for square grids as there was not enough time 
        to implement this with rectangular grids. 
        """
        
        r = np.sum(self.mesh**2, axis=0)
        r[r<self.soft**2] = self.soft**2
        r+= self.soft**2
        r = np.sqrt(r)

        g = 1/(4*np.pi*r)
        h_x,h_y,h_z = self.ngrid//2, self.ngrid//2, self.ngrid//2

        try:
            if self.ndim ==2: 
                g[h_x:, :h_y] = np.flip(g[:h_x,:h_y],axis=0)
                g[:,h_y:] = np.flip(g[:,:h_y],axis=1)

            elif self.ndim == 3:
                g[h_x:, :h_y, :h_z] = np.flip(g[:h_x, :h_y, :h_z],axis=0)
                g[:, h_y:, :h_z] = np.flip(g[:, :h_y, :h_z],axis=1)
                g[:, :, h_z:] = np.flip(g[:, :, :h_z],axis=2)
        except: 
            if self.ndim ==2:
                g[h_x:, :h_y+1] = np.flip(g[:h_x+1,:h_y+1],axis=0)
                g[:,h_y:] = np.flip(g[:,:h_y+1],axis=1)

            elif self.ndim ==3:
                g[h_x:, :h_y+1, :h_z+1] = np.flip(g[:h_x+1, :h_y+1, :h_z+1],axis=0)
                g[:, h_y:, :h_z+1] = np.flip(g[:, :h_y+1, :h_z+1],axis=1)
                g[:, :, h_z:] = np.flip(g[:, :, :h_z+1],axis=2)

        self.g = g
        
    def potential(self):
        """
        Getting the potential of the grid by convoluting the density 
        of the grid and the green function of the grid. 
        
        Recentering the particle once it is done to account for the fact that the particles 
        were all attracted to themselves and started to drift towards the origin.
        """
        ffD = np.fft.rfftn(self.densities)
        ffG = np.fft.rfftn(self.g)
        ffV = ffD*ffG
        V = np.fft.irfftn(ffV)

        #Centre the particle back to where it started 
        for i in range(self.ndim):
            #Perform shift
            V = 0.5*(np.roll(V,1,axis=i)+V)
        
        #If the boundary conditions is not periodic, set the potential to 0
        if self.boundaryType == 'Non-Periodic':
            pad = tuple((slice(1, -1),)) * self.ndim
            V = np.pad(V[pad],1)
            
        return V 
        
            
    def forces_mesh(self): 
        """
        Define function to get the forces take the gradient of the potential. 
        To take the gradient, use the central difference.
        """
        shape = [self.ndim]
        shape.extend([self.ngrid]*self.ndim)
        fmesh = np.zeros(shape,dtype=float)

        #get the potential
        V = self.potential()

        for i in range(self.ndim):
            fmesh[i] = 0.5*(np.roll(V,1,axis=i)-np.roll(V,-1,axis=i))
        
        #Multiply by gravitational constant and the densities to get the force
        fmesh = -fmesh*self.densities*self.G
        return fmesh 
        
    
    def forces_pctls(self):
        """
        Interpolating the forces using the inverse scheme of the NGP density scheme. 
        """
        fxyz = np.moveaxis(self.forces_mesh(),0,-1)
        f_final = fxyz[self.mesh_modified]
         
        return f_final
    
    def totalEnergy(self):
        """
        Function to compute the total energy of the system
        E = V + 0.5*m*v**2
        """
        K = 0.5*np.sum(self.mass*self.velocityP**2)
        P = -0.5*np.sum(np.sum(self.potential())*self.densities)
        return K+P
    

    def evolve_system(self,nsteps=1,file_save=None,file_save_pos=None):
        """
        Evolves the system and saves the energy along with position for further 
        analysis. The evolution uses the leapfrog method.
        
        Inputs:
            - nsteps (int): how many steps per evolution before it saves a result
            - file_save (file): File in which we want to save the total energy of the system
            - file_save_pos (file,int): First entry is the file in which we want to save the 
                                        energy is, second entry is the number of particles we 
                                        want to track.
        """
        for i in range(nsteps):
            #particle mesh method to calculate force
            F = self.forces_pctls()

            #leapfrog method to update position and velocity 
            self.velocityP += F*self.dt/self.mass
            self.positionP += self.velocityP*self.dt
            self.positionP = self.positionP% self.ngrid
            
            #Updates the grid again with the new positions
            self.ngp()

        energy = self.totalEnergy()
        
        if file_save != None:
            file_save.write(f"{energy}\n")
            file_save.flush()
        if file_save_pos != None:
            file_save_pos[0].write(f"{self.positionP[:file_save_pos[1]]}\n")
            file_save_pos[0].flush()
            
        return energy,self.positionP
