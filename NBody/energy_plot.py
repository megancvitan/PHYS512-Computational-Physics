import numpy as np
import matplotlib.pyplot as plt

#Give file name with energy data
filename = "Question3_Non_Periodic_Energy_3D.txt"
niter = 50

def generate_energy_plots(filename, niter):
    file = open("./energy_data/{}".format(filename),'r')
    energies = [float(file.readline()) for i in range(niter)]
    fig,ax = plt.subplots(1,1,figsize=(20,12))
    ax.set_xlabel("Frame", fontsize=15)
    ax.set_ylabel("Energy", fontsize=15)
    ax.set_title(filename, fontsize=15)
    ax.tick_params(labelsize=15)
    ax.plot(energies, '.', color='black')
    
generate_energy_plots(filename, niter)