import numpy as np 
from scipy.ndimage import gaussian_filter

#Defining functions for the noise models and whitening

#Defining a model for the noise
def noise_model(strain, window):
    #Multiply the strain by the window and get a normalization factor
    window_strain = strain*window
    n = np.sqrt(np.mean(window**2))
    #Find power spectrum from the |FT|^2 
    pow_spec = np.abs(np.fft.rfft(window_strain/n))**2
    #Perform gaussian smoothing
    pow_spec = gaussian_filter(pow_spec,1)

    return pow_spec

#Defining a function for whitening
def whitening(data, noise, window):
    #Get a normalization factor
    n = np.sqrt(np.mean(window**2))
    
    #Return the whitened strain or template using the window
    return np.fft.rfft(window*data)/(np.sqrt(noise)*n)