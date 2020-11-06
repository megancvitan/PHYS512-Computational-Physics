import numpy as np 
import noise

#Function to plot noise model
def plot_noise(event, strain_H, strain_L, t_H, _L, window, i, fig, ax):
    # i is the loop argument, fig and ax are the the plots argument

    #Power spectrum for L and H
    pow_H = noise.noise_model(strain_H, window)
    pow_L = noise.noise_model(strain_L, window)
    
    #Plotting for H and L
    ax[i,0].semilogy(pow_H)
    ax[i,0].set_ylabel(event+'_Han')
    ax[i,1].semilogy(pow_L)
    ax[i,1].set_ylabel(event+'_Liv')
    
    if i==3:
        ax[i,0].set_xlabel('frequency (Hz)')
        ax[i,1].set_xlabel('frequency (Hz)')


#Function to get the matched filters
def find_mf(event, strain_H, strain_L, t_H, t_L, window):

    #Power spectrum for L and H from noise.py
    pow_H = noise.noise_model(strain_H, window)
    pow_L = noise.noise_model(strain_L, window)

    #Whiten the data with the whitening() in noise.py
    A_white_H = noise.whitening(t_H, pow_H, window)
    A_white_L = noise.whitening(t_L, pow_L, window)
    d_white_H = noise.whitening(strain_H, pow_H, window)
    d_white_L = noise.whitening(strain_L, pow_L, window)

    #Run and return the matched filters on H and L
    mf_H = np.fft.fftshift(np.fft.irfft(np.conj(A_white_H)*d_white_H))
    mf_L = np.fft.fftshift(np.fft.irfft(np.conj(A_white_L)*d_white_L))

    return mf_H, mf_L

#Function to plot the matched filters
def plot_mf(event, strain_H, strain_L, t_H, t_L, window, Time, i, fig, ax):
    #Collect the matched filters with previously defined function
    mf_H, mf_L = find_mf(event, strain_H, strain_L, t_H, t_L, window)
    
    #Plotting the matched filters
    ax[i,0].plot(Time, mf_H)
    ax[i,0].set_ylabel(event+'_Han')
    ax[i,1].plot(Time, mf_L)
    ax[i,1].set_ylabel(event+'_Liv')
    
    if i==3:
        ax[i,0].set_xlabel('time (s)')
        ax[i,1].set_xlabel('time (s)')
        
