import numpy as np 
import noise
import waves

#Function to find the SNRs including H, L, H+L
def SNR(event, strain_H, strain_L, t_H, t_L, window):
    #Get the noise for L and H
    pow_H = noise.noise_model(strain_H, window)
    pow_L = noise.noise_model(strain_L, window)
    
    #Get the whitened strains
    A_white_H = noise.whitening(t_H, pow_H, window)
    A_white_L = noise.whitening(t_L, pow_L, window)

    #Get the matched filters
    mf_H, mf_L = waves.find_mf(event, strain_H, strain_L, t_H, t_L, window)

    #Compute and return the SNRs
    SNR_H = np.abs(mf_H*np.fft.fftshift(np.fft.irfft(np.sqrt(np.conj(A_white_H)*A_white_H))))
    SNR_L = np.abs(mf_L*np.fft.fftshift(np.fft.irfft(np.sqrt(np.conj(A_white_L)*A_white_L))))
    SNR_HL = np.sqrt(SNR_H**2 + SNR_L**2)

    return SNR_H, SNR_L, SNR_HL

#Function to plot the SNRs
def plot_SNR(event, strain_H, strain_L, t_H, t_L, window, Time, i, fig, ax):
    #Collect the SNRs from the SNR() function above
    S_H, S_L, S_HL = SNR(event, strain_H, strain_L, t_H, t_L, window)

    #Ploting the results
    ax[i,0].plot(Time, S_H)
    ax[i,0].set_ylabel(event+'_Han')
    ax[i,1].plot(Time, S_L)
    ax[i,1].set_ylabel(event+'_Liv')
    ax[i,2].plot(Time, S_HL)
    ax[i,2].set_ylabel(event+'_H+L')
    
    if i==3:
        ax[i,0].set_xlabel('time (s)')
        ax[i,1].set_xlabel('time (s)')
        ax[i,2].set_xlabel('time (s)')

# This function computes the theoretical SNRs through an analytical approach
def theory_SNR(event, strain_H, strain_L, t_H, t_L, window):
    #Get the noise for L and H
    pow_H = noise.noise_model(strain_H, window)
    pow_L = noise.noise_model(strain_L, window)
    
    #Get the whitened strains
    A_white_H = noise.whitening(t_H, pow_H, window)
    A_white_L = noise.whitening(t_L, pow_L, window)

    #Computing the analytical SNR with the whitened model
    #This is done by taking A_white_H in real space
    theory_SNR_H = np.abs(np.fft.irfft(A_white_H))
    theory_SNR_L = np.abs(np.fft.irfft(A_white_L))
    theory_SNR_HL = np.sqrt(theory_SNR_H**2 + theory_SNR_L**2)

    return theory_SNR_H, theory_SNR_L, theory_SNR_HL

#Function to plot the theoretically found SNR
def plot_theory_SNR(event, strain_H, strain_L, t_H, t_L, window, Time, i, fig, ax):
    #Using theory_SNR() to get the data 
    S_H, S_L, S_HL = theory_SNR(event, strain_H, strain_L, t_H, t_L, window)

    #Plotting the data
    ax[i,0].plot(Time, S_H)
    ax[i,0].set_ylabel(event+'_Han')
    ax[i,1].plot(Time, S_L)
    ax[i,1].set_ylabel(event+'_Liv')
    ax[i,2].plot(Time, S_HL)
    ax[i,2].set_ylabel(event+'_H+L')
    
    if i==3:
        ax[i,0].set_xlabel('time (s)')
        ax[i,1].set_xlabel('time (s)')
        ax[i,2].set_xlabel('time (s)')

def gauss(x, u, A, sig):
    return A*np.exp(-(x-u)**2/(2.0*sig**2))
