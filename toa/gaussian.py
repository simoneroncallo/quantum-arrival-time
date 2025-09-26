import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

class GaussianTrain():
    """ Build, manage and visualize a superposition of Gaussian wave packets. """
    def __init__(self, X, P, S, M, nospread = False):
        """ Initialize Gaussian wave packets.
        Parameters:
        X : array-like
            Initial positions of the packets.
        P : array-like
            Initial momenta of the packets.
        S : array-like
            Standard deviations of the packets.
        M : array-like
            Masses of the packets.
        nospread : bool, optional
            Disables packet spreading (default: False).
        """
        self.X = np.array(X, dtype = float)
        self.P = np.array(P, dtype = float)
        self.S = np.array(S, dtype = float)
        self.M = np.array(M, dtype = float)
        self.N = np.size(self.X)
        self.nospread = nospread

    def packet(self, tVar, xVar, idx):
        """ Evaluate single Gaussian wave packet.
        Parameters:
        tVar : float or array-like
            Time variable.
        xVar : float or array-like
            Position variable.
        idx : int
            Packet index.
        Returns:
        wavepckt : complex or np.ndarray
            Value of the Gaussian packet at (tVar, xVar).
        """

        t = np.array(tVar)
        x = np.array(xVar)
        if self.nospread:
            t = np.zeros_like(tVar)
            x = xVar - self.P[idx]*tVar/self.M[idx]
            
        wavepckt = 1/( np.pi**(1/4)*np.sqrt( self.S[idx] + 1j*t/(self.S[idx]*self.M[idx]) ) ) # Amplitude
        wavepckt *= np.exp( -1*( (x - self.X[idx] - self.P[idx]*t/self.M[idx])**2 )/\
                          ( (2*self.S[idx]**2)*(1 + 1j*t/( self.M[idx]*self.S[idx]**2 ) ) ) )
        wavepckt *= np.exp( 1j*self.P[idx]*( x - self.X[idx] - self.P[idx]*t/(2*self.M[idx]) ) )
        return wavepckt

    def superposition(self, tVar, xVar):
        """ Evaluate superposition of Gaussian wave packets.
        Parameters:
        tVar : float or array-like
            Time variable.
        xVar : float or array-like
            Position variable.
        Returns:
        trainpckt : complex or np.ndarray
            Value of the superposition at (tVar, xVar).
        """
        trainpckt = 0
        for idx in range(self.N):
            trainpckt += self.packet(tVar, xVar, idx)
        return trainpckt/np.sqrt(self.N)

    def visualize(self, numPoints, tLim, xLim, xDtc = 0):
        """ Visualize wave function with density and fixed-position plotw.
        Parameters:
        numPoints : int
            Number of points in each dimension (x and t).
        tLim : tuple or array-like
            Temporal domain.
        xLim : tuple or array-like
            Spatial domain.
        xDtc : float
            Detector position
        """
        t = np.linspace(tLim[0], tLim[1], numPoints)
        x = np.linspace(xLim[0], xLim[1], numPoints)
        T, X = np.meshgrid(t, x)
        density = np.abs(self.superposition(T, X))**2
        densityDtc = np.abs(self.superposition(t, xDtc))**2

        plt.rcParams.update({'font.size': 12})
        plt.rcParams["font.family"] = "serif"
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi = 400,\
                                            constrained_layout = True) # Figure
        pcm = ax1.pcolormesh(T, X, density, shading='auto', cmap='viridis')
        fig.colorbar(pcm, ax = ax1, label="Density")
        ax1.set_xlabel("Time")
        ax1.set_ylabel("Position")
        ax1.title.set_text("Density")
        ax2.plot(t, densityDtc, color = 'tab:blue', marker = '.', linestyle = '-',\
                 markersize = 1, linewidth = 0.5)
        ax2.set_xlabel("Time")
        ax2.title.set_text(f"Density at x = {xDtc}")
        
        plt.show()