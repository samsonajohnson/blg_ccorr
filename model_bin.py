import numpy as np
import matplotlib.pyplot as plt
import pysynphot
from astropy.io import fits
import ipdb

import scipy.signal as sig
from scipy.interpolate import interp1d

def syn_rebin(wave,specin,wavenew):
    #S taken form astrobetter: 
    #http://www.astrobetter.com/blog/2013/08/12/python-tip-re-sampling-spectra-with-pysynphot/
    spec = pysynphot.spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = pysynphot.spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = pysynphot.observation.Observation(spec, filt, binset=wavenew, force='taper')
    return obs.binflux


def john_rebin(wold, sold, wnew, z):
    #S taken from JohnJohn
    """ Define the left and right limits of the first Wnew pixel. Keep in mind
    that pixels are labled at their centers. """    
    dw_new = wnew[1]-wnew[0]
    w_left = wnew.min() - 0.5 * dw_new
    w_right = w_left + dw_new
    Nsub = 10. # use 10 sub pixels for the template across new pixel 
    """ Create a finely sampled 'sub-grid' across a pixel. We'll later integrate
    the stellar spectrum on this fine grid for each pixel  """
    wfine_sub = np.linspace(w_left, w_right, Nsub, endpoint=False)
    Npixels = len(wnew) #Number of pixels in the new wavelength scale
    """ Copy an individual pixel subgrid Npixels times, aka 'tiling' """
    wsub_tile = np.tile(wfine_sub, Npixels)
    """ Create a separate array that ensures that each pixel subgrid is 
    properly incremented, pixel-by-pixel following the Wnew wavelength scale"""
    step = np.repeat(np.arange(Npixels), Nsub) * dw_new
    wfine = wsub_tile + step #Finely sampled new wavelength scale
    dsub = wfine[1] - wfine[0]
    wfine += dsub/2. # Center each subgrid on the pixel
    ifunction = interp1d(wold*(1+z), sold) #Create Doppler-shifted spline-interpolation function
    sfine = ifunction(wfine) #Calculate interpolated spectrum 
    sfine_blocks = sfine.reshape([Npixels,Nsub]) #Prepare for integration
    snew = np.sum(sfine_blocks, axis=1)/Nsub #Efficient, vector-based integration! 
    return snew

def numconv(y, kern):
    #S from JohnJohn again
    ynew = y.copy()
    lenex = 10
    ynew = extend(ynew, lenex)
    new = fftconvolve(ynew, kern, mode='full')
    new /= kern.sum()
    nel = len(kern)+lenex*2
    return new[nel/2.:len(new)-nel/2.+1]


def extend(arr, nel, undo=False):
    #S from minerva_dopcode.py, needed for numconv?
    new = arr.copy()
    if undo:
        return new[nel:len(new)-nel]
    else:
        new = np.append(np.append(np.zeros(nel)+arr[0], new), np.zeros(nel)+arr\
[-1])
        return new

if __name__ == '__main__':

    hdu = fits.open('./t05500_g+0.5_m10p04_hr.fits')
    data = hdu[0].data
    waves = np.linspace(2500,9000,325001)
    wavenew = np.linspace(2500.2,8999.8,32499)
    
    ipdb.set_trace()
    
    john_flux = john_rebin(waves,data,wavenew,0)
    syn_flux = syn_rebin(waves,data,wavenew)
