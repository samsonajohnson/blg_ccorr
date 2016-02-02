import numpy as np
import matplotlib.pyplot as plt
import scipy 
import scipy.interpolate
from astropy.io import fits
import ipdb

#ipdb.set_trace()
# Only doing order with Ca lines right now
numpoints = 2000
min = 8475e-10
max = 8700e-10
blgmax = 8600e-10
blg = fits.open('./blg0966red_multi.fits')[0].data[6,4,:]
hr = fits.open('./hr4963red_multi.fits')[0].data[6,4,:]
# Powers for exponentially spaced points 
exppowers = np.linspace(np.log(min),np.log(max),num=numpoints,dtype='Float64')
expspace = np.exp(exppowers)
bexppowers = np.linspace(np.log(min),np.log(max),num=numpoints,dtype='Float64')
bexpspace = np.exp(bexppowers)

blgstart = 8456.28384384e-10
blgstep = 0.12100745e-10
blgwave = np.linspace(blgstart,blgstart+2047.*blgstep,2048)
hrstart = 8456.30225620e-10
hrstep = 0.12100766e-10
hrwave = np.linspace(hrstart,hrstart+2047.*hrstep,2048)

blgspline = scipy.interpolate.interp1d(blgwave,blg,kind='cubic')
hrspline = scipy.interpolate.interp1d(hrwave,hr,kind='cubic')

ipdb.set_trace()
autocorr = np.correlate(blgspline(expspace),blgspline(expspace))
ccorr = np.correlate(blgspline(expspace),hrspline(expspace))
plt.plot(np.arange(len(autocorr)),autocorr,'b')
plt.show()
plt.plot(np.arange(len(ccorr)),ccorr,'b')
plt.show()
ipdb.set_trace()

