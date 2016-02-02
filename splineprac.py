import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate
from astropy.io import fits
import ipdb

ipdb.set_trace()
data = fits.open('./blg0966red_multi.fits')[0].data[6,:,:]
seg = data[0,:]
#note that iraf gets wavelength info from header, could do that here
x = np.arange(seg.shape[0])

newx = np.linspace(x.min(),x.max(),2*x.shape[0])

pltmin = 0
pltmax = 2047
#linint = scipy.interpolate.interp1d(x[pltmin:pltmax],seg[pltmin:pltmax],kind='linear')
#spline = scipy.interpolate.interp1d(x[pltmin:pltmax],seg[pltmin:pltmax],kind='cubic')

#plt.figure()
#plt.subplot(131)
plt.plot(x[pltmin:pltmax],seg[pltmin:pltmax])
plt.axis([pltmin-20,pltmax+20,-20,seg.max()+100])
plt.show()
plt.subplot(132)
plt.plot(newx[2*(pltmin+1):2*(pltmax-1)],linint(newx[2*(pltmin+1):2*(pltmax-1)]))
plt.axis([pltmin-20,pltmax+20,-20,600])
plt.subplot(133)
resids = spline(newx[2*(pltmin+1):2*(pltmax-1):2]) - seg[pltmin:-2]
print resids.shape
plt.plot(x[pltmin:pltmax],resids[pltmin:pltmax],'r.')
plt.axis([pltmin-20,pltmax+20,resids.min(),resids.max()])
plt.show()



