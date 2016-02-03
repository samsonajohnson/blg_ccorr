import numpy as np
import matplotlib.pyplot as plt
import scipy 
import scipy.interpolate
from astropy.io import fits
import ipdb
import pickle


class multi_spec:
    def __init__(self,fits_name,pixels = 2048):
        self.hdu = fits.open(fits_name)
        self.spec_string = ''
        self.name = self.hdu[0].header['OBJECT']
        for key in self.hdu[0].header.keys():
            if 'WAT2' in key:
                #there is an issue in my method here of reading in the header
                #that if a line ends or starts with a space character (' '),
                #it is ignored. i'm not really sure why, but i think this work
                #work around will do.
                if len(self.hdu[0].header[key]) == 67:
                    self.spec_string += self.hdu[0].header[key] + ' '
                else:
                    self.spec_string += self.hdu[0].header[key] 
        self.data = self.hdu[0].data[6,:,:]
        self.cs_data=[]
        self.wavelens = []
        
        #dredge the wavelegnth ranges from the spectrum header
        for order in range(self.data.shape[0]):
            try:
                start = float(self.spec_string.split('spec')[order+2].\
                                  split()[5])
                step = float(self.spec_string.split('spec')[order+2].\
                                 split()[6])
            except:
                ipdb.set_trace()
            self.wavelens.append(np.linspace(start,start+2047.*step,pixels))
        self.wavelens = np.array(self.wavelens)

        #continuum subtract the data, fitting a cubic
        for order in range(self.data.shape[0]):
            coeffs = np.polyfit(self.wavelens[order],self.data[order],5)
            self.cs_data.append(self.data[order]\
                                    -(coeffs[5]\
                                          +coeffs[4]*self.wavelens[order]\
                                          +coeffs[3]*self.wavelens[order]**2\
                                          +coeffs[2]*self.wavelens[order]**3\
                                          +coeffs[1]*self.wavelens[order]**4\
                                          +coeffs[0]*self.wavelens[order]**5))
            #a catch to see if the continuum subtraction is still returning 
            #weird results
            if False in self.cs_data[order]==self.cs_data[order]:
                ipdb.set_trace()
            if float('Inf') in self.cs_data[order]:
                ipdb.set_trace()
        self.cs_data = np.array(self.cs_data)
        
        #define splines list, should probably go up higher
        self.splines = []
#        for order in range(self.data.shape[0]):
#            self.spline(order)
#        ipdb.set_trace()

    def spline(self,order):
        self.splines.append(scipy.interpolate.interp1d(\
                self.wavelens[order],\
                    self.cs_data[order],\
                    kind = 'cubic'))

    def pickle(self):
        f = file('./test.pk','wb')
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def unpickle(self):
        with file('test_file', 'rb') as f:
            return pickle.load(f)

                    
        
        
                                
                                

def corr(spec1, spec2, min=None, max=None):
    return np.correlate(spline1(expspace[indmin:indmax]),\
                            spline2(expspace),mode='same')


if __name__ == '__main__':
    #next steps
    # do derivation, figure the velocity conversion
    # loop 
    # read paper, ultimately want to sum signal over all orders

    #Only doing order with Ca lines right now
    numpoints = 20000.


    #read in the fits files from the stars we want
    blg = multi_spec('./blg0966red_multi.fits')
    hr = multi_spec('./hr4963red_multi.fits')
    hd = multi_spec('./hd142527red_multi.fits')

    #set star1 and star2 
    star1 = blg
    star2 = hd

    #initializing the arrays(lists at first)
    exppowers_arr = []
    expspace_arr = []

    #for all the orders
    for order in range(star1.data.shape[0]):
        #we need to be sure the range of wavelengths we sample is included
        #in the spline of the order from either star, e.g. it tries to
        #evaluate the spline at a point outside the range, which is a nogo

        #find the minimum of the lowest wavelengths in this order from 
        #either star
        minwl = np.max([star1.wavelens[order][0],star2.wavelens[order][0]])\
            +.00001
        #now the bigger of the minimums
        if star1.name == 'HR4963' and order == 3:
           ipdb.set_trace()
        maxwl = np.min([star1.wavelens[order][-1],star2.wavelens[order][-1]])\
            -.00001

        #linspace between the log of min/max
        exppowers_arr.append(np.linspace(np.log(minwl),np.log(maxwl),\
                                             num=numpoints,dtype='Float64'))
        #now back to no-log space, making these sampling points exponentially
        #spaced
        expspace_arr.append(np.exp(exppowers_arr[order]))
        
        print (np.log(blg.wavelens[order].max())-\
                   np.log(blg.wavelens[order].min()))/2048
        print  (exppowers_arr[order].max()-exppowers_arr[order].min())\
            /len(exppowers_arr[order])*3e8
    ipdb.set_trace()
    #convert to numpy arrays, easiser to deal with
    exppowers_arr = np.array(exppowers_arr)
    expspace_arr = np.array(expspace_arr)

#    print (np.log(blg.wavelens[order].max())-np.log(blg.wavelens[order].min()))/2048
#    print 'our step size'
    lnstep= (exppowers_arr[3].max()-exppowers_arr[3].min())/len(exppowers_arr[3])
    vstep = 3e8*lnstep
#    print lnstep,vstep


    t_order = 20

#    blgspline = scipy.interpolate.interp1d(blg.wavelens[t_order],blg.cs_data[t_order],kind='cubic')
#    hrspline = scipy.interpolate.interp1d(hr.wavelens[t_order],hr.cs_data[t_order],kind='cubic')
#    hdspline = scipy.interpolate.interp1d(hd.wavelens[t_order],hd.cs_data[t_order],kind='cubic')


    s1spline = scipy.interpolate.interp1d(star1.wavelens[t_order],star1.cs_data[t_order],kind='cubic')
    s2spline = scipy.interpolate.interp1d(star2.wavelens[t_order],star2.cs_data[t_order],kind='cubic')

    indmin = 499
    indmax = 19499
    ipdb.set_trace()
    autocorr = np.correlate(s1spline(expspace_arr[t_order][indmin:indmax]),\
                                s1spline(expspace_arr[t_order]),mode='same')
    ccorr = np.correlate(s1spline(expspace_arr[t_order][indmin:indmax]),\
                             s2spline(expspace_arr[t_order]),mode='same')
    plt.plot((np.arange(len(autocorr))-len(autocorr)/2)*vstep/1000.,autocorr,'b')
    plt.xlabel('Velcoity [km/s]')
    plt.show()
    plt.plot((np.arange(len(ccorr))-len(ccorr)/2)*vstep/1000.,ccorr,'b')
    plt.xlabel('Velcoity [km/s]')
    pltmax = np.max(ccorr[9499:10499])+.1*np.max(ccorr[9499:10499])
    pltmin = np.min(ccorr[9499:10499])-.1*np.min(ccorr[9499:10499])
    plt.axis([-500*vstep/1000, 500*vstep/1000,pltmin,pltmax])
    plt.show()
    print autocorr.argmax()
    print ccorr.argmax()
    ipdb.set_trace()


