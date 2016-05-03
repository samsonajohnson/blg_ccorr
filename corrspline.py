import numpy as np
import matplotlib.pyplot as plt
import scipy 
import scipy.interpolate
from astropy.io import fits
import ipdb
import pickle
import dill
import pysynphot
#import model_bin as mb
#import spec_fit as sf

class highres_spec:
    def __init__(self,fits_name):
        self.hdu = fits.open(fits_name)
        self.data = self.hdu[0].data
        #S hard coding in the wavelengths
        self.wavelens = np.linspace(2500.,9000.,325001)
    

    def corr_set(self,wavelens):
        #S returns a spline sampled at the given wavelengths
        model = mb.john_rebin(self.wavelens,self.data,wavelens,0.)
        coeffs = np.polyfit(wavelens,model,5)
        cs_model = model - np.polyval(coeffs,wavelens)
        spline = scipy.interpolate.interp1d(\
                wavelens,\
                    cs_data,\
                    kind = 'cubic')    
        return spline


class phe_spec:
    def __init__(self,specfits,wavefits,minwave=False,maxwave=False):
        spec_hdu = fits.open(specfits)
        self.data = spec_hdu[0].data
        wave_hdu = fits.open(wavefits)
        vac_waves = wave_hdu[0].data
        sigma_2 = (1.e4/vac_waves)**2
        f = 1.0 + 0.05792105/(238.0185-sigma_2)+0.00167917/(57.362-sigma_2)
#        ipdb.set_trace()
        self.wavelens = vac_waves/f
        if minwave and maxwave:

            low_inds = np.where(self.wavelens>minwave)[0]
            high_inds = np.where(self.wavelens[low_inds]<maxwave)[0]
            minds = low_inds[high_inds]
            self.wavelens = self.wavelens[minds]
            self.data = self.data[minds]

        
        
        

class multi_spec:
    def __init__(self,fits_name,pixels = 2048,spline=False):
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
        self.errs = self.hdu[0].data[2,:,:]
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
        if not spline:
            return
        #continuum subtract the data, fitting a cubic
#        ipdb.set_trace()
        for order in range(self.data.shape[0]):
            coeffs = np.polyfit(self.wavelens[order],\
                                    self.data[order]/np.sum(self.data[order]),\
                                    5)
            temp_cs_arr = (self.data[order]/np.sum(self.data[order])-\
                               np.polyval(coeffs,self.wavelens[order]))
#
#                               -(coeffs[5]\
#                                     +coeffs[4]*self.wavelens[order]\
#                                     +coeffs[3]*self.wavelens[order]**2\
#                                     +coeffs[2]*self.wavelens[order]**3\
#                                     +coeffs[1]*self.wavelens[order]**4\
#                                     +coeffs[0]*self.wavelens[order]**5))
            self.cs_data.append(temp_cs_arr)
                                          
            #a catch to see if the continuum subtraction is still returning 
            #weird results
            if np.isnan(np.min(self.cs_data[order])):
                ipdb.set_trace()
            if np.isinf(np.max(self.cs_data[order])) or \
                    np.isinf(np.min(self.cs_data[order])):
                ipdb.set_trace()
        self.cs_data = np.array(self.cs_data)
        
        #define splines list, should probably go up higher
        self.splines = []
        self.cs_splines = []
        for order in range(self.data.shape[0]):
            print('Working on spline '+str(1+order)+' of '+self.name)
            self.splines.append(self.spline(order))
            self.cs_splines.append(self.cs_spline(order))




    def spline(self,order):
        return scipy.interpolate.interp1d(\
                self.wavelens[order],\
                    self.data[order],\
                    kind = 'cubic')

    def cs_spline(self,order):
        return scipy.interpolate.interp1d(\
                self.wavelens[order],\
                    self.cs_data[order],\
                    kind = 'cubic')

    def pickle(self):
        f = file('./test.pk','wb')
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def unpickle(self):
        with file('test_file', 'rb') as f:
            return pickle.load(f)

                    
        
        
                                
                                

def mycorr(spec1, spec2, min=None, max=None):
    conv = np.convolve(spec1,spec2,mode='same')
    squaresum1 = np.sum(spec1**2)
    squaresum2 = np.sum(spec2**2)
    return conv/np.sqrt(squaresum1*squaresum2)
    pass


if __name__ == '__main__':
    #next steps
    # do derivation, figure the velocity conversion
    # loop 

    numpoints = 2048.
    edgnore = 1
    indmin = 50
    indmax = 50

    use_prev = raw_input('Load past multi_specs? (y or n)\n')
    #read in the fits files from the stars we want
    if use_prev == 'n':
        blg = multi_spec('./blg0966red_multi.fits')
        hr = multi_spec('./hr4963red_multi.fits')
        hd = multi_spec('./hd142527red_multi.fits')
        save_new = raw_input('save these mutli_specs?(y or n)\n')
        if save_new == 'y':
            ipdb.set_trace()
            for star in [blg,hr,hd]:
                dillfile = open('./'+star.name+'_norm2.pkl','wb')
                dill.dump(star,dillfile)
                dillfile.close()
    else:
        suffix = '_maxdiv'
        blg = dill.load(open('./BLG0966'+suffix+'.pkl','rb'))
        hr = dill.load(open('./HR4963'+suffix+'.pkl','rb'))
        hd = dill.load(open('./HD142527'+suffix+'.pkl','rb'))



            
    ipdb.set_trace()
    #set star1 and star2 
    # the star with the smaller spectra
    star1 = hr
    # the reference?
    star2 = blg

    #Allen's astro quantities, or similar
    #

    #initializing the arrays(lists at first)
    exppowers_arr = []
    expspace_arr = []

    #for all the orders
    for order in range(star1.data.shape[0]):
        #we need to be sure the range of wavelengths we sample is included
        #in the spline of the order from either star, e.g. it tries to
        #evaluate the spline at a point outside the range, which is a nogo

        #find the maximum of the lowest wavelengths in this order from 
        #either star
        minwl = np.max([star1.wavelens[order][0],star2.wavelens[order][0]])\
            +.00001
        #now the smaller of the maximums
        maxwl = np.min([star1.wavelens[order][-1],star2.wavelens[order][-1]])\
            -.00001

        #linspace between the log of min/max
        exppowers_arr.append(np.linspace(np.log(minwl),np.log(maxwl),\
                                             num=numpoints,dtype='Float64'))
        #now back to no-log space, making these sampling points exponentially
        #spaced
        expspace_arr.append(np.exp(exppowers_arr[order]))
        lnstep= (exppowers_arr[order].max()-exppowers_arr[order].min())/len(exppowers_arr[order])
        vstep = 3e8*lnstep
        print lnstep,vstep
        
#        print (np.log(blg.wavelens[order].max())-\
#                   np.log(blg.wavelens[order].min()))/2048
#        print  (exppowers_arr[order].max()-exppowers_arr[order].min())\
#            /len(exppowers_arr[order])*3e8


    exppowers_arr = np.array(exppowers_arr)
    expspace_arr = np.array(expspace_arr)

#    print (np.log(blg.wavelens[order].max())-np.log(blg.wavelens[order].min()))/2048
#    print 'our step size'
#    lnstep= (exppowers_arr[3].max()-exppowers_arr[3].min())/len(exppowers_arr[3])
#    vstep = 3e8*lnstep
#    print lnstep,vstep



    


    

    autocorr_arr = []
    ccorr_arr = []

    #going to ignore the 100 points on either end to avoid these 
    #bad points


    ipdb.set_trace()
    


    for order in range(star1.data.shape[0]):
        print 'working on corrs for order: '+str(order+1)
        autocorr = mycorr(\
            star1.cs_splines[order](expspace_arr[order]\
                                        [edgnore+indmin:-indmax-edgnore]),\
                star1.cs_splines[order](expspace_arr[order]\
                                            [edgnore:-edgnore])[::-1])

        ccorr = mycorr(\
            star1.cs_splines[order](expspace_arr[order]\
                                        [edgnore+indmin:-indmax-edgnore]),\
                star2.cs_splines[order](expspace_arr[order]\
                                            [edgnore:-edgnore])[::-1])
        
        autocorr_arr.append(autocorr)
        ccorr_arr.append(ccorr)
    
    #attempt at zucker
    ipdb.set_trace()
    sum = np.zeros(len(ccorr_arr[0]))
    for order in range(star1.data.shape[0]):
        sum += ccorr_arr[order]
    ave = sum/float(star1.data.shape[0])
    peak = (ave.argmax()-len(ave)/2)*vstep/1000.
    plt.plot((np.arange(len(ave))-len(ave)/2)*vstep/1000.,ave,\
                 label='peak: %.3f km/s'%(peak))
    plt.xlabel('Velcoity [km/s]')
    plt.title('Simple Average: '+star1.name+' X '+star2.name)
    pltindmin = len(ave)/2-indmin
    pltindmax = len(ave)/2+indmax
    pltmax = np.max(ave[pltindmin:pltindmax])+.1*\
        np.max(ave[pltindmin:pltindmax])
    pltmin = np.min(ave[pltindmin:pltindmax])-.1*\
        np.min(ave[pltindmin:pltindmax])
    plt.axis([-indmin*vstep/1000, indmax*vstep/1000,pltmin,pltmax])
    plt.legend()
    plt.show()
    ipdb.set_trace()

    prod = np.ones(len(ccorr_arr[0]))
    ipdb.set_trace()
    for order in range(star1.data.shape[0]):
        prod *= (1. - (ccorr_arr[order])**2)
    maxlike = 1. - (prod**(1./star1.data.shape[0]))

    peak = (maxlike.argmax()-len(maxlike)/2)*vstep/1000.
    plt.plot((np.arange(len(maxlike))-len(maxlike)/2)*vstep/1000.,maxlike,\
                 label='peak: %.3f km/s'%(peak))
    plt.title('Maximum Likelihood: '+star1.name+' X '+star2.name)
    plt.xlabel('Velcoity [km/s]')

    pltindmin = len(maxlike)/2-indmin
    pltindmax = len(maxlike)/2+indmax
    pltmax = np.max(maxlike[pltindmin:pltindmax])+.1*\
        np.max(maxlike[pltindmin:pltindmax])
    pltmin = np.min(maxlike[pltindmin:pltindmax])-.1*\
        np.min(maxlike[pltindmin:pltindmax])

    plt.axis([-indmin*vstep/1000, indmax*vstep/1000,pltmin,pltmax])
    plt.legend()
    plt.show()
    
    ipdb.set_trace()
    for order in (np.arange(star1.data.shape[0]/2)):
        plt.plot((np.arange(len(autocorr_arr[order]))\
                      -len(autocorr_arr[order])/2)*vstep/1000.,\
                     autocorr_arr[order],label=order)
    plt.xlabel('Velcoity [km/s]')
#    pltmax = np.max(ccorr[9499:10499])+.1*np.max(ccorr[9499:10499])
#    pltmin = np.min(ccorr[9499:10499])-.1*np.min(ccorr[9499:10499])
#    plt.axis([-500*vstep/1000, 500*vstep/1000,pltmin,pltmax])
    plt.legend()
    plt.show()
    for order in (np.arange(star1.data.shape[0]/2)):#+star1.data.shape[0]/2):
        plt.plot((np.arange(len(ccorr_arr[order]))\
                      -len(ccorr_arr[order])/2)*vstep/1000.,\
                     ccorr_arr[order],label=order)
    plt.xlabel('Velcoity [km/s]')
    plt.legend()
#    pltmax = np.max(ccorr[9499:10499])+.1*np.max(ccorr[9499:10499])
#    pltmin = np.min(ccorr[9499:10499])-.1*np.min(ccorr[9499:10499])
#    plt.axis([-500*vstep/1000, 500*vstep/1000,pltmin,pltmax])
#    plt.axis
    plt.show()
        
    """
    use_prev = raw_input('Load existing?(y or n)\n')
    if use_prev == 'n':
        for order in range(star1.data.shape[0]):
            print('Working on order '+str(1+order)+'...\n')
            #this is where all the time is going on. these splines! tried to 
            #figure a way to save and recall them but have had no luck
            s1spline = scipy.interpolate.interp1d(\
                star1.wavelens[order],star1.cs_data[order],kind='cubic')
            s2spline = scipy.interpolate.interp1d(\
                star2.wavelens[order],star2.cs_data[order],kind='cubic')
    
            #going to ignore the 100 points on either end to avoid these 
            #bad points
            edgnore = 50
            indmin = 499
            indmax = 19499
            #ipdb.set_trace()
            autocorr = np.correlate(\
                s1spline(expspace_arr[order][edgnore+indmin:\
                                                   indmax-edgnore+1]),\
                    s1spline(expspace_arr[order][edgnore:\
                                                       -edgnore]),mode='same')
            ccorr = np.correlate(\
                s1spline(expspace_arr[order][edgnore+indmin:\
                                                   indmax-edgnore+1]),\
                    s2spline(expspace_arr[order][edgnore:\
                                                       -edgnore]),mode='same')

            autocorr_arr.append(autocorr)
            ccorr_arr.append(ccorr)

#            autocorr_arr = np.array(autocorr_arr)
#            ccorr_arr = np.array(ccorr_arr)
            #save the resulting CCF's
            np.save('./autocorr_arr_cs_data.npy',autocorr_arr)
            np.save('./ccorr_arr_cs_data.npy',ccorr_arr)
    """
    

    ipdb.set_trace()
    #if use_prev == 'y':
    #    autocorr_arr = np.load('./autocorr_arr_cs_data.npy')
    #    ccorr_arr = np.load('./ccorr_arr_cs_data.npy')

    plt.plot((np.arange(len(autocorr))-len(autocorr)/2)*vstep/1000.,autocorr,'b')
    plt.xlabel('Velcoity [km/s]')
    plt.show()
    plt.plot((np.arange(len(ccorr))-len(ccorr)/2)*vstep/1000.,ccorr,'b')
    plt.xlabel('Velcoity [km/s]')
    pltmax = np.max(ccorr[9499:10499])+.1*np.max(ccorr[9499:10499])
    pltmin = np.min(ccorr[9499:10499])-.1*np.min(ccorr[9499:10499])
#    plt.axis([-500*vstep/1000, 500*vstep/1000,pltmin,pltmax])
    plt.show()
    print autocorr.argmax()
    print ccorr.argmax()
    ipdb.set_trace()


