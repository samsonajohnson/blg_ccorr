import numpy as np
import matplotlib.pyplot as plt
import scipy 
import scipy.interpolate
from astropy.io import fits
import ipdb
import pickle
import dill

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
#        ipdb.set_trace()
        for order in range(self.data.shape[0]):
            coeffs = np.polyfit(self.wavelens[order],self.data[order]\
                                    /np.sum(self.data[order]),5)
            temp_cs_arr = (self.data[order]/np.sum(self.data[order])\
                               -(coeffs[5]\
                                     +coeffs[4]*self.wavelens[order]\
                                     +coeffs[3]*self.wavelens[order]**2\
                                     +coeffs[2]*self.wavelens[order]**3\
                                     +coeffs[1]*self.wavelens[order]**4\
                                     +coeffs[0]*self.wavelens[order]**5))
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

                    
        
        
                                
                                

def corr(spec1, spec2, min=None, max=None):
    conv = np.convolve(spec1,spec2,mode='same')
    squaresum1 = np.sum(spec1**2)
    squaresum2 = np.sum(spec2**2)
    return conv/np.sqrt(squaresum1*squaresum2)
    pass


if __name__ == '__main__':
    #next steps
    # do derivation, figure the velocity conversion
    # loop 
    # read paper, ultimately want to sum signal over all orders

    #Only doing order with Ca lines right now
    numpoints = 20000.

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
                dillfile = open('./'+star.name+'_norm.pkl','wb')
                dill.dump(star,dillfile)
                dillfile.close()
    else:
        suffix = '_maxdiv'
        blg = dill.load(open('./BLG0966'+suffix+'.pkl','rb'))
        hr = dill.load(open('./HR4963'+suffix+'.pkl','rb'))
        hd = dill.load(open('./HD142527'+suffix+'.pkl','rb'))

            
    ipdb.set_trace()
    #set star1 and star2 
    star1 = blg
    # the star with the smaller spectra
    star2 = hr

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
        
#        print (np.log(blg.wavelens[order].max())-\
#                   np.log(blg.wavelens[order].min()))/2048
#        print  (exppowers_arr[order].max()-exppowers_arr[order].min())\
#            /len(exppowers_arr[order])*3e8


    exppowers_arr = np.array(exppowers_arr)
    expspace_arr = np.array(expspace_arr)

#    print (np.log(blg.wavelens[order].max())-np.log(blg.wavelens[order].min()))/2048
#    print 'our step size'
    lnstep= (exppowers_arr[3].max()-exppowers_arr[3].min())/len(exppowers_arr[3])
    vstep = 3e8*lnstep
#    print lnstep,vstep



    


    

    autocorr_arr = []
    ccorr_arr = []

    #going to ignore the 100 points on either end to avoid these 
    #bad points
    edgnore = 50
    indmin = 499
    indmax = 19499

    ipdb.set_trace()
    for order in range(star1.data.shape[0]):
 
        autocorr = corr(\
            star1.cs_splines[order](expspace_arr[order]\
                                        [edgnore+indmin:indmax-edgnore+1]),\
                star1.cs_splines[order](expspace_arr[order]\
                                            [edgnore:-edgnore+1])[::-1])

        ccorr = corr(\
            star1.cs_splines[order](expspace_arr[order]\
                                        [edgnore+indmin:indmax-edgnore+1]),\
                star2.cs_splines[order](expspace_arr[order]\
                                            [edgnore:-edgnore+1])[::-1])
#        autocorr = np.correlate(\
#            star1.cs_splines[order](expspace_arr[order][edgnore+indmin:\
#                                             indmax-edgnore+1]),\
#                star1.cs_splines[order](expspace_arr[order][edgnore:\
#                                                 -edgnore]),mode='same')
#        ccorr = np.correlate(\
#            star1.cs_splines[order](expspace_arr[order][edgnore+indmin:\
#                                             indmax-edgnore+1]),\
#                star2.cs_splines[order](expspace_arr[order][edgnore:\
#                                                 -edgnore]),mode='same')
        
        autocorr_arr.append(autocorr)
        ccorr_arr.append(ccorr)
    
    #attempt at zucker
    ipdb.set_trace()
    sum = np.zeros(len(ccorr_arr[0]))
    for order in range(star1.data.shape[0]):
        sum += ccorr_arr[order]#/np.median(ccorr_arr[order])
#        plt.plot((np.arange(len(sum))-len(sum)/2)*vstep/1000.,sum/order)
#        plt.show()
    ave = sum/float(star1.data.shape[0])
    plt.plot((np.arange(len(ave))-len(ave)/2)*vstep/1000.,ave)
    plt.xlabel('Velcoity [km/s]')
    plt.show()
    ipdb.set_trace()

    prod = np.ones(len(ccorr_arr[0]))
    ipdb.set_trace()
    for order in range(star1.data.shape[0]):
        print np.max(prod)
        prod *= (1. - (ccorr_arr[order])**2)
    maxlike = 1. - (prod**(1./star1.data.shape[0]))

    plt.plot((np.arange(len(ave))-len(ave)/2)*vstep/1000.,maxlike)
    plt.xlabel('Velcoity [km/s]')
    plt.show()
    
    ipdb.set_trace()
    for order in (np.arange(star1.data.shape[0]/2)):#+star1.data.shape[0]/2):
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


