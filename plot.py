import numpy as np
import corrspline as corr
#For rebinning of models                                                        
import model_bin as mb
import matplotlib
import ipdb
import matplotlib.pyplot as plt
import dill
import scipy
from scipy import optimize


if __name__ == '__main__':
    # get the full highres model spectrum class
    full_model = corr.highres_spec('./t05500_g+4.0_p00p00_hrplc.fits')
    # the blg spectrum class
    blg = corr.multi_spec('./blg0966red_multi.fits')
    # the simultaneous fit parameters, z and ip width (just a gaussian)
    p0 = [0.0001,1.]
    # the empty dictionary for putting the trimmed indices
    blg.inds={}
    # the pixels to sjip at the beginnning of each order
    skipf = 48
    # the length of chunk we are using
    cklen = 400
    # empty dict for data
    blg.fit_data = {}
    # list of orders we want to fit
    blg.fit_orders = [2,3,]#4,5,6,7,8,9,10,11,12]#np.arange(len(blg.data)-19)+2                                                                              
    # number of chunks we will be fitting based on cklen and skipf
    blg.fit_chunks = np.arange((len(blg.data[0])-skipf)/cklen)

    snr1 = []
    snr2 = []
    for order in range(len(blg.data)):#blg.fit_orders:
        print np.median(blg.data[order]), np.sqrt(np.median(blg.data[order])),\
            np.median(blg.errs[order])
        snr2.append(np.median(blg.data[order])/np.median(blg.errs[order]))
        snr1.append(np.sqrt(np.median(blg.data[order])))
#        plt.plot(blg.wavelens[order],blg.data[order],'.',zorder=2)
#        plt.errorbar(blg.wavelens[order],blg.data[order],blg.errs[order],\
#                         zorder=1,fmt='')
#    plt.plot(snr1,'o',label='$\sqrt{ median(data)}$')
    plt.plot(snr2,'ko',ms=8)#,label='$median(data)/median(noise)$')
    plt.xlabel('Order')
    plt.ylabel('SNR')
    plt.grid(True)
    matplotlib.rcParams.update({'font.size': 22})
    plt.show()
