#Just need for spectra classes
import corrspline as corr
#For rebinning of models
import model_bin as mb
import ipdb
import numpy as np
import matplotlib.pyplot as plt
import dill
import scipy
from scipy import optimize


def fit_model(p,full_model,wavelens):
    #S make a chunk of the model shifted by doppler factor z
    model = mb.john_rebin(full_model.wavelens,full_model.data,wavelens,p[0])
    #S continuum subtraction
    coeffs = np.polyfit(wavelens,model,5)
    #S do we need to add 1? will just come out in the wash. might want to boost
    #S the observed though
    cs_model = model - np.polyval(coeffs,wavelens) + 1.
    return p[1]*model/np.sum(model) + p[2]


def err_func(p,model,wavelens,spec):
    #S i think the fitter minimizes the square of this, no need to square here?
    return (fit_model(p,model,wavelens) - spec)
    


if __name__ == '__main__':
    full_model = corr.highres_spec('./t05500_g+0.5_m10p04_hr.fits')
    suffix = '_norm2'
    blg = dill.load(open('./BLG0966'+suffix+'.pkl','rb'))
    hr = dill.load(open('./HR4963'+suffix+'.pkl','rb'))
    hd = dill.load(open('./HD142527'+suffix+'.pkl','rb'))
    ipdb.set_trace()
    order = 3
    ind = np.where(blg.data[order]!=0)[0]
    p0 = [0.0001666,250000.,46.]
#    p0 = [0.0001666,2500000./2,-500.]#4
    out=scipy.optimize.leastsq(err_func,p0,args=\
                                   (full_model,blg.wavelens[order],\
                                        blg.data[order])\
                                   ,full_output=1)
    p1 = out[0]
    covar=out[1]
    mydict=out[2]
    message=out[3]
    ier=out[4]
    resids=mydict['fvec']
    chisq = np.sum(resids**2)
    degs_frdm=len(blg.wavelens[order])-len(p1)
    red_chisq = chisq/degs_frdm
    print 'Fitter status:',ier,' Message: ',message
    i=0
    for u in p1:
        print 'Param: ', i+1, ': ',u,' +/-',np.sqrt(covar[i,i])
        i+=1
    print 'Chisq: ',chisq,' Reduced Chisq: ',red_chisq


    plt.plot(blg.wavelens[order],blg.data[order],zorder=1)
    plt.plot(blg.wavelens[order],fit_model(p1,full_model,blg.wavelens[order]),\
                 'r',zorder=2)
    plt.show()
    ipdb.set_trace()
