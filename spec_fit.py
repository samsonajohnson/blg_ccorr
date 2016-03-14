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


def cfit_model(p,full_model,wavelens):
    #S make a chunk of the model shifted by doppler factor z
    model = mb.john_rebin(full_model.wavelens,full_model.data,wavelens,p[0])
    #S continuum subtraction
#    coeffs = np.polyfit(wavelens,model,5)
    #S do we need to add 1? will just come out in the wash. might want to boost
    #S the observed though
#    cs_model = model - np.polyval(coeffs,wavelens) 
    convmodel = mb.numconv(model,gaussian(p))
    factor = np.polyval(p[:1:-1],wavelens)
    return factor*convmodel/np.max(convmodel) 

def fit_model(p,full_model,wavelens):
    #S make a chunk of the model shifted by doppler factor z
    model = mb.john_rebin(full_model.wavelens,full_model.data,wavelens,p[0])
    #S continuum subtraction
    coeffs = np.polyfit(wavelens,model,5)
    #S do we need to add 1? will just come out in the wash. might want to boost
    #S the observed though
    cs_model = model - np.polyval(coeffs,wavelens) 
#    return p[1]*model/np.max(model) 
    factor = np.polyval(p[:0:-1],wavelens)
    return factor*model/np.max(model) 

def gaussian(p):
    xip = np.arange(-10,10+.25,.25)
    return np.exp(-(xip/p[1])**2.)

def err_func(p,model,wavelens,spec):
    #S i think the fitter minimizes the square of this, no need to square here?
#    conv_model = mb.numconv(fit_model(p,model,wavelens),gaussian(p))
    #S to include sigma clipped points, and points greater than zero
    zeroinds = np.where(spec>0)[0]
    med = np.median(spec[zeroinds])
    std = np.std(spec[zeroinds])
    siginds = np.where(np.abs(spec[zeroinds]-med)<3*std)[0]
    inds = zeroinds[siginds]
#    inds = zeroinds
    fitspec = spec[inds]
#    ipdb.set_trace()
#    return (cfit_model(p,model,wavelens)-\
#                fitspec)/np.sqrt(np.abs(fitspec))
    return (cfit_model(p,model,wavelens)[inds] -\
                fitspec)/np.sqrt(np.abs(fitspec))
    


if __name__ == '__main__':
#    full_model = corr.highres_spec('./t05500_g+0.5_m10p04_hr.fits')
    full_model = corr.highres_spec('./t05500_g+4.0_p00p00_hrplc.fits')
    suffix = '_norm2'
    blg = dill.load(open('./BLG0966'+suffix+'.pkl','rb'))
    hr = dill.load(open('./HR4963'+suffix+'.pkl','rb'))
    hd = dill.load(open('./HD142527'+suffix+'.pkl','rb'))
#    ipdb.set_trace()

    rv=[]
    rver = []
    """
    plt.plot(np.arange(len(blg.data)),np.sqrt(np.median(blg.data,axis=1)),'o')
    plt.xlabel('Order')
    plt.ylabel(r'$\sqrt{{\bf MEDIAN}(N)}$')
    plt.show()
    ipdb.set_trace()
    """
    ipdb.set_trace()
    offset = 250
    params = []
    for order in np.arange(len(blg.data)-28)+2:

#    while True:
#        order = int(raw_input('Order to fit: '))

        ind = np.where(blg.data[order]>0)[0]        

#        coeffs = np.polyfit(blg.wavelens[order][ind],blg.data[order][ind],2,\
#                                w=1/np.sqrt(blg.data[order][ind]))
#        cscorr = np.polyval(coeffs,blg.wavelens[order])
#        fit_data = (blg.data[0]-cscorr+median)

        p0 = [0.00014,1.,np.median(blg.data[order]),0.]#,0.]
#        p0 = [0.00014,np.median(blg.data[order]),0.]#,0.]

        out=scipy.optimize.leastsq(err_func,p0,args=\
                                       (full_model,blg.wavelens[order],\
#                                            fit_data),\
                                        blg.data[order]),\
                                        full_output=1,maxfev=1000)

        message=out[3]
        ier=out[4]
        print 'Fitter status:',ier,' Message: ',message        
        p1 = out[0]
        if out[1] == None:
            ipdb.set_trace()
            continue
        jacob=out[1]
        mydict=out[2]

        resids=mydict['fvec']
#        covar = np.std(resids)**2*jacob
        chisq = np.sum(resids**2)
        degs_frdm=len(blg.wavelens[order])-len(p1)
        red_chisq = chisq/degs_frdm

        i=0
        for u in p1:
            print 'Param: ', i+1, ': ',u,' +/-'#,np.sqrt(covar[i,i])
            i+=1
        print 'Chisq: ',chisq,' Reduced Chisq: ',red_chisq
        params.append(p1)
        rv.append(p1[0]*2.99e8)
#        rver.append(np.sqrt(covar[0,0])*3e8)
#    ipdb.set_trace()

        plt.plot(np.arange(len(blg.data[0])),offset*order+\
                     blg.data[order],'b',zorder=1)
#        plt.plot(np.arange(len(blg.data[0])),\
#                     np.polyval(coeffs,blg.wavelens[order])-coeffs[0]+\
#                                    offset*order,'y',zorder=2)
        plt.plot(np.arange(len(blg.data[0])),offset*order+\
                     cfit_model(p0,full_model,blg.wavelens[order]),\
                     'g',zorder=2)
        plt.plot(np.arange(len(blg.data[0])),offset*order+\
                     cfit_model(p1,full_model,blg.wavelens[order]),\
                     'r',zorder=2)

        """
        plt.plot(blg.wavelens[order],blg.data[order]+offset*order,'b',zorder=1)
        plt.plot(blg.wavelens[order],offset*order+fit_model(p0,full_model,blg.wavelens[order]),\
                     'g',zorder=2)
        plt.plot(blg.wavelens[order],offset*order+fit_model(p1,full_model,blg.wavelens[order]),\
                     'r',zorder=2)
        
        """
    plt.xlabel('Bin')
    plt.ylabel('Arbitrary')
    plt.show()
    ipdb.set_trace()
    rv = np.array(rv)/1000.
    rver = np.array(rver)/1000.
    plt.plot(np.arange(len(rv)),rv,'o')
    
    plt.errorbar(np.arange(len(rv)),rv,yerr=rver)
    plt.show()
    ipdb.set_trace()
