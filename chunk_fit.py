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
import mixmod


def cfit_model(p,full_model,wavelens):
    #S make a chunk of the model shifted by doppler factor z
    ipmodel = mb.numconv(full_model.data,gaussian(p))
    factor = np.polyval(p[:1:-1],full_model.wavelens)
    unbin_model = factor*ipmodel/np.max(ipmodel) 
#    model = mb.john_rebin(full_model.wavelens,full_model.data,wavelens,p[0])
    model = mb.john_rebin(full_model.wavelens,unbin_model,wavelens,p[0])
    return model


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
    return np.exp((-(xip)**2.)/2.036)

def err_func(p,model,star,order):
    wavelens = star.wavelens[order]
    inds = star.inds[str(order)]
    fitspec = star.data[order][inds]
    return (cfit_model(p,model,wavelens)[inds] -\
                fitspec)/np.sqrt(np.abs(fitspec))

def chunk_func(p,model,star,order,chunk):
    #S get all inds for the chunk
    ckinds = np.arange(star.cklen)+star.skipf+star.cklen*chunk
    #S get all wavelens and data for the chunk
    ckwaves = star.wavelens[order][ckinds]
    ckdata = star.data[order][ckinds]
    ckerrs = star.errs[order][ckinds]
    
    #S get the trimmed indices
    inds = star.inds['o'+str(order)+'_c'+str(chunk)]

    #S get the wavelengths in the model relevant to the order
    lowwave = ckwaves.min()-30.
    highwave = ckwaves.max()+30.
    low = np.where(full_model.wavelens > lowwave)[0]
    high = np.where(full_model.wavelens[low] < highwave)[0]
    # the relevant model inds for the order                             
    minds = low[high]
    
    # get the model wavelens and ipdata for the order                  
    modwaves = full_model.wavelens[minds]
    #S convolve the model over the relevant wavelengths
    modck = mb.numconv(full_model.data[minds],gaussian(p))
    #S add some wavelength scaling
    factor = np.polyval(p[:0:-1],modwaves)
    unbinmodel = factor*modck/np.max(modck)
    #S rebin the convolved model, and 
    tmodel = mb.john_rebin(modwaves,unbinmodel,ckwaves,p[0])

    return tmodel

    
    
def chunk_err_func(p,model,star,order,chunk):
    #S get all inds for the chunk
    ckinds = np.arange(star.cklen)+star.skipf+star.cklen*chunk
    #S get all wavelens and data for the chunk
    ckwaves = star.wavelens[order][ckinds]
    ckdata = star.data[order][ckinds]
    ckerrs = star.errs[order][ckinds]
    
    #S get the trimmed indices
    inds = star.inds['o'+str(order)+'_c'+str(chunk)]
    tmodel = chunk_func(p,model,star,order,chunk)

    #s calculate the errors
    errs = (tmodel[inds] - ckdata[inds])/(ckerrs[inds])
    return errs

if __name__ == '__main__':
#    full_model = corr.highres_spec('./t05500_g+0.5_m10p04_hr.fits')
    full_model = corr.highres_spec('./t05500_g+4.0_p00p00_hrplc.fits')
    blg = corr.multi_spec('./blg0966red_multi.fits')

    rv=[]
    rvma=[]
    rver = []
    """
    plt.plot(np.arange(len(blg.data)),np.sqrt(np.median(blg.data,axis=1)),'o')
    plt.xlabel('Order')
    plt.ylabel(r'$\sqrt{{\bf MEDIAN}(N)}$')
    plt.show()
    ipdb.set_trace()
    """
    ipdb.set_trace()
    params = []
    blg.inds = {}
    params_list=[]

    blg.skipf = 48
    blg.cklen = 400
    blg.fit_orders = np.arange(18)+2#[2,3,4,5,6,7,8,9,10,11,12]#np.arange(len(blg.data)-19)+2
    blg.fit_chunks = np.arange((len(blg.data[0])-blg.skipf)/blg.cklen)
    blg.fit_data = {}

    for order in blg.fit_orders:
        for chunk in blg.fit_chunks:
#            ipdb.set_trace()
            c_inds = np.arange(blg.cklen)+blg.skipf+blg.cklen*chunk
            chunk_data = blg.data[order][c_inds]
            chunk_waves = blg.wavelens[order][c_inds]
            
            # sigma clip from the top
            zinds = np.where(chunk_data>0.)[0]
            inds = zinds
            ct = len(inds)
            ran = np.arange(len(chunk_data))
            dct = 1
            while dct>0:#oldlen != newlen:
                med = np.median(chunk_data[inds])
                std = np.std(chunk_data[inds])
                inds = zinds[np.where(chunk_data[zinds]<med+3*std)[0]]
                newct = len(inds)
                dct = ct-newct
                ct = newct
                

#            ipdb.set_trace()
            blg.inds['o'+str(order)+'_c'+str(chunk)] = inds

            p0 = [0.000166,np.median(blg.data[order]),0.,0.]
#            ipdb.set_trace()
            out=scipy.optimize.leastsq(chunk_err_func,p0,args=\
                                           (full_model,blg,order,chunk),\
                                           full_output=1,maxfev=1000)
        
            message=out[3]
            ier=out[4]
            print 'Fitter status:',ier,' Message: ',message        
            p1 = out[0]
            
            #        if out[1] == None:
            #            ipdb.set_trace()
            #            continue
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
            rv.append(p1[0]*2.99e8)
            rvma.append(((p1[0]+1.)**2-1.)*2.99e8/((p1[0]+1.)**2+1.))

            params_list.append(p1)
            inds = blg.inds['o'+str(order)+'_c'+str(chunk)]
            
#            ipdb.set_trace()
            """
#            plt.plot(blg.wavelens[order][c_inds][inds],#offset*order+\
#                         cfit_model(p0,full_model,blg.wavelens[order])[c_inds][inds],\
#                         'g',zorder=2)
            plt.plot(blg.wavelens[order][c_inds][inds],#offset*order+\
                         cfit_model(p1,full_model,blg.wavelens[order])[c_inds][inds],\
                         'r',zorder=2,linewidth=2)
            plt.plot(blg.wavelens[order][c_inds][inds],#offset*order+\
                         blg.data[order][c_inds][inds],'b',zorder=1)
            
            plt.errorbar(blg.wavelens[order][c_inds][inds],\
                             blg.data[order][c_inds][inds],\
                             'b',zorder=1)
            """
        """
        plt.plot(blg.wavelens[order][inds],\
                     blg.data[order][inds],'b',zorder=1)
#        plt.plot(blg.wavelens[order][inds],\
#                     cfit_model(p0,full_model,blg.wavelens[order])[inds],\
#                     'g',zorder=2)
        plt.plot(blg.wavelens[order][inds],\
                     cfit_model(p1,full_model,blg.wavelens[order])[inds],\
                     'r',zorder=2,linewidth=2)
        """
        """
        plt.plot(blg.wavelens[order],blg.data[order]+offset*order,'b',zorder=1)
        plt.plot(blg.wavelens[order],offset*order+fit_model(p0,full_model,blg.wavelens[order]),\
                     'g',zorder=2)
        plt.plot(blg.wavelens[order],offset*order+fit_model(p1,full_model,blg.wavelens[order]),\
                     'r',zorder=2)
        
        """


    barycorr = -7204.6
    ipdb.set_trace()
    rv = (np.array(rv)+barycorr)/1000.
#    rvma = (np.array(rvma)+barycorr)/1000.
#    rver = np.array(rver)/1000.
    plt.plot(np.arange(len(rv)),rv,'o',label='Individual order RV')
#    plt.plot(np.arange(len(rvma)),rvma,'o',label='Individual order RVMA')
#    trv = np.concatenate([rv[0:5],rv[7:11],rv[13:]])

    mp0 = np.array([42.,5.,.01,3.])
    max = scipy.optimize.fmin(mixmod.maxlikelihood,mp0,args=(rv,))

    all_rv = max[0]

    low = all_rv-max[1]
    low1 = all_rv-max[1]*max[3]
    high = all_rv+max[1]
    high1 = all_rv+max[1]*max[3]

    plt.plot([0,len(blg.fit_chunks)*len(blg.fit_orders)],[all_rv,all_rv],'r-',label='Mean')
    plt.axhspan(low,high,color='y',alpha=0.5,lw=0,label='one sigma')
    plt.axhspan(low1,high1,color='r',alpha=0.5,lw=0,label='one sigma*gamma')
    plt.xlabel('"chunk"')
    plt.ylabel('RV [km s$^{-1}$]')
    plt.legend()
    plt.show()


    ipdb.set_trace()
    while True:
        xfit = np.linspace(0,100,1000)
        plt.plot(xfit,mixmod.modelpdf(max,xfit),label='Total pdf')
        plt.plot(xfit,(1-max[2])*mixmod.modelpdf([max[0],max[1],0.,max[3]],xfit),\
                     label='"Good" pdf')
        plt.plot(xfit,max[2]*mixmod.modelpdf([max[0],max[1],1.,max[3]],xfit),\
                     label='"Bad" pdf')
        plt.legend()
        bins = int(raw_input('bins(90?):'))
        plt.hist(rv,bins=bins,normed=1)
        plt.show()
#        ipdb.set_trace()

    ipdb.set_trace()
    t =5

    
    
