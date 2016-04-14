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
    return np.exp(-(xip/p[1])**2.)

def err_func(p,model,star,order):
    wavelens = star.wavelens[order]
    inds = star.inds[str(order)]
    fitspec = star.data[order][inds]
    return (cfit_model(p,model,wavelens)[inds] -\
                fitspec)/np.sqrt(np.abs(fitspec))

def chunk_err_func(p,model,star,order,chunk):
    offset = 48
    c_length = 400
    #S these are indices that are included in the chunk
    c_inds = np.arange(offset+chunk*c_length,offset+(chunk+1)*c_length)
    wavelens = star.wavelens[order][c_inds]
    #S these are indices OF THE CHUNK that have been sigma clipped
    inds = star.inds[str(order)+'_'+str(chunk)]
    fitspec = star.data[order][c_inds][inds]

    return (cfit_model(p,model,wavelens)[inds] -\
                fitspec)/np.sqrt(np.abs(fitspec))
    


if __name__ == '__main__':
#    full_model = corr.highres_spec('./t05500_g+0.5_m10p04_hr.fits')
    full_model = corr.highres_spec('./t05500_g+4.0_p00p00_hrplc.fits')
    suffix = '_norm2'
    blg = dill.load(open('./BLG0966'+suffix+'.pkl','rb'))
#    hr = dill.load(open('./HR4963'+suffix+'.pkl','rb'))
#    hd = dill.load(open('./HD142527'+suffix+'.pkl','rb'))
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
    params = []
    blg.inds = {}
    params_list=[]

    offind = 48
    offset = 250
    chunk_length = 400
    for order in [2,3,4,5]:
#    for order in np.arange(len(blg.data)-28)+2:
        for chunk_ind in  np.arange((len(blg.data[order])-offind)/chunk_length):
#            ipdb.set_trace()
            c_inds = np.arange(offind+chunk_ind*chunk_length,\
                                   offind+(chunk_ind+1)*chunk_length)
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
            blg.inds[str(order)+'_'+str(chunk_ind)] = inds

            p0 = [0.00014,1.,np.median(blg.data[order]),0.,0.]
#            ipdb.set_trace()
#            chunk_err_func(p0,full_model,blg,order,chunk_ind)
            out=scipy.optimize.leastsq(chunk_err_func,p0,args=\
                                           (full_model,blg,order,chunk_ind),\
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
            params_list.append(p1)
            inds = blg.inds[str(order)+'_'+str(chunk_ind)]
            
#            ipdb.set_trace()
            """
            plt.plot(blg.wavelens[order][c_inds][inds],#offset*order+\
                         cfit_model(p0,full_model,blg.wavelens[order])[c_inds][inds],\
                         'g',zorder=2)
            plt.plot(blg.wavelens[order][c_inds][inds],#offset*order+\
                         cfit_model(p1,full_model,blg.wavelens[order])[c_inds][inds],\
                         'r',zorder=2,linewidth=2)
            plt.plot(blg.wavelens[order][c_inds][inds],#offset*order+\
                         blg.data[order][c_inds][inds],'b',zorder=1)
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
    plt.xlabel('Wavelength')
    plt.ylabel('Counts')#''Arbitrary')
#    plt.axis([8400,9000,0,500])
    plt.show()

    barycorr = 24167.
    ipdb.set_trace()
    rv = (np.array(rv)+barycorr)/1000.
#    rver = np.array(rver)/1000.
    plt.plot(np.arange(len(rv)),rv,'o',label='Individual order RV')
    trv = np.concatenate([rv[0:5],rv[7:11],rv[13:]])
    all_rv = 49.791 + barycorr/1000.
    plt.plot([0,32],[all_rv,all_rv],'r-',label='RV from simul. fit')
#    plt.errorbar(np.arange(len(rv)),rv,yerr=rver)
    plt.xlabel('"Order"')
    plt.ylabel('RV [km s$^{-1}$]')
    plt.legend()
    plt.show()
    ipdb.set_trace()
