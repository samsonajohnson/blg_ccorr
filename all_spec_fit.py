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
#    coeffs = np.polyfit(wavelens,model,5)
    #S do we need to add 1? will just come out in the wash. might want to boost
    #S the observed though
#    cs_model = model - np.polyval(coeffs,wavelens) + 1.
    factor = np.polyval(p[:0:-1],wavelens)
    return factor*model/np.max(model)

def cfit_model(p,full_model,wavelens):
    #S make a chunk of the model shifted by doppler factor z
    model = mb.john_rebin(full_model.wavelens,full_model.data,wavelens,p[0])
    #S continuum subtraction
#    coeffs = np.polyfit(wavelens,model,5)
    #S do we need to add 1? will just come out in the wash. might want to boost
    #S the observed though
#    cs_model = model - np.polyval(coeffs,wavelens) + 1.
    factor = np.polyval(p[:1:-1],wavelens)
    convmodel = mb.numconv(model,gaussian(p))
    return factor*convmodel/np.max(convmodel)


def gaussian(p):
    xip = np.arange(-10,10+.25,.25)
    return np.exp(-(xip/p[1])**2.)


def all_err_func(p,model,star):
    all_errs = np.array([])
    j=0
    print p
    for order in np.arange(len(blg.data)-19)+2:
        """ old sigclip
        zind = np.where(star.data[order]>0.)[0]
        med=np.median(star.data[order][zind])
        std=np.std(star.data[order][zind])
        siginds = np.where(np.abs(star.data[order][zind]-med)<3*std)[0]
        """
        """ new sigclip
        zind = np.where(star.data[order]>0.)[0]
        oldlen = len(zind)
        newlen = -1
        tempinds = zind
        while oldlen != newlen:
            med=np.median(star.data[order][tempinds])
            std=np.std(star.data[order][tempinds])
            siginds = np.where(np.abs(star.data[order][tempinds]-med)<3*std)[0]
            tempinds = zind[siginds]
            newlen = len(tempinds)
        inds = zind[tempinds]
        """
#        temp_params = [p[0],p[1],p[2*j+2],p[2*j+3]]
        temp_params = [p[0],p[1],p[3*j+2],p[3*j+3],p[3*j+4]]
        j+=1
        inds = blg.inds[str(order)]
        all_errs=np.concatenate([all_errs,
                                 (cfit_model(temp_params,model,\
                                                star.wavelens[order])[inds]-\
                                    star.data[order][inds])])
    #S i think the fitter minimizes the square of this, no need to square here?

#    conv_model = mb.numconv(fit_model(p,model,wavelens),gaussian(p))
    return np.array(all_errs)

def err_func(p,model,wavelens,spec):
    #S i think the fitter minimizes the square of this, no need to square here?
#    conv_model = mb.numconv(fit_model(p,model,wavelens),gaussian(p))
    return (fit_model(p,model,wavelens) - spec)
    


if __name__ == '__main__':
#    full_model = corr.highres_spec('./t05500_g+0.5_m10p04_hr.fits')
    full_model = corr.highres_spec('./t05500_g+4.0_p00p00_hrplc.fits')
    suffix = '_norm2'
    blg = dill.load(open('./BLG0966'+suffix+'.pkl','rb'))
#    hr = dill.load(open('./HR4963'+suffix+'.pkl','rb'))
#    hd = dill.load(open('./HD142527'+suffix+'.pkl','rb'))

#    while True:
#        order = int(raw_input('Order to fit: '))
#    ipdb.set_trace()
    p0 = [0.00014,1.]
        #    p0 = [0.0001666,np.median(blg.data[order]),0.]
        #    p0 = [0.0001666,250000.,46.]
        #    p0 = [0.0001666,2500000./2,-500.]#4
    blg.inds={}
    for order in np.arange(len(blg.data)-19)+2:
        p0.append(np.median(blg.data[order]))
        p0.append(0.)
        p0.append(0.)
        """
        zind = np.where(blg.data[order]>0.)[0]
        med=np.median(blg.data[order][zind])
        std=np.std(blg.data[order][zind])
        siginds = np.where(np.abs(blg.data[order][zind]-med)<3*std)[0]
        blg.inds[str(order)] = zind[siginds]
        """
        zinds = np.where(blg.data[order]>0.)[0]
        inds = zinds
        ct = len(inds)

        ran = np.arange(len(blg.data[order]))
        dct = 1

        while dct>0:#oldlen != newlen:
            med = np.median(blg.data[order][inds])
            std = np.std(blg.data[order][inds])
            inds = zinds[np.where(blg.data[order][zinds]<med+3*std)[0]]
            newct = len(inds)
            dct = ct-newct
            ct = newct
#            plt.plot(np.arange(len(blg.data[order])),blg.data[order])
#            plt.plot(inds,blg.data[order][inds])
#            plt.plot([0,2048],[med+3*std,med+3*std])
#            plt.show()

        blg.inds[str(order)] = inds

 
    all_err_func(p0,full_model,blg)
    out=scipy.optimize.leastsq(all_err_func,p0,args=\
                                   (full_model,blg)\
                                   ,full_output=1,maxfev=1000)
        
    p1 = out[0]
    jacob=out[1]
    mydict=out[2]
    message=out[3]
    ier=out[4]
    resids=mydict['fvec']
    covar = np.std(resids)**2*jacob
    chisq = np.sum(resids**2)
    degs_frdm=len(blg.wavelens[order])-len(p1)
    red_chisq = chisq/degs_frdm
    print 'Fitter status:',ier,' Message: ',message
    i=0
    for u in p1:
        print 'Param: ', i+1, ': ',u,' +/-',np.sqrt(covar[i,i])
        i+=1
    print 'Chisq: ',chisq,' Reduced Chisq: ',red_chisq

    ipdb.set_trace()
    bins = np.arange(len(blg.data[0]))
    offset = 250
    j=0
 #   try:
#    for order in np.arange((len(p1)-2)/2)+2:
    for order in np.arange((len(p1)-2)/3)+2:
        """
        zind = np.where(star.data[order]>0.)[0]
        oldlen = len(zind)
        newlen = -1
        tempinds = zind
        while oldlen != newlen:
            med=np.median(star.data[order][tempinds])
            std=np.std(star.data[order][tempinds])
            siginds = np.where(np.abs(star.data[order][tempinds]-med)<3*std)[0]
            tempinds = zind[siginds]
            newlen = len(tempinds)
        inds = zind[tempinds]
        """
        inds = blg.inds[str(order)]

#        temppars0 = [p0[0],p0[1],p0[2*j+2],p0[2*j+3]]
#        temppars1 = [p1[0],p1[1],p1[2*j+2],p1[2*j+3]]
        temppars0 = [p0[0],p0[1],p0[3*j+2],p0[3*j+3],p0[3*j+4]]
        temppars1 = [p1[0],p1[1],p1[3*j+2],p1[3*j+3],p1[3*j+4]]
        j+=1
        plt.plot(bins[inds],offset*order+blg.data[order][inds],'b',zorder=1)
        plt.plot(bins[inds],offset*order+\
                     cfit_model(temppars0,full_model,\
                                    blg.wavelens[order])[inds],\
                     'g',zorder=2)
        plt.plot(bins[inds],offset*order+\
                     cfit_model(temppars1,full_model,\
                                    blg.wavelens[order])[inds],\
                     'r',zorder=2)
#    except IndexError as e:
#        print e
    ipdb.set_trace()
    plt.show()
    print str(p1[0]*3e8)+'+/-'+str(np.sqrt(covar[0,0])*3e8)
    ipdb.set_trace()
