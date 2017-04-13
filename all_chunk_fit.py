#butler 1996 sec 3

#Just need for spectra classes
import corrspline as corr
#For rebinning of models
import model_bin as mb
import ipdb
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
import ghip
import lmfit


def cfit_model(p,model_ck,wavelens,order):
    factor = np.polyval(p[:1:-1],full_model.wavelens)
    unbinmodel = factor*ipmodel/np.max(ipmodel)
    model = mb.john_rebin(full_model.wavelens,full_model.data,wavelens,p[0])
    return model


def gaussian(amp,ipwidth):
    xip = np.arange(-50,50+.25,.25)
    return amp*np.exp(-(xip)**2./ipwidth)


def all_err_func(p,model,star,plot_check=False):
    z = p['z'].value
    amp = p['amp'].value
    ipwidth = p['ipwidth'].value
    all_errs = np.array([])
    all_edge_errs = np.array([])
    j=0
    skipf = 104#148
    cklen = 368#300
#    print p[0], p[1]

    #S convolve the entire model spectrum, want to trim to just relevant wave-
    #S lengths, as well as only when the ip changes. 
#    ipdb.set_trace()
    full_model.ipdata = mb.numconv(full_model.data,gaussian(amp,ipwidth))
    for order in star.fit_orders: 
        for ck in star.fit_chunks:
            #S indices for all the wavelengths in a chunk, for rebinning 
            #S purposes
            ckinds = np.arange(cklen) + skipf + ck*cklen
            #S wavelengths we want to rebin for in the chunk
            ckwaves = star.wavelens[order][ckinds]
            ckdata = star.data[order][ckinds]
            ckerrs = star.errs[order][ckinds]

            # get the wavelengths in the model relevant to the order
            lowwave = ckwaves.min()-10.
            highwave = ckwaves.max()+10.
            low = np.where(full_model.wavelens > lowwave)[0]
            high = np.where(full_model.wavelens[low] < highwave)[0]
            # the relevant model inds for the order
            minds = low[high]

            #S get the model wavelens and ipdata for the order
            modwaves = full_model.wavelens[minds]
            modck = full_model.ipdata[minds]
            #S the index 'zero' point
            c0 = p['c0_o'+str(order)+'_c'+str(ck)].value
            c1 = p['c1_o'+str(order)+'_c'+str(ck)].value
            c2 = p['c2_o'+str(order)+'_c'+str(ck)].value
#            j = (order-np.min(star.fit_orders))*5+ck
#            temp_params = [p[0],p[1],p[3*j+2],p[3*j+3],p[3*j+4]]

            factor = np.polyval([c2,c1,c0],modwaves)
            unbinmodel = factor*modck/np.max(modck)
            tmodel = mb.john_rebin(modwaves,unbinmodel,ckwaves,z)
            
            inds = blg.inds['o'+str(order)+'_c'+str(ck)]
            star.fit_data['o'+str(order)+'_c'+str(ck)] = tmodel[inds]
            
            if plot_check:
                #ipdb.set_trace()
                plt.plot(ckwaves[inds],ckdata[inds],zorder=2)
                plt.plot(ckwaves[inds],tmodel[inds],linewidth=2,zorder=2)

            errs = (tmodel[inds] - ckdata[inds])/(ckerrs[inds])#/1.75)

            if ck != star.fit_chunks[0]:
                edge_errs = (star.fit_data['o'+str(order)+'_c'+str(ck-1)][-1]-\
                     star.fit_data['o'+str(order)+'_c'+str(ck)][0])

            else: 
                edge_errs = 0.
#            all_errs=np.concatenate([all_errs,errs])
            all_errs=np.concatenate([all_errs,errs,[edge_errs]])
            all_edge_errs=np.concatenate([all_edge_errs,[edge_errs]])

        if plot_check:
            print('All errors, and the shape')
            print all_errs, np.shape(all_errs)
            print('All edge errors, and the shape')
            print all_edge_errs, np.shape(all_edge_errs)
            plt.plot(blg.wavelens[order],blg.data[order],zorder=1)
#            plt.errorbar(blg.wavelens[order],blg.data[order],blg.errs[order],
#                         zorder=1)
#            imsave_path = '/Users/samsonjohnson/Desktop/spitzer/blg_ccorr'+\
#                '/images0412/'+'orders'+str(min(blg.fit_orders))+'_'+\
#                str(max(blg.fit_orders))+'/order'+str(order)+'cklen'+\
#                str(cklen)+'.pdf'
#            plt.savefig(imsave_path)
            plt.show()
    
    if plot_check:
        prevpts = 0
        prevmax = 0
        for order in blg.fit_orders:
            for ck in blg.fit_chunks:
                x = np.arange(len(blg.inds['o'+str(order)+'_c'+str(ck)]))+\
                    prevpts
                prevpts += len(x)
                tdata = blg.fit_data['o'+str(order)+'_c'+str(ck)]
                adata = blg.data[order][np.arange(400)+skipf+ck*cklen]\
                    [blg.inds['o'+str(order)+'_c'+str(ck)]]
                errs = blg.errs[order][np.arange(400)+skipf+ck*cklen]\
                    [blg.inds['o'+str(order)+'_c'+str(ck)]]/1.75
                plt.plot(x,prevmax+np.cumsum(((tdata-adata)/errs)**2),\
                             label=str(order)+' '+str(ck))
#                ipdb.set_trace()
                prevmax = np.cumsum(((tdata-adata)/errs)**2)[-1] + prevmax
        plt.plot(np.arange(prevpts))
        plt.axis([0,prevpts+1000,0,prevpts+1000])
        plt.legend()
        plt.show()
    return np.array(all_errs)

    

def err_func(p,model,wavelens,spec):
    #S i think the fitter minimizes the square of this, no need to square here?
#    conv_model = mb.numconv(fit_model(p,model,wavelens),gaussian(p))
    return (fit_model(p,model,wavelens) - spec)
    


if __name__ == '__main__':
    # get the full highres model spectrum class
#    full_model = corr.highres_spec('./t05500_g+4.5_p00p00_hrplc.fits')
#    full_model = corr.highres_spec('./t06500_g+4.0_p00p00_hrplc.fits')
    # the blg spectrum class
    blg = corr.multi_spec('./blg0966red_multi.fits')
    full_model = corr.phe_spec('./lte05500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits','./WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',minwave=min(blg.wavelens[-1])-500.,maxwave=max(blg.wavelens[0])+500.)
#    full_model = corr.phe_spec('./lte06300-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits','./WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',minwave=min(blg.wavelens[-1])-500.,maxwave=max(blg.wavelens[0])+500.)
    # the simultaneous fit parameters, z and ip width (just a gaussian)
#    p0 = [0.0001,1.]
#    ipdb.set_trace()
    p0 = lmfit.Parameters()
    p0.add('z',value=0.00016)
    p0.add('amp',value=-1.)
    p0.add('ipwidth',value=20.,min=1.,max=50.)
    # the empty dictionary for putting the trimmed indices
    blg.inds={}
    # the pixels to sjip at the beginnning of each order
    skipf = 104#148
    # the length of chunk we are using
    cklen = 368#300
    # empty dict for data
    blg.fit_data = {}
    # list of orders we want to fit
    blg.fit_orders = np.arange(8)+2#[2,3]#,4,5,6]#,7,8,9,10,11,12]#np.arange(len(blg.data)-19)+2
    # number of chunks we will be fitting based on cklen and skipf
    blg.fit_chunks = np.arange((len(blg.data[0])-skipf)/cklen)
    print(blg.fit_chunks)

    # go through and clip 'bad' points. here, bad is points that are three 
    # sigma above the median counts in a chunk, and points lower than zero
    for order in blg.fit_orders:
        # temp for order data
        odata = blg.data[order]
        # a check for plotting
#        plt.plot(blg.wavelens[order],blg.data[order],zorder=1)

        for ck in blg.fit_chunks:
            # temp for data in a chunk
            ckdata = odata[ck*cklen+skipf:(ck+1)*cklen+skipf]
            # getting initial guess for fit parameters, assmuming a flat 
            # wavelength dependence
            p0.add('c0_o'+str(order)+'_c'+str(ck),value=np.median(ckdata))
            p0.add('c1_o'+str(order)+'_c'+str(ck),value=0.0)
            p0.add('c2_o'+str(order)+'_c'+str(ck),value=0.0)
#            p0.append(np.median(ckdata))
#            p0.append(0.)
#            p0.append(0.)

            # indices greater than zero
            zinds = np.where(ckdata>0.)[0]
            # set inds to zinds for intial count of bad inds
            inds = zinds
            # get the count of inds, as w are watching for when the change 
            # goes to zero
            ct = len(inds)
            
            # zero index indices for the chunk
            ran = np.arange(len(ckdata))
            # the change in inds count, initialize at 1 so we enter the loop
            dct = 1
            
            # while there is a change in the number of inds being counted
            while dct>0:#oldlen != newlen:
                med = np.median(ckdata[inds])
                std = np.std(ckdata[inds])
                inds = zinds[np.where(ckdata[zinds]<med+3*std)[0]]
                newct = len(inds)
                dct = ct-newct
                ct = newct
            # make the list of inds a dict entry for blg for later reference
            blg.inds['o'+str(order)+'_c'+str(ck)] = inds# + skipf + cklen*ck
            # initialize the length of the fit_data as zeros of len(inds)
            blg.fit_data['o'+str(order)+'_c'+str(ck)] = np.zeros(len(inds))



    
    """
    for debugging
    tma = mb.john_rebin(full_model.wavelens,full_model.data,blg.wavelens[3],0.00017)
    plt.plot(tma/np.max(tma),'.')
    plt.plot(blg.data[3]/np.max(blg.data[3]),'.')
    for ind in [4]:
        plt.plot(blg.inds['o3_c'+str(ind)],tma[blg.inds['o3_c'+str(ind)]]/np.max(tma),zorder=2)
        plt.plot(blg.inds['o3_c'+str(ind)],blg.data[3][blg.inds['o3_c'+str(ind)]]/np.max(blg.data[3]),zorder=2)

    """    

#    ipdb.set_trace()
    all_err_func(p0,full_model,blg)
#    print('Skipf: '+str(skipf)+', cklen: '+str(cklen))
#    print('Orders being fit:')
#    print(blg.fit_orders)
    result = lmfit.minimize(all_err_func,p0,args=(full_model,blg),method='leastsq')
    lmfit.report_fit(result,show_correl=False)
#    ipdb.set_trace()
    all_err_func(result.params,full_model,blg,plot_check=True)

#    ipdb.set_trace()
    bins = np.arange(len(blg.data[0]))
    offset = 250
    j=0

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
