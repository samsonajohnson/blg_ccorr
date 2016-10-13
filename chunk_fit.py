#Just need for spectra classes
import corrspline as corr
#For rebinning of models
import model_bin as mb
import ipdb
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import dill
import scipy
from scipy import optimize
import mixmod
import lmfit

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

def gaussian(ipwidth):
    xip = np.arange(-50,50+.25,.25)
    return np.exp((-(xip)**2.)/ipwidth)

def err_func(p,model,star,order):
    wavelens = star.wavelens[order]
    inds = star.inds[str(order)]
    fitspec = star.data[order][inds]
    return (cfit_model(p,model,wavelens)[inds] -\
                fitspec)/np.sqrt(np.abs(fitspec))

def chunk_func(p,model,star,order,chunk):
    z = p['z'].value
    ipwidth = p['ipwidth'].value
    c0 = p['c0'].value
    c1 = p['c1'].value
    c2 = p['c2'].value

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
    modck = mb.numconv(full_model.data[minds],gaussian(ipwidth))
    #S add some wavelength scaling
    factor = np.polyval([c2,c1,c0],modwaves)
    unbinmodel = factor*modck/np.max(modck)
    #S rebin the convolved model, and 
    tmodel = mb.john_rebin(modwaves,unbinmodel,ckwaves,z)

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
    errs = (tmodel[inds] - ckdata[inds])/(ckerrs[inds])#/4.75)
    return errs

if __name__ == '__main__':
#    full_model = corr.highres_spec('./t05500_g+4.0_p00p00_hrplc.fits')
#    full_model = corr.highres_spec('./t06500_g+4.0_p00p00_hrplc.fits')
    blg = corr.multi_spec('./blg0966red_multi.fits')
#    full_model = corr.phe_spec('./lte05500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits','./WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',minwave=min(blg.wavelens[-1])-500.,maxwave=max(blg.wavelens[0])+500.)
    full_model = corr.phe_spec('./lte06300-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits','./WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',minwave=min(blg.wavelens[-1])-500.,maxwave=max(blg.wavelens[0])+500.)
    ipdb.set_trace()
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

    blg.skipf = 48#100
    blg.cklen = 400#278
    numcks = 8
    if 'lte' in full_model.name:
        blg.fit_orders = np.arange(numcks)
    else:
        blg.fit_orders = np.arange(numcks)+2
#    blg.fit_orders = [3]
    blg.fit_chunks = np.arange((len(blg.data[0])-blg.skipf)/blg.cklen)
    blg.fit_data = {}
    table_lines = []
    for order in blg.fit_orders:
 #       ipdb.set_trace()
        for chunk in blg.fit_chunks:
#            ipdb.set_trace()
            c_inds = np.arange(blg.cklen)+blg.skipf+blg.cklen*chunk
            chunk_data = blg.data[order][c_inds]
            chunk_waves = blg.wavelens[order][c_inds]
#            print (chunk_waves[-1]-chunk_waves[0])/chunk_waves[-1]
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
                


            blg.inds['o'+str(order)+'_c'+str(chunk)] = inds
            p0 = lmfit.Parameters()
            p0.add('z',value=0.000163)

            if full_model.name == 't05500':
                p0.add('ipwidth',value=1.60,min=1,max=50,vary=False)
            if full_model.name == 't06500':
                p0.add('ipwidth',value=1.01,min=1,max=50,vary=False)
            if full_model.name == 'lte05500':
                p0.add('ipwidth',value=32.17,vary=False)
            if full_model.name == 'lte06300':
                p0.add('ipwidth',value=17.56,vary=False)
            p0.add('c0',value=np.median(blg.data[order]))
            p0.add('c1',value=0.0)
            p0.add('c2',value=0.0)

            result = lmfit.minimize(chunk_err_func,p0,args=(full_model,blg,order,chunk),method='leastsq')
            print lmfit.fit_report(result,show_correl=False)
            
            rv.append(result.params['z'].value*2.99e8)

            rver.append(result.params['z'].stderr*2.99e8)
#            rvma.append(((p1[0]+1.)**2-1.)*2.99e8/((p1[0]+1.)**2+1.))

#            params_list.append(p1)
            inds = blg.inds['o'+str(order)+'_c'+str(chunk)]



    barycorr = -7204.6
    #    ipdb.set_trace()
    rv = (np.array(rv)+barycorr)/1000.
    #    rvma = (np.array(rvma)+barycorr)/1000.
    rver = np.array(rver)/1000.
    rverrms = np.sqrt(np.sum(rver**2)/len(rver))
    num_plots = len(blg.fit_orders)
    colormap = plt.cm.nipy_spectral

    #   ipdb.set_trace()
    #    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])


    for order in blg.fit_orders:
        for chunk in blg.fit_chunks:
            rec_rv = rv[order+chunk]
            rec_rv_err = rver[order+chunk]
            c_inds = np.arange(blg.cklen)+blg.skipf+blg.cklen*chunk
            chunk_waves = blg.wavelens[order][c_inds]

            table_str = '%i & %i & %.2f & %.2f & %0.1f & %.1f\\ \n'%(order,chunk,np.min(chunk_waves),np.max(chunk_waves),rec_rv,rec_rv_err)
            table_lines.append(table_str)

    ipdb.set_trace()
    with open('tablefile.tex','w') as tfile:
        tfile.write('\begin{tabular}\n')
        tfile.write('\hline \n')
        tfile.write('Order&Chunk&$\lambda_{min}$&$\lambda_{max}$&$RV$&$\sigma_{RV}$\\ \n')
        tfile.write('\hline \n')
        for line in table_lines:
            tfile.write(line)
        tfile.write('\end{tabular}\n')
    rvfig, rvax1 = plt.subplots()
    for or_num in blg.fit_orders:
        lowind  = or_num*len(blg.fit_chunks)
        upind = lowind + len(blg.fit_chunks)
        print lowind, upind
        rvax1.plot(np.arange(lowind,upind),rv[lowind:upind],'s',label='Chunk from order ' + str(or_num),color=colormap(np.linspace(.9,.05,num_plots)[or_num]),ms=10,zorder=2)
        rvax1.errorbar(np.arange(len(rv)),rv,yerr=rver,linestyle='None')
#    plt.plot(np.arange(len(rvma)),rvma,'o',label='Individual order RVMA')
#    trv = np.concatenate([rv[0:5],rv[7:11],rv[13:]])

    mp1 = np.array([42.,5.,.01,20.])

    max1 = scipy.optimize.fmin(mixmod.maxlikelihood,mp1,args=(rv,))
    print max1
    mu1 = max1[0]
    sigma = max1[1]
    beta = max1[2]
    gamma = max1[3]
    mp2 = [42.,20.]

#    ipdb.set_trace()
    max2 = scipy.optimize.fmin(mixmod.consmaxlike,mp2,args=(rv,rver),\
                                   maxiter=10000)
    all_rv = max1[0]

#    ipdb.set_trace()
    low = mu1 - sigma
    low1 = mu1 - sigma*gamma
    high = mu1 + sigma
    high1 = mu1 + sigma*gamma
    

    rvax1.plot([0,len(blg.fit_chunks)*len(blg.fit_orders)],[all_rv,all_rv],'r-',label='Mean',zorder=1,linewidth=2)

    rvax1.axhspan(low,high,color='k',alpha=0.25,lw=0,label='1$\sigma$')
    rvax1.axhspan(low1,low,color='k',alpha=0.10,lw=0,label='1$\sigma*\gamma$')
    rvax1.axhspan(high,high1,color='k',alpha=0.10,lw=0)
    rvax1.set_xlabel('"Chunk"',fontsize=22)
    rvax1.set_ylabel('RV [km s$^{-1}$]',fontsize=22)
    rvax1.tick_params(axis='both', which='major', labelsize=16)
 #   plt.title(full_model.name+', '+str(max1))
#    plt.legend()

    new_labels = []
    for order in blg.fit_orders:
        label_str = '%.2f'%np.median(blg.wavelens[order])
        new_labels.append(label_str)
#    ipdb.set_trace()
    rvax2 = rvax1.twiny()
    new_ticks = np.arange(len(blg.fit_orders))*len(blg.fit_chunks)+2
    rvax2.set_xlim(rvax1.get_xlim())
    rvax2.set_xticks(new_ticks)
    rvax2.set_xticklabels(new_labels)
    rvax2.tick_params(axis='both', which='major', labelsize=16)
    rvax2.set_xlabel(r"Central wavelength of order [$\AA$]",fontsize=22)
#    matplotlib.rcParams.update({'font.size': 22})
    plt.show()


    ipdb.set_trace()
    while True:
        xfit = np.linspace(0,100,1000)
        plt.plot(xfit,mixmod.modelpdf(max1,xfit),label='Total pdf')
#        plt.plot(xfit,mixmod.
        plt.plot(xfit,(1-max1[2])*mixmod.modelpdf([max1[0],max1[1],0.,\
                                                      max1[3]],xfit),\
                     label='"Good" pdf')
        plt.plot(xfit,max1[2]*mixmod.modelpdf([max1[0],max1[1],1.,\
                                                  max1[3]],xfit),\
                     label='"Bad" pdf')
        plt.legend()
        bins = int(raw_input('bins(90?):'))
        plt.hist(rv,bins=bins,normed=1)
        plt.show()
#        ipdb.set_trace()

    ipdb.set_trace()
    t =5

    
    
