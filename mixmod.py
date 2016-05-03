import numpy as np
import scipy
from scipy import optimize
from scipy import special
import ipdb

def maxlikelihood(p,x):
    """
    Need to look into nedler mead minimizer, found paper from scipy citation.
    """

    mu1,sigma,beta,gamma = p
    # start likelihood as identity 
    lhood = 1.
    # for each point, multiply it to the likelihood
    for d in x:
        lhood *= (1./np.sqrt(sigma**2*2.*np.pi))*\
            ((beta/gamma)*(np.exp(-(mu1-d)**2/(2*(gamma*sigma)**2))) +\
                 (1.-beta)*(np.exp(-(mu1-d)**2/(2*sigma**2))))
    return -lhood

def modelpdf(p,x):
    mu,sigma,beta,gamma = p
    return (1./np.sqrt(sigma**2*2.*np.pi))*\
        ((beta/gamma)*(np.exp(-(mu-x)**2/(2*(gamma*sigma)**2))) +\
             (1.-beta)*(np.exp(-(mu-x)**2/(2*sigma**2))))
    
def conslikelihood(p,x,err):
    mu,sigmax = p

    try:
        np.log(sigmax/err)
    except:
        ipdb.set_trace()
    return (scipy.special.erf(np.abs(mu-x))/(np.sqrt(2)*err) - \
                scipy.special.erf(np.abs(mu-x))/(np.sqrt(2)*sigmax))/\
                (2.*np.abs(mu-x)*np.log(sigmax/err))

def conslikelihood2(p,x,err):
    mu = p
    return (scipy.special.erf(np.abs(mu-x))/(np.sqrt(2)*err))/(np.abs(mu-x))

def consmaxlike(p,x,err):
    lhood = 1.
    if p[1]<err.max():
        p[1]=err.max()
    for ind in np.arange(len(x)):
        lhood *= conslikelihood(p,x[ind],err[ind])

    return -lhood

if __name__ == '__main__':
    ipdb.set_trace()
