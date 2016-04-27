import numpy as np
import scipy
from scipy import optimize

def maxlikelihood(p,x):
    """
    Need to look into nedler mead minimizer, found paper from scipy citation.
    """

    mu,sigma,beta,gamma = p
    # start likelihood as identity 
    lhood = 1.
    # for each point, multiply it to the likelihood
    for d in x:
        lhood *= (1./np.sqrt(sigma**2*2.*np.pi))*\
            ((beta/gamma)*(np.exp(-(mu-d)**2/(2*(gamma*sigma)**2))) +\
                 (1.-beta)*(np.exp(-(mu-d)**2/(2*sigma**2))))
    return -lhood

def modelpdf(p,x):
    mu,sigma,beta,gamma = p
    return (1./np.sqrt(sigma**2*2.*np.pi))*\
        ((beta/gamma)*(np.exp(-(mu-x)**2/(2*(gamma*sigma)**2))) +\
             (1.-beta)*(np.exp(-(mu-x)**2/(2*sigma**2))))
    
