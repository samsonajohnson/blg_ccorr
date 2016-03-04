#Just need for spectra classes
import corrspline as corr
#For rebinning of models
import model_bin as mb
import ipdb
import numpy as np
import matplotlib.pyplot as plt

def model_chunk(wavelens,p):
    full_model = corr.highres_spec('./t05500_g+0.5_m10p04_hr.fits')
    model = mb.john_rebin(full_model.wavelens,full_model.data,wavelens,p[0])
    #S make a chunk of the model shifted by doppler factor z
    coeffs = np.polyfit(wavelens,model,5)
    cs_model = model - np.polyval(coeffs,wavelens) + 1.
    return p[1]*model/np.sum(model) + p[2]



if __name__ == '__main__':
    ipdb.set_trace()

