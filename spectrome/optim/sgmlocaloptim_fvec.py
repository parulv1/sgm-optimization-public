import numpy as np
from ..forward import ntf_local as nt
from scipy.stats import pearsonr
from ..utils import functions

def local_corr(x,data,fvec):

    htotal = np.empty((len(fvec[:]),1), dtype=complex)
    for i in range(len(fvec)):
        htotal[i] = nt.ntf_local(x,fvec[i]) 
    
    spectrum = np.abs(htotal)

    filtered = functions.mag2db(spectrum[:,0])

    corrs = -pearsonr(data, filtered)[0]
    
    return corrs