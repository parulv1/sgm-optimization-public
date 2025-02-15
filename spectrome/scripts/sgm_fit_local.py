import sys
sys.path.append("../../")

from spectrome.utils import path
import numpy as np
from spectrome.optim import sgmlocaloptim_fvec
from spectrome.optim import sgmlocaloptim_pearsonp_fvec
from spectrome.brain import Brain
from scipy.optimize import basinhopping
from spectrome.stability import localstability

import time
import xarray as xr

import itertools
import multiprocessing
from functools import partial
import pickle as pkl


# Load data
ind_psd_xr = xr.open_dataarray('../data/individual_psd_reordered_matlab.nc')
ind_psd = ind_psd_xr.values

# Load frequency vector
fvec = ind_psd_xr["frequencies"].values

print("Loaded data", flush= True)

# Initial guess for parameters gei, gii, taue, taui
allx0 = np.array(
    [
        0.4,
        1,
        17,
        17,
    ]
)

# Parameter bounds
bnds1 = ((0.001,0.7), (0.001,1.9), (5.0,30.0), (5.0,30.0))
bnds2 = ((0.001,0.5), (0.001,1.9), (5.0,30.0), (5.0,30.0))
bnds3 = ((0.001,0.4), (0.001,1.9), (5.0,30.0), (5.0,30.0))
bnds4 = ((0.001,0.3), (0.001,1.9), (5.0,30.0), (5.0,30.0))
bnds5 = ((0.001,0.25), (0.001,1.9), (5.0,30.0), (5.0,30.0))

nsubs = ind_psd.shape[2]
nroi = range(ind_psd.shape[0])

print('nsubs:',nsubs)
print('nrois:',nroi)

subid = [int(sys.argv[1])]
print('subid:', subid)

paramlist = list(itertools.product(subid,nroi))

print("paramlist len:", len(paramlist))
print('subid:', subid)

# Function for calculating how many times the parameters hit the bounds 
def flagout(optres,bnd):
    return sum(optres["x"][i] in [bnd[i][0], bnd[i][1]] for i in range(4))
    
# Optimization function
def optsgm(dat,fvec,it,bnds):
    i = it[0]
    j = it[1]
    
    F_ind_db = 10*np.log10(dat[j,:,i]) # roi x freq x sub
    
    opt_res1 = basinhopping(
    sgmlocaloptim_fvec.local_corr,
    x0=allx0,
    minimizer_kwargs={"args": (F_ind_db, fvec),"bounds":bnds},
    niter=2000,
    T=0.1,
    stepsize=4,
    seed=24,
    niter_success=200,
    disp=False,
    )
    p1 = -opt_res1["fun"]
    
    if flagout(opt_res1,bnds) >=1:
        opt_res2 = basinhopping(
        sgmlocaloptim_fvec.local_corr,
        x0=allx0,
        minimizer_kwargs={"args": (F_ind_db, fvec),"bounds":bnds},
        niter=2000,
        T=0.1,
        stepsize=6,
        seed=24,
        niter_success=200,
        disp=False,
        )
        p2 = -opt_res2["fun"]
    else:
        p2 = -100
        
    
    if p1>=p2:
        opt_res = opt_res1
    else:
        opt_res = opt_res2
    
    gei = opt_res["x"][0]
    gii = opt_res["x"][1]
    tau_e = opt_res["x"][2]
    tau_i = opt_res["x"][3]
    
    f = flagout(opt_res,bnds)
    
    pval = sgmlocaloptim_pearsonp_fvec.local_corr2(opt_res["x"], F_ind_db, fvec)
    
    opt_par = [gei, gii, tau_e, tau_i, -opt_res["fun"], i, j, f, pval]
    
    return opt_par

# Main optimization function that also checks for stability. If unstable, run optimization again with different parameter bounds.
def optsgm_st(dat,fvec,it):
    
    res = optsgm(dat,fvec,it,bnds1)

    # create spectrome brain:
    brain = Brain.Brain()
    brain.ntf_params["gei"] = res[0]
    brain.ntf_params["gii"] = res[1]
    brain.ntf_params["tau_e"] = res[2]/1000
    brain.ntf_params["tau_i"] = res[3]/1000

    
    st = localstability.local_stability(brain.ntf_params)
    
    if st>0:
        res = optsgm(dat,fvec,it,bnds2)

    brain.ntf_params["gei"] = res[0]
    brain.ntf_params["gii"] = res[1]
    brain.ntf_params["tau_e"] = res[2]/1000
    brain.ntf_params["tau_i"] = res[3]/1000

    st = localstability.local_stability(brain.ntf_params)

    if st>0:
        res = optsgm(dat,fvec,it,bnds3)

    brain.ntf_params["gei"] = res[0]
    brain.ntf_params["gii"] = res[1]
    brain.ntf_params["tau_e"] = res[2]/1000
    brain.ntf_params["tau_i"] = res[3]/1000

    st = localstability.local_stability(brain.ntf_params)

    if st>0:
        res = optsgm(dat,fvec,it,bnds4)
        
    brain.ntf_params["gei"] = res[0]
    brain.ntf_params["gii"] = res[1]
    brain.ntf_params["tau_e"] = res[2]/1000
    brain.ntf_params["tau_i"] = res[3]/1000

    st = localstability.local_stability(brain.ntf_params)

    if st>0:
        res = optsgm(dat,fvec,it,bnds5)

    return res



if __name__ == '__main__':
    start = time.time()
    #Generate processes equal to the number of cores
    pool = multiprocessing.Pool(2)

    # Distribute the parameter sets evenly across the cores
    print("Starting optimization of regional data", flush = True)
    func = partial(optsgm_st, ind_psd, fvec)
    res  = pool.map(func,paramlist)
    pool.close()
    res2 = np.array(res)
    filename = f'../results/filename_{subid}.p'
    with open(filename, 'wb') as f:
        pkl.dump(res2, f)

    print("Finished data optimization", flush=True)

    end = time.time()

    print(f"Runtime of the optimization is {(end - start)/3600} hours",flush=True)