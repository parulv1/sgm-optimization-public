import sys
sys.path.append("../../")

from spectrome.utils import functions, path
import numpy as np
import xarray as xr

from spectrome.optim import sgmglobaloptim, sgmglobaloptim_pearson
from scipy.optimize import dual_annealing
from spectrome.brain import Brain
from spectrome.stability import localstability


import time

import itertools
import multiprocessing
from functools import partial

start = time.time()


# define data directory
data_dir = path.get_data_path()

# cortical areas with MEG collected + source reconstructed
rois_with_MEG = np.arange(0,68)


## Load MEG:
## individual connectomes, this is a Nregion x Nregion x Nsubjects array:
ind_conn_xr = xr.open_dataarray('../data/individual_connectomes_reordered.nc')
ind_conn = ind_conn_xr.values

# ind_psd_xr = xr.open_dataarray(data_dir + '/individual_psd_reordered_smooth.nc')
ind_psd_xr = xr.open_dataarray('../data/individual_psd_reordered_matlab.nc')
ind_psd = ind_psd_xr.values

fvec = ind_psd_xr["frequencies"].values

nsubs = ind_psd.shape[2]
paramlist = range(nsubs)

# Parameter bounds for optimization
v_lower = 5
v_upper = 20
bnds1 = ((5.0,30.0), (5.0,200.0), (0.1,1.0), (v_lower,v_upper), (0.001,0.7), (0.001,2.0), (5.0,30.0))
bnds2 = ((5.0,30.0), (5.0,200.0), (0.1,1.0), (v_lower,v_upper), (0.001,0.5), (0.001,1.5), (5.0,30.0))
bnds3 = ((5.0,30.0), (5.0,200.0), (0.1,1.0), (v_lower,v_upper), (0.001,0.4), (0.001,1.5), (5.0,30.0))
bnds4 = ((5.0,30.0), (5.0,200.0), (0.1,1.0), (v_lower,v_upper), (0.001,0.3), (0.001,1.5), (5.0,30.0))


# Initial guesses taue, taui, alpha, speed, gei, gii, tauC


def inguess(i):
    if i==0:
        allx0 = np.array([15.0, 10.0, 1.0, 5, 0.3, 0.6, 6.0])
    if i==1:
        allx0 = np.array([25.0, 80.0, 0.5, 10, 0.2, 0.1, 15.0])
    if i==2:
        allx0 = np.array([6.0, 150.0, 0.1, 18, 0.1, 1.2, 25.0])
    return allx0

# Function for calculating how many times the parameters hit the bounds 
def flagout(optres,bnd):
    return sum(optres["x"][i] in [bnd[i][0], bnd[i][1]] for i in range(7))
    
# Optimization function
def optsgm(cdk,psd,rois_with_MEG,fvec,s,bnds):
    print('optsgm', os.getpid(), s)

    C_ind = cdk[:,:,s] # grab current subject's individual connectome
    F_ind = psd[:,:,s] # grab current subject's MEG

    F_ind_db = 10*np.log10(F_ind)

    data_dir = path.get_data_path()
    # create spectrome brain:
    brain = Brain.Brain()
    brain.add_connectome(data_dir) # grabs distance matrix
    # re-ordering for DK atlas and normalizing the connectomes:
    brain.reorder_connectome(brain.connectome, brain.distance_matrix)
    brain.connectome = C_ind # re-assign connectome to individual connectome
    brain.bi_symmetric_c()
    brain.reduce_extreme_dir()

    opt_res0 = dual_annealing(
        sgmglobaloptim.global_corr,
        x0=inguess(0),
        bounds=bnds,
        args=(brain,F_ind_db,F_ind,rois_with_MEG,fvec),
        maxiter=500,
        initial_temp=5230.0,
        seed=24,
        visit=2.62,
        no_local_search=False
    )

    # print('sub:',s,'opt_res0:',opt_res0["message"])

    opt_res1 = dual_annealing(
        sgmglobaloptim.global_corr,
        x0=inguess(1),
        bounds=bnds,
        args=(brain,F_ind_db,F_ind,rois_with_MEG,fvec),
        maxiter=500,
        initial_temp=5230.0,
        seed=24,
        visit=2.62,
        no_local_search=False
    )

    # print('sub:',s,'opt_res1:',opt_res1["message"])

    opt_res = opt_res0 if opt_res0["fun"] < opt_res1["fun"] else opt_res1

    opt_res2 = dual_annealing(
        sgmglobaloptim.global_corr,
        x0=inguess(2),
        bounds=bnds,
        args=(brain,F_ind_db,F_ind,rois_with_MEG,fvec),
        maxiter=500,
        initial_temp=5230.0,
        seed=24,
        visit=2.62,
        no_local_search=False
    )


    if opt_res2["fun"] < opt_res["fun"]:
        opt_res = opt_res2


    print('printing optim result')

    print(opt_res)

    tau_e = opt_res["x"][0]
    tau_i = opt_res["x"][1]
    alpha = opt_res["x"][2]
    speed = opt_res["x"][3]
    gei = opt_res["x"][4]
    gii = opt_res["x"][5]
    tauC = opt_res["x"][6]

    f = flagout(opt_res,bnds)

    
    
    rcorr, spcorr = sgmglobaloptim_pearson.global_corr(opt_res["x"], brain, F_ind_db, F_ind, rois_with_MEG, fvec)
    

    return [
        tau_e,
        tau_i,
        alpha,
        speed,
        gei,
        gii,
        tauC,
        -opt_res["fun"],
        rcorr,
        spcorr,
        s,
        f,
        opt_res["status"],
        opt_res["success"],
    ]

# Main optimization code. It checks for stability. If unstable, run optimization again with lower bounds
def optsgm_st(cdk,psd,rois_with_MEG,fvec,s):
    
    res = optsgm(cdk,psd,rois_with_MEG,fvec,s,bnds1)
    
    data_dir = path.get_data_path()
    # create spectrome brain:
    brain = Brain.Brain()
    brain.add_connectome(data_dir) # grabs distance matrix
    # re-ordering for DK atlas and normalizing the connectomes:
    brain.reorder_connectome(brain.connectome, brain.distance_matrix)
    brain.bi_symmetric_c()
    brain.reduce_extreme_dir()
    brain.ntf_params["tau_e"] = res[0]/1000
    brain.ntf_params["tau_i"] = res[1]/1000
    brain.ntf_params["gei"] = res[4]
    brain.ntf_params["gii"] = res[5]
    
    st = localstability.local_stability(brain.ntf_params)
    
    if st>0:
        res = optsgm(cdk,psd,rois_with_MEG,fvec,s,bnds2)

    brain.ntf_params["tau_e"] = res[0]/1000
    brain.ntf_params["tau_i"] = res[1]/1000
    brain.ntf_params["gei"] = res[4]
    brain.ntf_params["gii"] = res[5]

    st = localstability.local_stability(brain.ntf_params)

    if st>0:
        res = optsgm(cdk,psd,rois_with_MEG,fvec,s,bnds3)

    brain.ntf_params["tau_e"] = res[0]/1000
    brain.ntf_params["tau_i"] = res[1]/1000
    brain.ntf_params["gei"] = res[4]
    brain.ntf_params["gii"] = res[5]

    st = localstability.local_stability(brain.ntf_params)

    if st>0:
        res = optsgm(cdk,psd,rois_with_MEG,fvec,s,bnds4)

    return res

if __name__ == '__main__':
    print('Starting main')

    #Generate processes equal to the number of cores
    pool = multiprocessing.Pool(36)

    print('doing param optimization')
    # Distribute the parameter sets evenly across the cores

    func = partial(optsgm_st, ind_conn, ind_psd, rois_with_MEG, fvec)
    res  = pool.map(func,paramlist)
    # pool.close()
    res2 = np.array(res)
    np.savetxt("../results/filename.csv", res2, delimiter=",",header="taue, taui, alpha, speed, gei, gii, tauC, r_tot, r_psd, r_sp, sub, flag, status, success")

    print("Finished Chang data optimization for SGM")

    end = time.time()
    print("T is 5230 and visit is 2.62 and maxiter is 500 and dual annealing with three initial conditions")

    print(f"Runtime of the global sgm optimization with simulated annealing for all subjects is {(end - start)/3600} hours")



