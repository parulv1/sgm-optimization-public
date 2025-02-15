import numpy as np
from sympy import *


def local_stability(parameters):

    tau_e = parameters["tau_e"]
    tau_i = parameters["tau_i"]
    # # speed = parameters["speed"]
    gei = parameters[
        "gei"
    ]  # excitatory-inhibitory synaptic conductance as ratio of E-E syn
    # gee = parameters["gee"] # excitatory-excitatory synaptic conductance
    gii = parameters[
        "gii"
    ]  # inhibitory-inhibitory synaptic conductance as ratio of E-E syn
    gee = 1

    fe = float(1/tau_e)
    fi = float(1/tau_i)
    s = Symbol('s')

    
    a = poly(
        ( (s * (s+fe)**2 * (s+fi)**2 + gee * fe**3 * (s+fi)**2) *  
        (s * (s+fe)**2 * (s+fi)**2 + gii * fi**3 * (s+fe)**2) + 
        gei**2 * fe**5 * fi**5 ) / (fe**5 * fi**5), 
        s
    )


    b = a.all_coeffs()
    roots = np.roots(b)

    st = 0

    for result in roots:
        if result.real>0:
            st += 1

    return st
