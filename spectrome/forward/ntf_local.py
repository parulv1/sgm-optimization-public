import numpy as np

def ntf_local(param,freq):
    gei = param[0]
    gii = param[1]
    tau_e = param[2]/1000
    tau_i = param[3]/1000
    w = 2 * np.pi * freq
    
    gee = 1
   
    Fe = np.divide(1 / tau_e ** 2, (1j * w + 1 / tau_e) ** 2)
    Fi = np.divide(1 / tau_i ** 2, (1j * w + 1 / tau_i) ** 2)
    
    Hed = (1 + (Fe * Fi * gei)/(tau_e * (1j * w + Fi*gii/tau_i)))/(1j * w + Fe*gee/tau_e + (Fe * Fi * gei)**2/(tau_e * tau_i * (1j * w + Fi * gii / tau_i)))
    
    Hid = (1 - (Fe * Fi * gei)/(tau_i * (1j * w + Fe*gee/tau_e)))/(1j * w + Fi * gii/tau_i + (Fe * Fi * gei)**2/(tau_e * tau_i * (1j * w + Fe*gee / tau_e)))


    Htotal = Hed + Hid
    
    return Htotal