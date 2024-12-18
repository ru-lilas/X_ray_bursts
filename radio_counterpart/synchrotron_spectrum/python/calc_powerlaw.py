from imports import *
from parameter import *
from calc import gamma_factor_min as gamma_min
from calc import M_ej,r_sh,l_sh

def calc_N0(p,gamma_min,eps_e,M_ej,v_sh,r_sh,l_sh):
    K = (0.5*M_ej*v_sh**2).to(u.erg)
    V_shell = (4*PI*r_sh**2*l_sh).to(u.cm**3)
    val = \
        (p-2)*eps_e*K*gamma_min**(p-2) \
        / (V_shell*mc2)
    return val.to(u.cm**(-3))

N0 = calc_N0(p,gamma_min,eps_e,M_ej,v_sh,r_sh,l_sh)
print(f'N0 = {N0:.2g}')