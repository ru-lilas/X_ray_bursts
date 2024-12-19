from imports import *
from parameter import zeta_e, v_sh

def calc_gamma_min(zeta,v_sh):
    val = (zeta*m_p*v_sh**2) / ( 2*mc2 )
    return val.decompose()

gamma_min = calc_gamma_min(zeta_e,v_sh)
print(f'gamma_min = {gamma_min:.2g}')