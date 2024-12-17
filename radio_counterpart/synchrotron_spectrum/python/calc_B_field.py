from imports import *
from calc import M_ej, v_sh, r_sh, l_sh, e_B

def calc_mag_B2(eps_B,M_ej,v_sh,r_sh,l_sh):
    B2 = ((eps_B*M_ej*v_sh**2) / (r_sh**2*l_sh)).to(u.erg*u.cm**(-3))
    return B2

B = np.sqrt(calc_mag_B2(e_B,M_ej,v_sh,r_sh,l_sh))
print(f'B = {B:.2g}')