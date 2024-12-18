from imports import *
from parameter import eps_B
from calc_shock import M_ej, v_sh, r_sh, l_sh

def calc_mag_B2(eps_B,M_ej,v_sh,r_sh,l_sh):
    B2 = ((eps_B*M_ej*v_sh**2) / (r_sh**2*l_sh)).to(u.erg*u.cm**(-3))
    return B2

B2 = calc_mag_B2(eps_B,M_ej,v_sh,r_sh,l_sh)
B = np.sqrt(B2)

print(f'B = {B:.2g}')