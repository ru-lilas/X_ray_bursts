from imports import *
from parameter import v_sh,t_x,t_rec,eta,t_rise
from calc_accretion_rate import Mdot

l_sh = (v_sh*t_x).to(u.cm)
M_ej = (eta * Mdot * t_rec).to(u.g)
r_sh = (v_sh*t_rise).to(u.cm)
print(
        f'l_sh = {l_sh:.2g}\n' \
        f'M_ej = {M_ej:.2g}\n' \
        f'r_sh = {r_sh:.2g}'
        )