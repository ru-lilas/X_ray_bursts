import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy import units as u
from astropy import constants as const
import numpy as np

# physical and mathematical constants
PI = np.pi
c = const.c.cgs
G = const.G.cgs
M = const.M_sun.cgs
e = const.e.gauss
m_e = const.m_e.cgs
m_p = const.m_p.cgs
sigma_T = const.sigma_T.cgs

v_sh = 0.3*c
t_rec = 200 * u.min
t_rise = 3 *u.min
tau = 10 * u.min
e_B = 0.01
def calc_L_Edd(mass):
    value = 4*PI*G*mass*m_p*c / sigma_T
    return value.to(u.erg / u.year)

def calc_Mdot_Edd(mass):
    L_Edd = calc_L_Edd(mass)
    value = L_Edd / c**2
    return value

def calc_mag_B(eps_B,ej_mass,v_sh,r_sh,t_decay):
    B2 = (eps_B*ej_mass*v_sh / (r_sh**2*t_decay)).to(u.erg/(u.cm)**3)
    print(f'B^2 = {B2:3g}')
    return np.sqrt(B2)

M_Edd_sol = calc_Mdot_Edd(M)
print(f'M_Edd_sol = {M_Edd_sol.to(u.g/u.s):3g} = {(M_Edd_sol/M).value:3g} M_sol {u.year**(-1)}')

eta_acc = 0.1
Mdot_4U1728 = eta_acc * calc_Mdot_Edd(1.4*M)
print(f'Typical Mdot of 4U 1728-34 is {Mdot_4U1728.to(u.g/u.s):3g}')

# mass of shock ejecta
M_ej = Mdot_4U1728 * t_rec
print(f'Ejecta mass = eta * {(M_ej).to(u.g):3g}')

# strength of magnetic field
r_sh = v_sh * t_rise
B = calc_mag_B(e_B,M_ej,v_sh,r_sh,tau)
print(B.value*u.G)

# li_const = [['c',c],['M_sol',M],['m_e',m_e],['m_p',m_p]]
# df_const = pd.DataFrame(li_const, columns=['Name', 'Value'])
# print(df_const)
# energy_density = u.def_unit('energy_density',u.erg / (u.cm)**3)
# print(f'm_e/m_p = {m_p/m_e :3g}')
# M_x = 1.0 * M
# e_Edd = 0.1
# e_B = 0.01
# zeta_e = 0.4 # 無衝突加熱における, 熱的陽子の, 熱的電子へのエネルギー供給率

# L_Edd = 10**38 * M_x  * u.erg / (M*u.s)
# Md_Edd = (L_Edd*0.1*M_x / (e_Edd*M*c**2)).to(u.g/u.s)

# print(f'v_sh = {v_sh:3g}')
# print(f't_rise = {t_rise:3g}')
# print(f'L_Edd = {L_Edd:3g}')
# print(f'Mdot_Edd = {Md_Edd:3g}')

# M_ej = Md_Edd * t_rec
# B2 = e_B * Md_Edd * t_rec / (v_sh * t_rise**2 * tau)
# B2 = M_ej * v_sh / (r_sh**2 * tau)
# B = np.sqrt(B2.cgs)
# print(f'{(r_sh/const.au).cgs:3g}\n{M_ej:3g}\n{B2:3g}')
# print(f'{B:3g}')
# print(v_sh * t_rise**2 * tau)
# # print(f'{(e*B).cgs}')
# print(Md_Edd * t_rec)
# print(np.sqrt(7636))