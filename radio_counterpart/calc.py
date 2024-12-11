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
m_u = const.u.cgs
sigma_T = const.sigma_T.cgs

v_sh = 0.3*c
t_rec = 200 * u.min
t_rise = 3 *u.min
t_x = 10 * u.s
tau = 10 * u.min
e_B = 0.01
def calc_L_Edd(mass):
    value = 4*PI*G*mass*m_p*c / sigma_T
    return value.to(u.erg / u.year)

def calc_Mdot_Edd(mass):
    L_Edd = calc_L_Edd(mass)
    value = L_Edd / c**2
    return value

def calc_mag_B2(eps_B,M_ej,v_sh,r_sh,l_sh):
    B2 = ((eps_B*M_ej*v_sh**2) / (r_sh**2*l_sh)).to(u.erg/(u.cm)**3)
    return B2

def calc_rho(M_ej,r_sh,v_sh,tau):
    value = M_ej / (4*PI*r_sh**2*v_sh*tau)
    return value


M_Edd_sol = calc_Mdot_Edd(M)
print(f'M_Edd_sol = {M_Edd_sol.to(u.g/u.s):3g} = {(M_Edd_sol/M).value:3g} M_sol {u.year**(-1)}')

eta_acc = 0.1
Mdot_4U1728 = eta_acc * calc_Mdot_Edd(1.4*M)
print(f'Typical Mdot of 4U 1728-34 is {Mdot_4U1728.to(u.g/u.s):3g}')

# thickness of shock ejecta
l_sh = (v_sh *t_x).to(u.cm)

# shock position from NS surface
r_sh = (v_sh * t_rise).to(u.cm)

# mass of shock ejecta
M_ej = Mdot_4U1728 * t_rec

# printing
print(f'M_ej / eta = {(M_ej).to(u.g):.2g}')
print(f'l_sh = {l_sh:.1g}, r_sh = {r_sh:.1g}')

# strength of magnetic field
B2 = calc_mag_B2(e_B,M_ej,v_sh,r_sh,l_sh)
B = np.sqrt(B2)
print(f'B^2 = {B2:.2g}, B= {B:.2g}')

# mass density of ambient medium
rho = calc_rho(M_ej,r_sh,v_sh,tau)
print(f'rho_amb = {rho.to(u.g/(u.cm)**3):3g}')

# number density of ambient medium
# assuming all media are H

n_H = rho/m_u
print(f'n_amb(H) = {n_H.to(u.cm**(-3)):.3g}')


# C = M_ej*v_sh**2 / (2*m_e *c**2)
# print(C)

# minimum Lorentz factor
zeta_e = 0.4
gamma_factor_min = (zeta_e * m_p * v_sh**2) / (2* m_e * c**2)
print(f'gamma_factor_min = {gamma_factor_min:3g}')

# synchrotron radiation
def calc_frequency(gamma_factor, B):
    value = (gamma_factor**2*e*B) / (2*PI*m_e*c)
    return value

nu_e = calc_frequency(gamma_factor_min,B)
print(nu_e.to(u.GHz))

# synchrotron cooling
u_B = B2/(8*PI)
print(f'u_B = {u_B:3g}')
def calc_synP(gamma_factor,u_B):
    # beta2 = 1 - gamma_factor**(-2)
    # print(beta2)
    value = (4*sigma_T*c*gamma_factor**2*u_B) / 3.0
    return value
cooling_rate = calc_synP(gamma_factor_min,u_B)
print(f'P = {cooling_rate:2g} = {cooling_rate.to(u.MeV/u.s):2g}')

cooling_time = (gamma_factor_min*m_e*c**2) / cooling_rate
print(f't_cool = {cooling_time.to(u.min):.3g}')

print(f'Thomson cross-section: {sigma_T:.2g}')