import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy import units as u
from astropy import constants as const
import numpy as np
from scipy.special import kv
from scipy.special import gamma
from scipy.integrate import quad
from fractions import Fraction
import numpy as np

FONTSIZE_DEFAULT = 18
OUTPUT_DIR = './' # output file's path

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
    B2 = ((eps_B*M_ej*v_sh**2) / (r_sh**2*l_sh)).to(u.erg*u.cm**(-3))
    return B2

def calc_rho(M_ej,r_sh,v_sh,tau):
    value = M_ej / (4*PI*r_sh**2*v_sh*tau)
    return value

def convert_P_to_F(power,distance):
    val = power / (4*PI*distance**2)
    return val.to(u.Jy)


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
print(f'l_sh = {l_sh:.2g}, r_sh = {r_sh:.2g}')

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

# maximum Lorentz factor
def calc_max_gamma(B):
    val = ((9.0*PI*e) / (10*sigma_T*B)).decompose()**(0.5)
    return val

print(f'gamma_max = {calc_max_gamma(B):.2g}')

ratio_ej_to_rest = M_ej*v_sh**2 / (2*m_e *c**2)
print(f'E_ej / m_e c^2 = {ratio_ej_to_rest.decompose()}')
def calc_synchrotron_coefficient(p,eps_e,min_gamma):
    val = (p-2)*min_gamma**(p-2)*eps_e*ratio_ej_to_rest
    return val

# minimum Lorentz factor
zeta_e = 0.4
gamma_factor_min = (zeta_e * m_p * v_sh**2) / (2* m_e * c**2)
print(f'gamma_factor_min = {gamma_factor_min:3g}')

# synchrotron radiation
def calc_frequency(gamma_factor, B):
    value = (gamma_factor**2*e*B) / (2*PI*m_e*c)
    return value

nu_e = calc_frequency(gamma_factor_min,B)
print(f'nu_e = {nu_e.to(u.GHz):.2g}')

# synchrotron cooling
u_B = B2/(8*PI)
print(f'u_B = {u_B:3g}')
def calc_synP(gamma_factor,u_B):
    # beta2 = 1 - gamma_factor**(-2)
    # print(beta2)
    value = (4*sigma_T*c*gamma_factor**2*u_B) / 3.0
    return value
cooling_rate = calc_synP(gamma_factor_min,u_B)
print(f'P = {cooling_rate:.2g} = {cooling_rate.to(u.keV/u.s):.2g}')

cooling_time = (gamma_factor_min*m_e*c**2) / cooling_rate
print(f't_cool = {cooling_time.to(u.min):.3g}')

# total number of electron swept by shock
N_e = (M_ej / m_u).decompose()
# N_e = N_etot / (PI*gamma_factor_min) # beaming effect
print(f'N_e = {N_e:.2g}')

#
P_tot = N_e * cooling_rate
print(f'P_tot = {P_tot:.2g}')

# distance from 4U 1728-34
d = 4.5*u.kpc
print(f'd = {d:.2g} = {d.to(u.cm):.2g}')

# flux estimation
flux = P_tot / (4*PI*d**2)
print(f'F = {flux.to(u.erg*u.cm**(-2)*u.s**(-1)):.2g}')
print(f'F_nu = {(flux/nu_e).to(u.mJy):.2g}')

print(f'Thomson cross-section: {sigma_T:.2g}')

# calculating power spectrum
# v = float(Fraction(5,3)) # order of kv
# def integrand(xi):
#     return kv(v,xi)

# def F(x):
#     val, err = quad(integrand,x,np.inf)
#     return x*val

# def calc_nu_c(gamma_factor,B):
#     E = e*B
#     print(E.to(u.erg*u.em**(-1)))
#     val = (3.0*gamma_factor**2*E) / (4.0*PI*m_e*c)
#     return val

# def P(x,gamma_factor,B):
#     # nu_c = calc_nu_c(gamma_factor,B)
#     val = (np.sqrt(3)*e**3*B*F(x)) / (m_e*c**2)
#     return val

# nu_c = calc_nu_c(gamma_factor_min,B)
# print(f'nu_c = {calc_nu_c(gamma_factor_min,B):.2g}')

# nu_range = np.linspace(1.0,1000,1000)
# # list_x = [list for list in (nu_range/nu_c) if list < 4.0]
# list_x = (nu_range/nu_c) 
# # list_x = [list for list in list_x if list < 4.0]
# # list_P = [P(nu,gamma_factor_min,B).value for nu in nu_range]
# list_P = [P(x,gamma_factor_min,B).value for x in list_x]
# list_I = list_P/nu_range
# df = pd.DataFrame({'nu':nu_range, 'P(nu)':list_P})
# print(df)

def calc_P_nu_noSSA(B,p,nu,eps_e,min_gamma):
    C = calc_synchrotron_coefficient(p,eps_e,min_gamma).decompose()
    eB = (e*B).to((u.erg/u.cm))
    mc2= (m_e * c**2).to(u.erg)
    
    f0 = ((np.sqrt(3) * e**2 * eB)/(mc2)).to(u.erg)
    f1 = C/(p+1)
    f2 = ((3.0 * eB) / (2.0*PI*m_e*c*nu))**((p-1)*0.5)
    g1 = gamma((3*p+19)/12)
    g2 = gamma((3*p-1)/12)

    val = f0*f1*f2*g1*g2
    return val
    # return val.to(u.erg)

def calc_SSA_alpha(p,B,C):
    eB = (e*B).to(u.erg/u.cm)
    mc2 = (m_e * c**2).to(u.erg)
    g1 = gamma((3*p+2)/12)
    g2 = gamma((3*p+22)/12)
    f1 = 3**(p+1) * e**6 * c**(p+4) * (eB)**p \
        / (16 * (2*PI)**(p+2) * (mc2)**(3*p+2))

p=2.5
eps_e = 0.4
min_nu = 100*u.MHz
max_nu = 100*u.GHz
list_nu = np.linspace(min_nu,max_nu,10000)
# print(convert_P_to_F(calc_P_nu_noSSA(B,p,max_nu,eps_e,gamma_factor_min),d))

list_F_nu=[convert_P_to_F(calc_P_nu_noSSA(B,p,nu,eps_e,gamma_factor_min),d).value \
            for nu in list_nu]
# print(list_P_nu)
df = pd.DataFrame({'nu(GHz)':list_nu, 'F(Jy)':list_F_nu})
print(df)
plt.rcParams["font.size"] = FONTSIZE_DEFAULT
fix, ax = plt.subplots(figsize = (8,8))
# ax.set_xlim()
ax.loglog()
ax.set_xlabel(r'$\nu$(GHz)', size=FONTSIZE_DEFAULT)
ax.set_ylabel(r'$F_\nu $(Jy)', size=FONTSIZE_DEFAULT)

ax.plot(list_nu,list_F_nu, ls='-', lw=4)
plt.savefig(OUTPUT_DIR+'syn-noSSA-p'+f'{int(p*10)}'+'.png')

# print(calc_P_nu_noSSA(10*44,B,2,nu_e))
# print(((np.sqrt(3)*e**2)/(m_e*c**2)).to(u.fm))
# print((e*B).to(u.erg/u.cm))
# print((((np.sqrt(3)*e**2)/(m_e*c**2))*e*B).to(u.erg))
# print((((2*PI*m_e*c)/(3*e*B))**(-0.5)).decompose())