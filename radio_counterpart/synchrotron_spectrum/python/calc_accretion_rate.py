from imports import *
from parameter import eps

def calc_L_Edd(mass):
    value = 4*PI*G*mass*m_p*c / sigma_T
    return value.to(u.erg / u.year)

def calc_Mdot_Edd(mass):
    L_Edd = calc_L_Edd(mass)
    value = L_Edd / c**2
    return value.to(u.g/u.s)

Mdot = eps * calc_Mdot_Edd(1.4*M)
print(f'Mdot = {Mdot:.3g}')