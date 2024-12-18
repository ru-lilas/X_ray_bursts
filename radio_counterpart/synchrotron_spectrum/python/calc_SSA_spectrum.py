from imports import *
from parameter import *
from calc import p, r_sh,d, l_sh
from calc_B_field import B
from calc_powerlaw import N0

RED = '#FF4B00'

# power spectrum for power-law distribution
def calc_spectrum_pl(p,B,N0,nu):
    eB = (e*B).to(u.erg/u.cm)
    r_e = (e**2 / mc2).to(u.cm)
    x = (nu/c).to(u.cm**(-1))
    G = gamma((3*p+19)/12)*gamma((3*p-1)/12)
    val = \
        (np.sqrt(3)*N0*r_e*eB) / (p+1) \
        *((2*PI*mc2*x)/(3*eB))**(-0.5*(p-1)) \
        *G
    return val

def calc_noSSA_flux(P,r_sh,l_sh,d):
    # print(f'F/P = {(l_sh*(r_sh/d)**2).to(u.cm)}')
    val = P*l_sh*(r_sh/d)**2
    return val

# absorption coefficient
def calc_SSA_alpha(p,B,N0,nu):
    f = (e*B).to(u.erg/u.cm)
    r_e = e**2 / mc2
    x = c/nu
    g1 = gamma((3*p+2)/12)
    g2 = gamma((3*p+22)/12)

    val = \
        (np.sqrt(3)*mc2**(p-1)*N0*r_e*f*x**2)/(8*PI) \
        *((3*f*x)/(2*PI*mc2**3))**(0.5*p) \
        *g1*g2
    return val

# Source function
def calc_SSA_source(p,nu,B):
    # TODO source functionの計算コードを埋める
    f  = (e*B).to(u.erg/u.cm)
    x = nu/c

    # G:ratio of gamma functions
    G = \
        ( gamma((3*p+19)/12)*gamma((3*p-1)/12) ) \
        *( gamma((3*p+2)/12) *gamma((3*p+22)/12) )
    val = \
    (2*x**2*mc2) / (p+1) \
    *np.sqrt(2*PI*mc2*x/(3*f)) \
    *G
    return val

#  the frequency at (Optical depth) = 1
def calc_nu0(p,N0,B,l_sh):
    eB = (e*B).to(u.erg/u.cm)
    r_e = (e**2 / (mc2)).to(u.cm)
    G = \
        ( gamma((3*p+19)/12)*gamma((3*p-1)/12) )
    val = \
        ((N0*eB*l_sh*np.sqrt(3)*r_e*(mc2)**(p-1))*G / (8*PI))**(2/(p+4)) \
        *((3*eB)/(2*PI*mc2**3))**(p/(p+4)) \
        *c
    return val.to(u.GHz)

nu0 = calc_nu0(p,N0,B,l_sh)
print(f'nu0 = {nu0:.2g}')

nu_min = 1.0 * u.GHz
nu_turn = nu0
nu_max = 1000 *u.GHz


nu_list_opthick = np.linspace(nu_min, nu_turn, 100)
nu_list_noSSA = np.linspace(nu_min,nu_max,1000)

power_noSSA = [calc_spectrum_pl(p,B,N0,nu) for nu in nu_list_noSSA]
flux_noSSA = [calc_noSSA_flux(P,r_sh,l_sh,d).to(u.mJy) for P in power_noSSA]

source_func = [calc_SSA_source(p,nu,B).to(u.mJy) for nu in nu_list_opthick]
# S_nu += [ \
#     calc_noSSA_flux(calc_spectrum_pl(p,B,N0,nu),r_sh,l_sh,d).to(u.mJy)\
#     for nu in nu_list_noSSA \
#     ]
flux_SSA = [(2*PI*S*(r_sh/d)**2).to(u.mJy) for S in source_func]
# print(flux_noSSA)
# df = pd.DataFrame({'nu(GHz)':nu_list_opthick, 'F_nu(mJy)':[F.value for F in F_nu]})
# print(df)

