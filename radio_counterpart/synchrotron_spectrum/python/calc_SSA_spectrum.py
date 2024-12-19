from imports import *
from parameter import *
from calc_shock import r_sh,l_sh
from calc_B_field import B
from calc_powerlaw import N0
from calc_Lorentz_factor import gamma_min

# power spectrum for power-law distribution
def calc_spectrum_pl(p,B,N0,nu):
    eB = (e*B).to(u.erg/u.cm)
    x = (nu/c).to(u.cm**(-1))
    G = gamma((3*p+19)/12)*gamma((3*p-1)/12)
    val = \
        (np.sqrt(3)*(N0*r_e*eB).to(u.erg*u.cm**(-3))) / (p+1) \
        *((2*PI*mc2*x)/(3*eB))**(-0.5*(p-1)) \
        *G
    return val

def calc_noSSA_flux(P,r_sh,l_sh,d):
    # print(f'F/P = {(l_sh*(r_sh/d)**2).to(u.cm)}')
    val = P*l_sh*(r_sh/d)**2
    return val.to(u.mJy)

# absorption coefficient
def calc_SSA_alpha(p,B,N0,nu):
    eB = (e*B).to(u.erg/u.cm)
    g1 = gamma((3*p+2)/12)
    g2 = gamma((3*p+22)/12)

    val = \
        ((np.sqrt(3) * N0 * r_e * eB) / (8*PI*mc2*(nu/c)**2)).to(u.cm**(-1)) \
        * ((3*eB) / (2*PI*mc2*(nu/c)))**(p/2) \
        *g1*g2
    return val.to(u.cm**(-1))

# Source function
def calc_SSA_source(p,nu,B):
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
    return val.to(u.erg*u.cm**(-2))

def calc_SSA_flux(S,r_sh,d):
    return (S*(r_sh/d)**2).to(u.mJy)


#  the frequency at (Optical depth) = 1
def calc_nu0(p,N0,B,l_sh):
    eB = (e*B).to(u.erg/u.cm)
    r_e = (e**2 / (mc2)).to(u.cm)
    G = \
        ( gamma((3*p+19)/12)*gamma((3*p-1)/12) )
    val = \
        ((N0*eB*l_sh*np.sqrt(3)*r_e*G) / (8*PI*mc2))**(2/(p+4)) \
        *((3*eB)/(2*PI*mc2))**(p/(p+4)) \
        *c
    return val.to(u.GHz)

def calc_tau(p,N0,B,l_sh,nu):
    eB = (e*B).to(u.erg/u.cm)
    r_e = (e**2 / (mc2)).to(u.cm)
    G = \
        ( gamma((3*p+19)/12)*gamma((3*p-1)/12) )
    val = \
        (np.sqrt(3)*N0*r_e*eB*l_sh) / ((nu/c)**2*mc2) \
        *(3*eB/(2*PI*mc2*(nu/c)))**(p/2) \
        *G
    return val
    

nu_a = calc_nu0(p,N0,B,l_sh) # SSA frequency

# calculate more concrete case
nu_sample = 100*u.GHz

P_crit      = calc_spectrum_pl(p,B,N0,nu_sample)
F_crit      = calc_noSSA_flux(P_crit,r_sh,l_sh,d)
alpha_crit  = calc_SSA_alpha(p,B,N0,nu_sample)
S_crit      = calc_SSA_source(p,nu_sample,B)
F_crit_SSA  = calc_SSA_flux(S_crit,r_sh,d)
print(
    f'P_crit = {P_crit:.2g} \n' \
    f'F_crit = {(F_crit):.2g} = {(F_crit.to(u.erg/u.cm**2)):.2g}\n'
    f'alpha_nu = {alpha_crit:.2g} \n'
    f'S_crit = {S_crit:.2g} \n'
    f'F_crit(SSA) = {F_crit_SSA:.2g}'
    )

print(f'nu_a = {nu_a:.2g}')

# nu_min  = e*B*gamma_min**2 / (2*PI*m_e*c) # frequency at minimum Lorentz factor
nu_min = 1*u.GHz
nu_max  = 1000 *u.GHz # arbitrary value

nu_SSA \
= np.linspace(nu_min, nu_a, 100) # frequency list
nu_noSSA \
= np.linspace(nu_a,nu_max,100) # frequency list

power_noSSA:list    = [calc_spectrum_pl(p,B,N0,nu) for nu in nu_noSSA]
flux_noSSA:list     = [calc_noSSA_flux(P,r_sh,l_sh,d) for P in power_noSSA]

source_SSA:list     = [calc_SSA_source(p,nu,B) for nu in nu_SSA]
flux_SSA:list       = [calc_SSA_flux(S,r_sh,d) for S in source_SSA]

df_noSSA    = pd.DataFrame({'nu':nu_noSSA, 'flux':flux_noSSA})
df_SSA      = pd.DataFrame({'nu':nu_SSA, 'flux':flux_SSA})

print(df_SSA)
