from imports import *
from parameter import *
from calc_shock import r_sh,l_sh
from calc_B_field import B
from calc_powerlaw import N0
from calc_Lorentz_factor import gamma_min

# power(erg cm-3 s-1 Hz-1)
def calc_power(p,B,N0,nu):
    eB = (e*B).to(u.erg/u.cm)
    x = (nu/c).to(u.cm**(-1))
    G = gamma((3*p+19)/12)*gamma((3*p-1)/12)

    val = \
        (np.sqrt(3)*(N0*r_e*eB).to(u.erg*u.cm**(-3))) / (p+1) \
        *((2*PI*mc2*x)/(3*eB))**(-0.5*(p-1)) \
        *G
    return val

# intensity(erg cm-2 s-1 Hz-1 str-1)
def calc_intensity(p,B,N0,l,nu):
    P = calc_power(p,B,N0,nu)
    val = P*l / (4*PI)
    return val

# absorption coefficient(cm-1)
def calc_alpha(p,B,N0,nu):
    eB  = (e*B).to(u.erg/u.cm)
    g1  = gamma((3*p+2)/12)
    g2  = gamma((3*p+22)/12)
    x   = nu/c

    val = \
        ((np.sqrt(3) * N0 * r_e * eB) / (8*PI*mc2*x**2)).to(u.cm**(-1)) \
        * ((3*eB) / (2*PI*mc2*x))**(p/2) \
        *g1*g2
    return val.to(u.cm**(-1))

# source function in optically-thick region(erg cm-2 s-1 Hz-1 str-1)
def calc_SSA_source_function(p,B,nu):
    eB  = (e*B).to(u.erg/u.cm)
    x   = (nu/c).to(u.cm**(-1))
    G   = \
        ( gamma((3*p+19)/12)*gamma((3*p-1)/12) ) \
        /( gamma((3*p+2)/12) *gamma((3*p+22)/12) )
    val = \
        (2*mc2*x**(5/2)) / (p+1) \
        *np.sqrt((2*PI*mc2) / (3*eB)) \
        *G
    return val

# optical depth
def calc_tau(p,B,N0,lsh,nu):
    a = calc_alpha(p,B,N0,nu)
    return a*lsh

# SSA frequency
def calc_nu_SSA(p,B,N0,lsh):
    eB  = (e*B).to(u.erg/u.cm)
    g1  = gamma((3*p+2)/12)
    g2  = gamma((3*p+22)/12)
    val = \
        ((np.sqrt(3) * N0 * r_e * eB *lsh * g1 * g2) / (8*PI*mc2))**(2/(p+4)) \
        *((3*eB) / (2*PI*mc2))**(p/(p+4)) \
        *c
    return val.to(u.GHz)

# converting intensity into flux at observer
def convert_into_flux(I,r,d):
    Omega = (r/d)**2 # solid angle
    return (I*Omega).to(u.mJy)

nu_sample   = 100*u.GHz
P_crit      = calc_power(p,B,N0,nu_sample)
I_crit      = calc_intensity(p,B,N0,l_sh,nu_sample)
alpha_crit  = calc_alpha(p,B,N0,nu_sample)
S_crit      = calc_SSA_source_function(p,B,nu_sample)
nu_SSA      = calc_nu_SSA(p,B,N0,l_sh)
print(
    f'P_crit = {P_crit:.2g}\n'
    f'I_crit = {I_crit:.2g}\n'
    f'a_crit = {alpha_crit:.2g}\n'
    f'S_crit = {S_crit:.2g}\n'
    f'nu_SSA = {nu_SSA:.2g}'
    )

# making dataframe
nu_min = 1.0*u.GHz
nu_max = 1000*u.GHz
nu_full:list    = np.linspace(nu_min,nu_max,1000)
nu_thick:list   = np.linspace(nu_min,nu_SSA,500)
nu_thin:list    = np.linspace(nu_SSA,nu_max,500)

mask_thick  = nu_full < nu_SSA
mask_thin   = nu_full > nu_SSA

intensity:list = [ \
    calc_intensity(p,B,N0,l_sh,nu).to(u.mJy) \
    for nu in nu_full[mask_thin] \
    ]
intensity_dashed:list = [ \
    calc_intensity(p,B,N0,l_sh,nu).to(u.mJy) \
    for nu in nu_full[mask_thick] \
    ]
source_func:list = [ \
    calc_SSA_source_function(p,B,nu).to(u.mJy) \
    for nu in nu_full[mask_thick] \
    ]
flux_thin:list = [ \
    convert_into_flux(I,r_sh,d) \
    for I in intensity \
    ]
flux_dashed:list = [ \
    convert_into_flux(I,r_sh,d) \
    for I in intensity_dashed \
    ]
flux_thick:list = [ \
    convert_into_flux(I,r_sh,d) \
    for I in source_func \
    ]
optical_depth:list = [ \
    calc_tau(p,B,N0,l_sh,nu) \
    for nu in nu_full
]

df_flux_thick   = pd.DataFrame({
    'nu':nu_full[mask_thick],
    'F' :flux_thick
})
df_flux_thin    = pd.DataFrame({
    'nu':nu_full[mask_thin],
    'F' :flux_thin
})
df_flux_dashed  = pd.DataFrame({
    'nu':nu_full[mask_thick],
    'F' :flux_dashed
})

df_thick = pd.DataFrame({
    'nu':nu_full[mask_thick],
    'S':source_func
})

df_thin = pd.DataFrame({
    'nu':nu_full[mask_thin],
    'I':intensity
})
