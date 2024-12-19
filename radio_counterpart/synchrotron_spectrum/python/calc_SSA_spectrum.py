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

nu_sample   = 100*u.GHz
P_crit      = calc_power(p,B,N0,nu_sample)
I_crit      = calc_intensity(p,B,N0,l_sh,nu_sample)
alpha_crit  = calc_alpha(p,B,N0,nu_sample)
S_crit      = calc_SSA_source_function(p,B,nu_sample)
print(
    f'P_crit = {P_crit:.2g}\n'
    f'I_crit = {I_crit:.2g}\n'
    f'a_crit = {alpha_crit:.2g}\n'
    f'S_crit = {S_crit:.2g}\n'
    )

# making dataframe
nu_min = 1.0*u.GHz
nu_max = 1000*u.GHz
nu_full:list = np.linspace(nu_min,nu_max,1000)

intensity:list = [ \
    calc_intensity(p,B,N0,l_sh,nu) \
    for nu in nu_full \
    ]
source_func:list = [ \
    calc_SSA_source_function(p,B,nu) \
    for nu in nu_full \
    ]

df = pd.DataFrame({
    'nu':nu_full,
    'I':intensity,
    'S':source_func
})

print(df)

# def calc_noSSA_flux(P,r_sh,l_sh,d):
#     # print(f'F/P = {(l_sh*(r_sh/d)**2).to(u.cm)}')
#     val = P*l_sh*(r_sh/d)**2
#     return val.to(u.mJy)

# # source flux
# def calc_SSA_source(p,nu,B):
#     f  = (e*B).to(u.erg/u.cm)
#     x = nu/c

#     # G:ratio of gamma functions
#     G = \
#         ( gamma((3*p+19)/12)*gamma((3*p-1)/12) ) \
#         *( gamma((3*p+2)/12) *gamma((3*p+22)/12) )
#     val = \
#     (8*PI*x**2*mc2) / (p+1) \
#     *np.sqrt(2*PI*mc2*x/(3*f)) \
#     *G
#     return val.to(u.erg*u.cm**(-2))
# def calc_SSA_flux(S,r_sh,d):
#     return (S*(r_sh/d)**2).to(u.mJy)


# #  the frequency at (Optical depth) = 1
# def calc_nu0(p,N0,B,l_sh):
#     eB = (e*B).to(u.erg/u.cm)
#     r_e = (e**2 / (mc2)).to(u.cm)
#     G = \
#         ( gamma((3*p+19)/12)*gamma((3*p-1)/12) )
#     val = \
#         ((N0*eB*l_sh*np.sqrt(3)*r_e*G) / (8*PI*mc2))**(2/(p+4)) \
#         *((3*eB)/(2*PI*mc2))**(p/(p+4)) \
#         *c
#     return val.to(u.GHz)

# def calc_tau(p,N0,B,l_sh,nu):
#     eB = (e*B).to(u.erg/u.cm)
#     r_e = (e**2 / (mc2)).to(u.cm)
#     G = \
#         ( gamma((3*p+19)/12)*gamma((3*p-1)/12) )
#     val = \
#         (np.sqrt(3)*N0*r_e*eB*l_sh) / ((nu/c)**2*mc2) \
#         *(3*eB/(2*PI*mc2*(nu/c)))**(p/2) \
#         *G
#     return val

# def calc_decayfactor(tau):
#     return (1-np.exp(-tau))

# def calc_pow_decayed(p,N0,B,l_sh,nu):
#     P = calc_spectrum_pl(p,B,N0,nu)
#     tau = calc_tau(p,N0,B,l_sh,nu)
#     f = calc_decayfactor(tau)
#     return P*f

# nu_a = calc_nu0(p,N0,B,l_sh) # SSA frequency

# # calculate more concrete case
# nu_sample = 1000*u.GHz

# P_crit      = calc_spectrum_pl(p,B,N0,nu_sample)
# F_crit      = calc_noSSA_flux(P_crit,r_sh,l_sh,d)
# alpha_crit  = calc_SSA_alpha(p,B,N0,nu_sample)
# tau_crit    = calc_tau(p,N0,B,l_sh,nu_sample)
# decay_crit  = calc_decayfactor(tau_crit)
# P_crit_dec  = P_crit*decay_crit
# F_crit_dec  = calc_noSSA_flux(P_crit_dec,r_sh,l_sh,d)
# S_crit      = calc_SSA_source(p,nu_sample,B)
# F_crit_SSA  = calc_SSA_flux(S_crit,r_sh,d)
# print(
#     f'P_crit = {P_crit:.2g} \n' \
#     f'P_crit_decayed = {P_crit_dec:.2g} \n'
#     f'F_crit = {(F_crit):.2g} = {(F_crit.to(u.erg/u.cm**2)):.2g}\n'
#     f'F_crit_decayed = {F_crit_dec:.2g} = {(F_crit_dec.to(u.erg/u.cm**2)):.2g}\n'
#     f'alpha_nu = {alpha_crit:.2g} \n'
#     f'S_crit = {S_crit:.2g} \n'
#     f'F_crit(SSA) = {F_crit_SSA:.2g}'
#     )

# print(f'nu_a = {nu_a:.2g}')

# # nu_min  = e*B*gamma_min**2 / (2*PI*m_e*c) # frequency at minimum Lorentz factor
# nu_min = 1*u.GHz
# nu_max  = 1000 *u.GHz # arbitrary value

# nu_SSA \
# = np.linspace(nu_min, nu_a, 100) # frequency list
# nu_noSSA \
# = np.linspace(nu_a,nu_max,100) # frequency list
# nu_full \
# = np.linspace(nu_min,nu_max,100)

# # source_noSSA:list   = \
# # [
# #     calc_spectrum_pl(p,B,N0,nu)*calc_decayfactor(calc_tau(p,N0,B,l_sh,nu)) \
# #     / calc_SSA_alpha(p,B,N0,nu) \
# #     for nu in nu_noSSA
# # ]
# flux_source:list = [calc_spectrum_pl(p,B,N0,nu)/calc_SSA_alpha(p,B,N0,nu) for nu in nu_full]
# # flux_source_SSA:list = [4*PI*calc_SSA_source(p,nu,B)*calc_decayfactor(calc_tau(p,N0,B,l_sh,nu)) for nu in nu_noSSA]
# flux_source_SSA:list = [calc_spectrum_pl(p,B,N0,nu)*calc_decayfactor(calc_tau(p,N0,B,l_sh,nu))/calc_SSA_alpha(p,B,N0,nu) for nu in nu_noSSA]
# flux_Pl:list = [4*PI*calc_SSA_source(p,nu,B)*calc_SSA_alpha(p,B,N0,nu)*l_sh for nu in nu_noSSA]

# # power_noSSA_nodecay:list    = [calc_spectrum_pl(p,B,N0,nu) for nu in nu_noSSA]
# # power_noSSA:list    = [calc_pow_decayed(p,N0,B,l_sh,nu) for nu in nu_noSSA]
# # intensity_noSSA_nodecay:list = [P*l_sh for P in power_noSSA_nodecay]
# # df = pd.DataFrame({'nu':nu_noSSA, 'I':intensity_noSSA_nodecay})
# # print(df)
# # flux_noSSA_nodecay:list     = [calc_noSSA_flux(P,r_sh,l_sh,d) for P in power_noSSA_nodecay]
# # flux_noSSA:list     = [calc_noSSA_flux(P,r_sh,l_sh,d) for P in power_noSSA]
# # source_noSSA:list   = [calc_SSA_source(p,nu,B)*calc_decayfactor(calc_tau(p,N0,B,l_sh,nu)) for nu in nu_noSSA]

# # source_SSA:list     = [calc_SSA_source(p,nu,B) for nu in nu_SSA]
# # flux_SSA:list       = [calc_SSA_flux(S,r_sh,d) for S in source_SSA]

# # df_noSSA    = pd.DataFrame({'nu':nu_noSSA, 'flux':flux_noSSA})
# # df_SSA      = pd.DataFrame({'nu':nu_SSA, 'flux':flux_SSA})

# # print(df_noSSA)
