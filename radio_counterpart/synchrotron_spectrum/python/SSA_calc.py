from imports import *
from calc import p, r_sh,d
from calc_B_field import B

RED = '#FF4B00'

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

nu_min = 1.0 * u.GHz
nu_max = 10.0 * u.GHz
nu_list = np.linspace(nu_min, nu_max, 100)

print(f'{(((r_sh/d)**2).decompose()):.2g}')

S_nu = [calc_SSA_source(p,nu,B).to(u.mJy) for nu in nu_list]
F_nu = [(2*PI*S*(r_sh/d)**2).to(u.mJy) for S in S_nu]
df = pd.DataFrame({'nu(GHz)':nu_list, 'F_nu(mJy)':[F.value for F in F_nu]})
print(df)

plt.rcParams["font.size"] = FONTSIZE_DEFAULT
fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlabel(r'$\nu$ (GHz)')
ax.set_ylabel(r'Flux (mJy)')
ax.loglog()
ax.grid(True, which='both')

ax.plot(df['nu(GHz)'],df['F_nu(mJy)'], ls='-', lw=3, color=RED)

plt.savefig('./output/SSA_plt.png', bbox_inches='tight')

