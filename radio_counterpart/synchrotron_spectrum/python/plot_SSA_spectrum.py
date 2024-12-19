from parameter import *
from calc_SSA_spectrum import plt,FONTSIZE_DEFAULT,RED,GREEN,BLUE,CYAN
from calc_SSA_spectrum import nu_list_noSSA, nu_list_opthick
from calc_SSA_spectrum import flux_noSSA, flux_SSA 

plt.rcParams["font.size"] = FONTSIZE_DEFAULT
fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlabel(r'$\nu$ (GHz)')
ax.set_ylabel(r'Flux (mJy)')
ax.loglog()
ax.grid(True, which='both')

ax.plot( \
    [nu.value for nu in nu_list_opthick], \
    [F.value for F in flux_SSA], \
    ls='-', lw=3, color=RED)

ax.plot( \
    [nu.value for nu in nu_list_noSSA], \
    [F.value for F in flux_noSSA], \
    ls='-', lw=3, color=RED)

FILENAME = \
    './output/SSA_plt' \
    f'-p{int(p*10):02}' \
    f'-h{int(eta*100):03}' \
    '.png'

plt.savefig(FILENAME, bbox_inches='tight')