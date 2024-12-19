from calc_SSA_spectrum import *

plt.rcParams["font.size"] = FONTSIZE_DEFAULT
fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlabel(r'$\nu$ (GHz)')
ax.set_ylabel(r'Flux (mJy)')
ax.loglog()
ax.grid(True, which='both')

# ax.plot( \
#     nu_SSA, \
#     [F.value for F in flux_SSA], \
#     ls='-', lw=3, color=RED)

# ax.plot( \
#     nu_noSSA, \
#     [F.value for F in flux_noSSA], \
#     ls='-', lw=3, color=RED)

# ax.plot( \
#     nu_noSSA, \
#     [F.value for F in flux_noSSA_nodecay], \
#     ls='--', lw=3, color=BLUE)

# ax.plot( \
#     nu_noSSA, \
#     [F.value for F in flux_noSSA_nodecay], \
#     ls='--', lw=3, color=BLUE)
# ax.plot( \
#     nu_noSSA, \
#     [F.value for F in source_noSSA], \
#     ls='-', lw=3, color=BLUE)
# ax.plot( \
#     nu_noSSA, \
#     [I.value for I in intensity_noSSA_nodecay], \
#     ls='--', lw=3, color=RED)
# ax.plot( \
#     [nu for nu in nu_thin[solid].value], \
#     [F.value for F in flux_noSSA_solid], \
#     ls='-', lw=3, color=RED
#     )
ax.plot( \
    nu_full, \
    [F.value for F in flux_source], \
    ls='-', lw=3, color=BLUE)
ax.plot( \
    nu_noSSA, \
    [F.value for F in flux_source_SSA], \
    ls='-', lw=3, color=RED)
ax.plot( \
    nu_noSSA, \
    [F.value for F in flux_Pl], \
    ls='-', lw=3, color=GREEN)


FILENAME = \
    './output/SSA_plt' \
    f'-p{int(p*10):02}' \
    f'-h{int(eta*100):03}' \
    '.png'

plt.savefig(FILENAME, bbox_inches='tight')