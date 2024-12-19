from calc_SSA_spectrum import *

plt.rcParams["font.size"] = FONTSIZE_DEFAULT
fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlabel(r'$\nu$ (GHz)')
ax.set_ylabel(r'Flux (mJy)')
ax.loglog()
ax.grid(True, which='both')

# I_mJy:list          = [I.to(u.mJy) for I in df_thin['I']]
# I_mJy_dashed:list   = [I.to(u.mJy) for I in ]
# S_mJy:list          = [S.to(u.mJy) for S in df_thick['S']]

df_plt_thick = pd.DataFrame({ 
    'F':[F.value for F in flux_thick],
    'nu':nu_full[mask_thick]
})
df_plt_thin = pd.DataFrame({ 
    'F':[F.value for F in flux_thin],
    'nu':nu_full[mask_thin]
})
df_plt_dashed = pd.DataFrame({ 
    'F':[F.value for F in flux_dashed],
    'nu':nu_full[mask_thick]
})

ax.plot( \
    df_plt_thin['nu'],
    df_plt_thin['F'],
    ls='-', lw=3,
    color=RED
)

ax.plot( \
    df_plt_thick['nu'],
    df_plt_thick['F'],
    ls='-', lw=3,
    color=BLUE
)

ax.plot( \
    df_plt_dashed['nu'],
    df_plt_dashed['F'],
    ls='--', lw=3,
    color=RED
)

FILENAME = \
    './output/plt_flux' \
    f'-p{int(p*10):02}' \
    f'-h{int(eta*100):03}' \
    '.png'

plt.savefig(FILENAME, bbox_inches='tight')