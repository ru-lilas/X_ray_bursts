from calc_SSA_spectrum import *

plt.rcParams["font.size"] = FONTSIZE_DEFAULT
fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlabel(r'$\nu$ (GHz)')
ax.set_ylabel(r'Intensity (mJy str$^{-1}$)')
ax.loglog()
ax.grid(True, which='both')

I_mJy:list = [I.to(u.mJy) for I in df['I']]
S_mJy:list = [S.to(u.mJy) for S in df['S']]
df_plt = pd.DataFrame({ 
    'I':[I.value for I in I_mJy],
    'S':[S.value for S in S_mJy],
    'nu':df['nu']
})

ax.plot( \
    df_plt['nu'],
    df_plt['I'],
    ls='-', lw=3,
    color=RED
)

ax.plot( \
    df_plt['nu'],
    df_plt['S'],
    ls='-', lw=3,
    color=BLUE
)

FILENAME = \
    './output/plt_intensity' \
    f'-p{int(p*10):02}' \
    f'-h{int(eta*100):03}' \
    '.png'

plt.savefig(FILENAME, bbox_inches='tight')