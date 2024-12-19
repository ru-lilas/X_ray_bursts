from imports import *

# observational value
v_sh = 0.3*c # shock velocity
t_rec = 200 * u.min # reccurence time of X-ray burst
t_rise = 3 *u.min # rising time of radio burst
t_x = 10 * u.s # duration time of X-ray burst
tau = 10 * u.min # decay time of radio burst
d = 4.5*u.kpc # distance from 4U 1728-34

# parameters
eta = float(input(">>eta:")) # M_ej / (Mdot t_rec)
eps = 0.1 # Mdot / Mdot(Edd)
eps_B = 0.01 
eps_e = float(input(">>eps_e:"))
p = float(input(">>power-law index:")) # power-law index
zeta_e = 0.4

# observation
flux90 = [
1.18,
1.16,
1.18,
1.10,
1.48,
1.18,
1.18,
]
flux55 = [
1.10,
1.30,
1.10,
0.98,
1.32,
1.02,
1.06,
]
band55 = [5.5 for _ in range(len(flux55))]
band90 = [9.0 for _ in range(len(flux55))]