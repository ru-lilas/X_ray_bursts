from imports import c,u

# observational value
v_sh = 0.3*c # shock velocity
t_rec = 200 * u.min # reccurence time of X-ray burst
t_rise = 3 *u.min # rising time of radio burst
t_x = 10 * u.s # duration time of X-ray burst
tau = 10 * u.min # decay time of radio burst
d = 4.5*u.kpc # distance from 4U 1728-34

# parameters
eta = 1.0 # M_ej / (Mdot t_rec)
eps = 0.1 # Mdot / Mdot(Edd)
eps_B = 0.01 
eps_e = 0.1
p = 2.5 # power-law index
zeta_e = 0.4
