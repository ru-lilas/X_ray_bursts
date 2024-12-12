import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy import units as u
from astropy import constants as const
from scipy.special import kv
from scipy.integrate import quad
from fractions import Fraction
import numpy as np

# physical and mathematical constants
PI = np.pi
c = const.c.cgs
G = const.G.cgs
M = const.M_sun.cgs
e = const.e.gauss
m_e = const.m_e.cgs
m_p = const.m_p.cgs
m_u = const.u.cgs
sigma_T = const.sigma_T.cgs

FONTSIZE_DEFAULT = 18
OUTPUT_DIR = './' # output file's path

v = float(Fraction(5,3)) # order of kv
def integrand(xi):
    return kv(v,xi)

def F(x):
    val, err = quad(integrand,x,np.inf)
    return x*val

def calc_nu_c(gamma_factor,B):
    val = (3.0*gamma_factor**2*e*B) / (4.0*PI*m_e*c)
    return val

def P(nu,gamma_factor,B):
    nu_c = calc_nu_c(gamma_factor,B)
    val = (e**3*B*F(nu/nu_c)) / (m_e*c**2)
    return val 