import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy import units as u
from astropy import constants as const
import numpy as np
from scipy.special import kv
from scipy.special import gamma
from scipy.integrate import quad
from fractions import Fraction
import numpy as np

FONTSIZE_DEFAULT = 18
OUTPUT_DIR = './' # output file's path
RED = '#FF4B00'
GREEN = '#03AF7A'
BLUE = '#005AFF'
CYAN = '#4DC4FF'

# physical and mathematical constants
PI = np.pi
c = const.c.cgs
G = const.G.cgs
M = const.M_sun.cgs
e = const.e.gauss
hbarc = const.hbar * c
a = const.alpha
m_e = const.m_e.cgs
m_p = const.m_p.cgs
m_u = const.u.cgs
mc2 = (m_e*c**2).to(u.erg)
sigma_T = const.sigma_T.cgs
r_e = (e**2 / mc2).to(u.cm)