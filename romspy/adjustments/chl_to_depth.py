import numpy as np
import math
import numba as nb
from numba.core.types import float64
from scipy.special import erf

"""
Vertical extrapolation of chlorophyll from the surface values
using Morel and Berthon (1989) parameterization

Reference: Morel and Berthon, Surface pigments, algal biomass
        profiles, and potential production of the euphotic layer:
        Relationships reinvestigated in view of remote-sensing 
        applications. Limnol. Oceanogr., 34, 1989, 1545-1562.

Input:

  Cpd : mean pigment concentration within the surface layer
        2D horizontal matrix (shape: yx). UNITS: (mg Chla/m3) if 
        convert=False (default), otherwise (micro mole/l).
  z   : vertical positions (m, negative) 3D matrix with shape zyx
  convert (optional argument) : convert from umol/l to mgChla/m3 if True
                                DO NOT CHANGE UNITS IF False, i.e. if
                                Cpd is already given in mgChla/m3

Output: has shape zyx

  chlo: chlorophyll concentration (mgC.m-3) 3D matrix (convert = True)
  chlo: chlorophyll concentration (mgChla.m-3) 3D matrix (convert = False)
"""

# ufunc version of math.erf
# Note: this function was introduced since using scipy.special.erf in
# the "run" function failed (in nopython mode) with this message:
# Untyped global name 'erf': Cannot determine Numba type of <class 'numpy.ufunc'>
@nb.vectorize([float64(float64)])
def erf_vect(x):
    return math.erf(x)

@nb.njit
def run(Cpd_2d,z,convert=False):
    nz, ny, nx = z.shape
    #Cpd = np.reshape(Cpd, (1,ny,nx))
    #Cpd = np.tile(Cpd,(nz,1,1))
    Cpd = np.empty((nz,ny,nx), dtype=np.float64)
    for i in range(nz):
        Cpd[i,:,:] = Cpd_2d
    if convert:
        # Mass balance for chla molecule: 893.5/660. [mg Chla (mg C)-1]
        chla_C  = 1.3538
        # Molecular mass of chlorophyll [mg/mol]:
        Mchla = 893.5
        Cpd = Mchla*Cpd/1000  # mg/m3
    # Total pigment content within the euphotic layer (Ctot)
    # expressed as the mean pigment concentration within
    # the surface layer (Cpd) (stratified waters)
    Ctot = 38*Cpd**0.425  # Eq. (2b) in Morel & Berthon (1989)
    #Ctot[Cpd>1] = 40.2*Cpd[Cpd>1]**0.507  # Eq. (2c) in Morel & Berthon (1989)
    Ctot = np.where(Cpd>1, 40.2*Cpd**0.507, Ctot)  # Eq. (2c) in Morel & Berthon (1989)
    # Depth of the euphotic layer:
    Ze = -568.2*Ctot**(-0.746)  # Eq. (1a) in Morel & Berthon (1989)
    #Ze[Ze<=-102] = -200.0*Ctot[Ze<=-102]**(-0.293)  # Eq. (1b) in Morel & Berthon (1989)
    Ze = np.where(Ze<=-102, -200.0*Ctot**(-0.293), Ze)  # Eq. (1b) in Morel & Berthon (1989)
    zeta = z/Ze
    # Mean pigment concentration within
    # the euphotic layer (Cze)
    # espressed as the mean pigment concentration within
    # the surface layer (Cpd) (stratified waters)
    Cze = 1.12*Cpd**0.803  # Eq. on top left of p. 1557 in Morel & Berthon (1989)
    # Vertical chlorophyll profile:
    lc = np.log10(Cpd)
    lc2 = lc*lc
    lc3=lc2*lc
    # Equations after (6) on p. 1557 in Morel & Berthon (1989):
    Cb = 0.768+0.087*lc-0.179*lc2-0.025*lc3
    Cmax = 0.299-0.289*lc+0.579*lc2
    zetamax = 0.600-0.640*lc+0.021*lc2+0.115*lc3
    dzeta = 0.710+0.159*lc+0.021*lc2
    # Recompute iteratively the euphotic layer depth:
    eps = 1000
    while eps>0.1:
        Zeold = Ze
        Ctot = -Ze*(Cze*Cb+Cze*Cmax*0.5*np.sqrt(np.pi)*dzeta*\
                (erf_vect((1-zetamax)/dzeta)-erf_vect(-zetamax/dzeta)))
        Ze = -568.2*Ctot**(-0.746)  # Eq. (1a) in Morel & Berthon (1989)
        #Ze[Ze<=-102] = -200*Ctot[Ze<=-102]**(-0.293)  # Eq. (1b) in Morel & Berthon (1989)
        Ze = np.where(Ze<=-102, -200*Ctot**(-0.293), Ze)  # Eq. (1b) in Morel & Berthon (1989)
        Ze = 0.5*(Ze+Zeold)
        zeta = z/Ze
        eps = np.max(np.abs(Zeold-Ze))
    # The second line is eq. (6) in Morel & Berthon (1989)
    C = 0.5*(1+np.tanh(2*(2*Ze-z)/Ze))*\
        Cze*(Cb+Cmax*np.exp(-((zeta-zetamax)/dzeta)**2))
    #C[C<0] = 0
    C = np.where(C<0, 0, C)
    # Compute chlo in mgC if convert is True:
    if convert:
        chlo = C/chla_C
        #chlo[Cpd<0.01] = 0  # WARN HF: I wouldn't do this!
        chlo = np.where(Cpd<0.01, 0, chlo)  # WARN HF: I wouldn't do this!
    else:
        # Otherwise, leave it in mg Chla m-3:
        chlo = C
    return chlo