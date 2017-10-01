#===============================================================================
# Parameters
#===============================================================================
import numpy as np
from moduleBuildWingBoxAndPostProc import *
from moduleCommon import *

class structtype():
    pass

import numpy as np
paraRead_wing = structtype()
paraRead_wing = loadParameters(paraRead_wing, 'inputAbaqus_wing.txt')
# Check with consistency in 01_Variablendefinition_Dom
# Define parameters
span=float(paraRead_wing.span) ##b = float(5)                    # Span [m]
c_0 = float(paraRead_wing.c_0)                       # Chordlength [m]
V = float(paraRead_wing.V)#1.5*30*np.sqrt(2)               # Freestream velocity [m/s]
rho = float(paraRead_wing.rho)                   # Density [kg/m^3] from U.S Standard Atmosphere Air Properties in SI Units
mu = float(paraRead_wing.mu)               # Dynamic viscosity [Ns/m^2]
Re = rho*V*c_0/mu
alpha_init = 8.0#12.0#7.0                # Initial angle of attack in deg
Re_scal = round(Re/c_0)         # Scaled Reynold's number for Xfoil
##alpha_min = 0.0                 # Start value of alpha for sweep                 
##alpha_max = 10.0                # End value of value for sweep
##step = 0.5                      # Increment by which alpha is swept
N = int(paraRead_wing.N)                          # Number of half-spanwise elements; not higher than 40, 25 recommended
##D = 0.05                        # Damping factor
##tol = 0.00001                   # Convergence tolerance
