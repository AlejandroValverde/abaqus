import math			# necessary for certain mathematical functions
import os			# necessary for creation of folders

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import string
import numpy as np
# from parameters_Dom import N, c_0#, b

from functions_Urban import *

#### PARAMETERS - 1
iSIM=1
PositionFrontSpar=2.000000e-01
PositionRearSpar=3.000000e-01
FaserWinkel=0
thickness_bucklingSpar=1.000000e-04

#===============================================================================
# Parameters DOM
#===============================================================================
import numpy as np
# Check with consistency in 01_Variablendefinition_Dom
# Define parameters
span=3.0 ##b = float(5)                    # Span [m]
c_0 = 0.5                       # Chordlength [m]
V = 60.0#1.5*30*np.sqrt(2)               # Freestream velocity [m/s]
rho = 1.225                     # Density [kg/m^3] from U.S Standard Atmosphere Air Properties in SI Units
mu = 1.789E-5                   # Dynamic viscosity [Ns/m^2]
Re = rho*V*c_0/mu
alpha_init = 8.0#12.0#7.0                # Initial angle of attack in deg
Re_scal = round(Re/c_0)         # Scaled Reynold's number for Xfoil
##alpha_min = 0.0                 # Start value of alpha for sweep                 
##alpha_max = 10.0                # End value of value for sweep
##step = 0.5                      # Increment by which alpha is swept
N = 25                          # Number of half-spanwise elements; not higher than 40, 25 recommended
##D = 0.05                        # Damping factor
##tol = 0.00001                   # Convergence tolerance

##############################################################
##############################################################

#### Variables init

start_dir = folder_path=os.getcwd()
#folder_path = 'D:/faselu/ftero_opt_local'
[folder_path,folder_path_ende]=os.path.split(start_dir)
xfoil_path = folder_path+'/xfoil6.96/bin/xfoil.exe'

##################################################
jobname = 'Wing-stat'
res_name = '/xy_result.txt'
res_name_spars = '/SparVorne.txt' 
res_name_spars_hinten ='/SparHinten.txt'

daempfung=[25000]#[500,600,700]
jobNumber=1


# Calculation of control points
# y_C = y_C_fc()
width_element = span/2/N
y = (2*N+1)
y = np.zeros(y)
y[0] = -span/2
for i in range(2*N):
    y[i+1] = y[0]+width_element*(i+1)        
n = len(y)
y.tolist()
##y remains nd-array
y_C = y[0:n-1]+0.5*width_element
##only takes part of y which including y[0] but excluding y[n-1]
##this plus operation doesn't work for the list datatype because list can only be concatenated with a list


# SI-units

##Spar_Angle = -30
# alf = 0 
Anz_Unterteilung = N                                                            # Check for consisteny with N in parameters
c = c_0                                                                         # chord length [m], Check for consisteny with c_0 in parameters
##rectangular wing, the chord lengths are the same along the spannwise direction
width = span/2                                                                     # width of the wing (half of the spanwidth) [m], check for consisteny with b in parameters
x_spar1 = PositionFrontSpar*c     #0.25*c     #0.16*c     # 0.21*c                                   # position of spar on chord
x_spar2 = PositionRearSpar*c#(PositionFrontSpar+PositionRearSpar+AngleRearSpar)*c     #0.35*c     #0.26*c     # 0.31*c                          # position of spar on chord
x_spar1m = PositionFrontSpar*c     #(PositionFrontSpar+AngleFrontSpar)*c                               # position of spar on chord
x_spar2m = PositionRearSpar*c#(PositionFrontSpar+PositionRearSpar)*c                               # position of spar on chord

#x_spar1 = 0.14*c
#x_spar1m = x_spar1
#x_spar2m = x_spar1 + 0.05
#x_spar2 = x_spar2m #+ 0.01

anzRibs = 6
t_stringer = 0.001
element_width = width/N

# Environmental data: ISA (international standard atmosphere)

# a = 340.29			# speed of sound [m/s]
##rho = 1.225    		# density [kg/m^3]  1.112 at 1000 m
# p_u = 0         	        # ambient pressure [N/m^2]=[Pa] 101325	
# Ma = 0.125
##V = 30*np.sqrt(2)                        # [m/s]

x_str1 = x_spar2+ (c-x_spar2)/4	                # position of spar on chord
x_str2 = x_spar2+ (c-x_spar2)/2           # position of spar on chord1
x_str3 = c-0.02*c #x_spar2+ 3*(c-x_spar2)/4