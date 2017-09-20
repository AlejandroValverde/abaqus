import math         # necessary for certain mathematical functions
import os           # necessary for creation of folders
import sys

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
from moduleBuildWingBoxAndPostProc import *

class structtype():
    pass

#===============================================================================
# 0-Variables, clean
#==============================================================================

# - N: Number of half-spanwise elements
# - width_element: Width of the elements in the spanwise direction
# - y, y_C: Control points in the spanwise direction
# - c_0, c: Chord length [m]
# - span: Wing span [m]
# - width: width of the wing (half of the span width) [m]
# - PositionFrontSpar: Dimensional-less (c) position of the front spar
# - PositionRearSpar: Dimensional-less (c) position of the rear spar
# - x_spar1: PositionFrontSpar*c
# - x_spar2: PositionRearSpar*c
# - anzRibs: Number of ribs
# - s_stringer: Stringer thickness
# - element_width: width/N
# - x_str1, x_str2, x_str3: Position of str on chord 1, 2, 3
# - x_spar2: 
# - n_skinpts: Number of points in the skin
# - skin_points_upX: Point vector for the skin 
# - skin_points_upY: Point vector for the skin
# - skin_points_dnX: Point vector for the skin
# - skin_points_dnY: Point vector for the skin
# - mod6_UP_r, mod6_UP_m, mod6_UP_f: Define nodes to find edges in the modeling part of the code, upper surface
# - mod6_DN_m: Define nodes to find edges in the modeling part of the code, lower surface
# - mod7_UP, mod7_DN:
# - mod8_UP, mod8_DN: 
# - x_surf, y_surf: Surface points
# - points_nose, points_box_up, points_box_low, points_trail1_low, points_trail1_up, 
#   points_trail2_low, points_trail2_up, points_trail3_low, points_trail3_up, 
#   points_trail4_low, points_trail4_up: Points defining sections of the airfoil p.e.: x=points_nose[i][0], ypoints_nose[i][1]
# - points_trail_rib_low, points_trail_rib_up: Points defining sections of the rib
# - rib_box_up, rib_box_low: 
# 

#===============================================================================
# 1-Parameters 
#===============================================================================
import numpy as np
paraRead_wing = structtype()
paraRead_wing = loadParameters(paraRead_wing, 'inputAbaqus_wing.txt')

#### PARAMETERS - 1
iSIM=1
PositionFrontSpar = float(paraRead_wing.PositionFrontSpar)
PositionRearSpar = float(paraRead_wing.PositionRearSpar)
FaserWinkel=0
thickness_bucklingSpar=1.000000e-04

# Check with consistency in 01_Variablendefinition_Dom
# Define parameters
span= float(paraRead_wing.span) ##b = float(5)                    # Span [mm]
c_0 = float(paraRead_wing.c_0)                       # Chordlength [mm]
V = float(paraRead_wing.V)#1.5*30*np.sqrt(2)               # Freestream velocity [m/s]
rho = float(paraRead_wing.rho)                     # Density [kg/m^3] from U.S Standard Atmosphere Air Properties in SI Units
mu = float(paraRead_wing.mu)                   # Dynamic viscosity [Ns/m^2]
Re = rho*V*c_0/mu
# alpha_init = 8.0#12.0#7.0                # Initial angle of attack in deg
Re_scal = round(Re/c_0)         # Scaled Reynold's number for Xfoil
##alpha_min = 0.0                 # Start value of alpha for sweep                 
##alpha_max = 10.0                # End value of value for sweep
##step = 0.5                      # Increment by which alpha is swept
N = int(paraRead_wing.N)       # Number of half-spanwise elements; not higher than 40, 25 recommended
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

#Abaqus model
WingModel = mdb.Model(name = 'Model-SparAngle-'+str(int(jobNumber))) 
model = mdb.models['Model-SparAngle-'+str(int(jobNumber))]

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

#===============================================================================
# 2-PointsAndSurfaces
#===============================================================================

#Import airfoil points
Cpxy = np.genfromtxt(start_dir+'/NACA0012_199points.dat', delimiter=None, skip_header=1)

##Cpxy is a 199*2 dimensional array
x = Cpxy[:, 0]
y = Cpxy[:, 1]
##x and y are 199*1 dimensional array
n_skinpts = len(x)
skin_points = [[0 for i in range(2)] for i in range(n_skinpts)]
for i in range(n_skinpts):
    # skin_points[i][0] = x[i]
    # skin_points[i][1] = y[i]
    skin_points[i][0] = x[i]*c_0
    skin_points[i][1] = y[i]*c_0

skin_points_upX = [[0] for i in range(n_skinpts/2)]
skin_points_upY = [[0] for i in range(n_skinpts/2)]
for i in range(n_skinpts/2):
	skin_points_upX[i] = x[i]*c_0
	skin_points_upY[i] = y[i]*c_0

#skin_points_upX.sort(reverse=True)
#skin_points_upY.sort(reverse=True)
#skin_points_upX.sort()
#skin_points_upY.sort()	
skin_points_upX = np.flipud(skin_points_upX) 
skin_points_upY = np.flipud(skin_points_upY) 
	
skin_points_dnX = [[0] for i in range(n_skinpts/2)]
skin_points_dnY = [[0] for i in range(n_skinpts/2)]
for i in range(n_skinpts/2):
	j = i+n_skinpts/2
	skin_points_dnX[i] = x[j]*c_0
	skin_points_dnY[i] = y[j]*c_0

#skin_points_dn = skin_points[n_skinpts/2:-1]
#skin_points_up = skin_points[0:n_skinpts/2]


##define nodes to find edges in 06_modelling...
mod6_UP_r = [0,0];
mod6_UP_m = [0,0];
mod6_UP_f = [0,0];
distSPARsectUP = x_spar2m - x_spar1
distSPARsectDN = x_spar2 - x_spar1m
for i in range(n_skinpts/2):	
	if i > 2 and mod6_UP_r[0] == 0:
		mod6_UP_r[0] = x[i]*c_0;
		mod6_UP_r[1] = y[i]*c_0;
	if x[i]*c_0 < x_spar2m - distSPARsectUP/2 and mod6_UP_m[0] == 0:
		mod6_UP_m[0] = x[i]*c_0;
		mod6_UP_m[1] = y[i]*c_0;
	if x[i]*c_0 < x_spar1/2 and mod6_UP_f[0] == 0:
		mod6_UP_f[0] = x[i]*c_0;
		mod6_UP_f[1] = y[i]*c_0;

mod6_DN_m = [0,0];
for i in range(n_skinpts/2+1,n_skinpts):	
	if x[i]*c_0 > x_spar1m + distSPARsectDN/2 and mod6_DN_m[0] == 0:
		mod6_DN_m[0] = x[i]*c_0;
		mod6_DN_m[1] = y[i]*c_0;


mod7_UP = [x_spar1, numpy.interp(x_spar1, skin_points_upX, skin_points_upY)]
mod7_DN = [x_spar1m, numpy.interp(x_spar1m, skin_points_dnX, skin_points_dnY)]

mod8_UP = [x_spar2m, numpy.interp(x_spar2m, skin_points_upX, skin_points_upY)]
mod8_DN = [x_spar2, numpy.interp(x_spar2, skin_points_dnX, skin_points_dnY)]

		
x_surf=[]
y_surf=[]
anzahl_surfacepoints = 0
for i in range(n_skinpts/2):
    if abs(skin_points[i+1][0]) <= x_str3 and abs(skin_points[i][0]) > x_str3:
        x_surf.append( skin_points[i][0] )
        y_surf.append( skin_points[i][1] )
        anzahl_surfacepoints=anzahl_surfacepoints+1
    elif abs(skin_points[i+1][0]) <= x_spar2 and abs(skin_points[i][0]) > x_spar2:
        x_surf.append( skin_points[i][0] )
        y_surf.append( skin_points[i][1] )
        anzahl_surfacepoints=anzahl_surfacepoints+1
    elif abs(skin_points[i+1][0]) <= x_str2 and abs(skin_points[i][0]) > x_str2:
        x_surf.append( skin_points[i][0] )
        y_surf.append( skin_points[i][1] )
        anzahl_surfacepoints = anzahl_surfacepoints+1
    elif abs(skin_points[i+1][0]) <= x_str1 and abs(skin_points[i][0]) > x_str1:
        x_surf.append( skin_points[i][0] )
        y_surf.append( skin_points[i][1] )
        anzahl_surfacepoints = anzahl_surfacepoints+1
    elif abs(skin_points[i+1][0]) <= x_spar1 and abs(skin_points[i][0]) > x_spar1:
        x_surf.append( skin_points[i][0] )
        y_surf.append( skin_points[i][1] )
        anzahl_surfacepoints = anzahl_surfacepoints+1
        x_surf.append( skin_points[i+2][0] )                                 
        y_surf.append( skin_points[i+2][1] )
        anzahl_surfacepoints = anzahl_surfacepoints+1

points_nose=[]
points_box_up=[]
points_box_low=[]
points_trail1_low=[]
points_trail1_up=[]
points_trail2_low=[]
points_trail2_up=[]
points_trail3_low=[]
points_trail3_up=[]
points_trail4_low=[]
points_trail4_up=[]

lss = 100#100#98		
for i in range(n_skinpts):
#   if (abs(skin_points[i][0]) < x_spar1 and i<lss) or (abs(skin_points[i][0]) < x_spar1m and i>=lss):
#       points_nose.append([skin_points[i][0], skin_points[i][1]])
    if (abs(skin_points[i][0]) <= x_spar1 and i<lss) or (abs(skin_points[i][0]) <= x_spar1m and i>=lss):
        points_nose.append([skin_points[i][0], skin_points[i][1]])
		
    elif abs(skin_points[i][0]) > x_spar1 and abs(skin_points[i][0]) <  x_spar2m and i<lss:#skin_points[i][1] > 0:
        points_box_up.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_spar1m and abs(skin_points[i][0]) <  x_spar2 and i>=lss:#skin_points[i][1] < 0:
        points_box_low.append([skin_points[i][0], skin_points[i][1]])
#    elif abs(skin_points[i][0]) >= x_spar1 and abs(skin_points[i][0]) <  x_spar2m and i<lss:#skin_points[i][1] > 0:
#        points_box_up.append([skin_points[i][0], skin_points[i][1]])
#        
#    elif abs(skin_points[i][0]) >= x_spar1m and abs(skin_points[i][0]) <  x_spar2 and i>=lss:#skin_points[i][1] < 0:
#        points_box_low.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_spar1 and abs(skin_points[i][0]) <  x_str1 and i<lss:#skin_points[i][1] > 0:
        points_trail1_up.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_spar1 and abs(skin_points[i][0]) <  x_str1 and i>=lss:#skin_points[i][1] < 0:
        points_trail1_low.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_str1 and abs(skin_points[i][0]) <  x_str2 and i<lss:#skin_points[i][1] > 0:
        points_trail2_up.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_str1 and abs(skin_points[i][0]) <  x_str2 and i>=lss:#skin_points[i][1] < 0:
        points_trail2_low.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_str2 and abs(skin_points[i][0]) <  x_str3 and i<lss:#skin_points[i][1] > 0:
        points_trail3_up.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_str2 and abs(skin_points[i][0]) <  x_str3 and i>=lss:#skin_points[i][1] < 0:
        points_trail3_low.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_str3 and i<lss:#skin_points[i][1] > 0:
        points_trail4_up.append([skin_points[i][0], skin_points[i][1]])
        
    elif abs(skin_points[i][0]) > x_str3 and i>=lss:#skin_points[i][1] < 0:
        points_trail4_low.append([skin_points[i][0], skin_points[i][1]])


		
points_trail_rib_low=[]
points_trail_rib_up=[]


for i in range(n_skinpts):
    if abs(skin_points[i][0]) > x_spar2m and i<lss:#skin_points[i][1] > 0:
        points_trail_rib_up.append([skin_points[i][0], skin_points[i][1]])       
    elif abs(skin_points[i][0]) > x_spar2 and i>=lss:#skin_points[i][1] < 0:
        points_trail_rib_low.append([skin_points[i][0], skin_points[i][1]])      


points_nose[0] = [x_spar1, numpy.interp(x_spar1, skin_points_upX, skin_points_upY)]
points_nose[-1] = [x_spar1m, numpy.interp(x_spar1m, skin_points_dnX, skin_points_dnY)]
points_box_up[0] = [x_spar2m, numpy.interp(x_spar2m, skin_points_upX, skin_points_upY)]
points_box_low[-1] = [x_spar2, numpy.interp(x_spar2, skin_points_dnX, skin_points_dnY)]

#===============================================================================
# 3-RibTrailOpen
#===============================================================================

## closed rib

rib_nose = [row[:] for row in points_nose]

# Box right

offset_distance=0.003#0.003#0.008
rib_box_up = [row[:] for row in points_box_up]                      
rib_box_up.insert(len(rib_box_up), points_nose[0])

rib_box_low = [row[:] for row in points_box_low]                      
rib_box_low.insert(0, points_nose[-1])
#

############################################################################################################################################################# Trail ###################################################################################################################################################################################################
## Rib Trail
rib_trail_up = [row[:] for row in points_trail_rib_up]                      
rib_trail_up.insert(len(rib_trail_up), points_box_up[0])
rib_trail_low = [row[:] for row in points_trail_rib_low]                      
rib_trail_low.insert(0, points_box_low[-1])

pointsSkin = structtype()

pointsSkin.upX = skin_points_upX
pointsSkin.upY = skin_points_upY
pointsSkin.dnX = skin_points_dnX
pointsSkin.dnY = skin_points_dnY

pointsSkin.rib_box_low = rib_box_low
pointsSkin.rib_box_up = rib_box_up

orthovector_trail_up = [[0 for x in xrange(2)] for x in xrange(len(rib_trail_up))]
for x in xrange(len(rib_trail_up)-1):
                rib_trail_up_offsets_x=rib_trail_up[x+1][0]-rib_trail_up[x][0]
                rib_trail_up_offsets_y=rib_trail_up[x+1][1]-rib_trail_up[x][1]
                orthovector_y=-sqrt(offset_distance*offset_distance/((rib_trail_up_offsets_y*rib_trail_up_offsets_y)/(rib_trail_up_offsets_x*rib_trail_up_offsets_x) + 1))                
                orthovector_x=(-rib_trail_up_offsets_y/rib_trail_up_offsets_x)*orthovector_y               
                orthovector_trail_up[x][0]=orthovector_x+rib_trail_up[x][0]
                orthovector_trail_up[x][1]=orthovector_y+rib_trail_up[x][1]

anzDel=8
for x in xrange(anzDel):
        del orthovector_trail_up[0]

del orthovector_trail_up[-1]

orthovector_trail_low = [[0 for x in xrange(2)] for x in xrange(len(rib_trail_low))]
for x in xrange(len(rib_trail_low)-1):
                rib_trail_low_offsets_x=rib_trail_low[x+1][0]-rib_trail_low[x][0]
                rib_trail_low_offsets_y=rib_trail_low[x+1][1]-rib_trail_low[x][1]
                orthovector_y=sqrt(offset_distance*offset_distance/((rib_trail_low_offsets_y*rib_trail_low_offsets_y)/(rib_trail_low_offsets_x*rib_trail_low_offsets_x) + 1))                
                orthovector_x=(-rib_trail_low_offsets_y/rib_trail_low_offsets_x)*orthovector_y            
                orthovector_trail_low[x][0]=orthovector_x+rib_trail_low[x][0]
                orthovector_trail_low[x][1]=orthovector_y+rib_trail_low[x][1]

del orthovector_trail_low[0]
for x in xrange(anzDel):
        del orthovector_trail_low[-1]

# airfoil thickness at rear spar
airfoilTH = rib_trail_up[-1][1] - rib_trail_low[0][1]
fairfoilTH = 10 #* (c_0/0.5)
airfoilDIST = 0.002 * (c_0/0.5) ### Alejandro

#Build part trail
WingModel.ConstrainedSketch(name='__profile__', sheetSize=2*c)
WingModel.sketches['__profile__'].Spline(points=rib_trail_up[:])    
#WingModel.sketches['__profile__'].Line(point1=rib_trail_up[-1], point2=[orthovector_trail_up[-1][0]+0.004,orthovector_trail_up[-1][1]])
WingModel.sketches['__profile__'].Line(point1=rib_trail_up[-1], point2=[rib_trail_up[-1][0]+airfoilDIST,rib_trail_up[-1][1]-airfoilTH/fairfoilTH])
WingModel.sketches['__profile__'].Line(point1=rib_trail_up[0], point2=rib_trail_low[-1])
WingModel.sketches['__profile__'].Spline(points=rib_trail_low[:])
#WingModel.sketches['__profile__'].Line(point1=rib_trail_low[0], point2=orthovector_trail_low[0])
#WingModel.sketches['__profile__'].Line(point1=orthovector_trail_up[-1], point2=orthovector_trail_low[0])
if x_spar2 == x_spar2m:
	#WingModel.sketches['__profile__'].Line(point1=rib_trail_low[0], point2=[orthovector_trail_up[-1][0],orthovector_trail_low[0][1]])
	#WingModel.sketches['__profile__'].Line(point1=[orthovector_trail_up[-1][0]+0.004,orthovector_trail_up[-1][1]], point2=[orthovector_trail_up[-1][0],orthovector_trail_low[0][1]])
	WingModel.sketches['__profile__'].Line(point1=rib_trail_low[0], point2=[rib_trail_low[0][0]+airfoilDIST,rib_trail_low[0][1]+airfoilTH/fairfoilTH])
	WingModel.sketches['__profile__'].Line(point1=[rib_trail_up[-1][0]+airfoilDIST,rib_trail_up[-1][1]-airfoilTH/fairfoilTH], point2=[rib_trail_low[0][0]+airfoilDIST,rib_trail_low[0][1]+airfoilTH/fairfoilTH])
else:
	#WingModel.sketches['__profile__'].Line(point1=rib_trail_low[0], point2=[orthovector_trail_up[-1][0]+x_spar2-x_spar2m,orthovector_trail_low[0][1]])
	#WingModel.sketches['__profile__'].Line(point1=[orthovector_trail_up[-1][0]+0.004,orthovector_trail_up[-1][1]], point2=[orthovector_trail_up[-1][0]+x_spar2-x_spar2m,orthovector_trail_low[0][1]])
	WingModel.sketches['__profile__'].Line(point1=rib_trail_low[0], point2=[rib_trail_low[0][0]+airfoilDIST,rib_trail_low[0][1]+airfoilTH/fairfoilTH])
	WingModel.sketches['__profile__'].Line(point1=[rib_trail_up[-1][0]+airfoilDIST,rib_trail_up[-1][1]-airfoilTH/fairfoilTH], point2=[rib_trail_low[0][0]+airfoilDIST,rib_trail_low[0][1]+airfoilTH/fairfoilTH])
	
WingModel.Part(dimensionality=THREE_D, name='Part-Trail', type=
    DEFORMABLE_BODY)
WingModel.parts['Part-Trail'].BaseShell(sketch=
    WingModel.sketches['__profile__'])
del WingModel.sketches['__profile__']

##create the uncutted trail part
LeftCutTrail = 4#3#17
RightCutTrail = 5#4#20
##define the position of the cut poins
#HoeheCutGerade = -0.01 # Koordinate
HoeheCutGerade = rib_trail_low[LeftCutTrail][1]+3*0.005* (c_0/0.5)#0.01 # Koordinate ### Alejandro
HoeheCutGesamt = rib_trail_low[LeftCutTrail][1]+3*0.0125* (c_0/0.5)#0.02 #0.0 # Koordinate ### Alejandro
BeginCutTrail=rib_trail_low[LeftCutTrail]
EndCutTrail=rib_trail_low[RightCutTrail]

WingModel.ConstrainedSketch(gridSpacing=0.006, name='__profile__', sheetSize=c, transform=WingModel.parts['Part-Trail'].MakeSketchTransform(
    sketchPlane=WingModel.parts['Part-Trail'].faces.findAt((rib_trail_up[0][0], (rib_trail_up[0][1]+rib_trail_low[-1][1])/2, 0.0), (0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, 
    sketchUpEdge=WingModel.parts['Part-Trail'].edges.findAt((
    rib_trail_up[0][0], (rib_trail_up[0][1]+rib_trail_low[-1][1])/2, 0.0), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))


WingModel.sketches['__profile__'].Line(point1=(BeginCutTrail[0], rib_trail_low[0][1]-HoeheCutGesamt), point2=(BeginCutTrail[0], HoeheCutGerade))
WingModel.sketches['__profile__'].Line(point1=(EndCutTrail[0], rib_trail_low[0][1]-HoeheCutGesamt), point2=(EndCutTrail[0], HoeheCutGerade))
WingModel.sketches['__profile__'].Line(point1=(EndCutTrail[0], rib_trail_low[0][1]-HoeheCutGesamt), point2=(BeginCutTrail[0], rib_trail_low[0][1]-HoeheCutGesamt))
WingModel.sketches['__profile__'].Arc3Points(point1=(BeginCutTrail[0], HoeheCutGerade), point2=(EndCutTrail[0], HoeheCutGerade), point3=((EndCutTrail[0]+BeginCutTrail[0])/2, HoeheCutGesamt))

WingModel.parts['Part-Trail'].CutExtrude(
    flipExtrudeDirection=OFF, sketch=
    WingModel.sketches['__profile__'], sketchOrientation=RIGHT,
    sketchPlane=WingModel.parts['Part-Trail'].faces.findAt((
    rib_trail_up[0][0], (rib_trail_up[0][1]+rib_trail_low[-1][1])/2, 0.0), ), sketchPlaneSide=SIDE1, 
    sketchUpEdge=WingModel.parts['Part-Trail'].edges.findAt((
    rib_trail_up[0][0], (rib_trail_up[0][1]+rib_trail_low[-1][1])/2, 0.0),))

del WingModel.sketches['__profile__']

#################################################################################################################################### new rib nose ##############################################################################
## Rib Nose
#LeftCut=62    # gleiche Zahl heisst ein Punkt abstand
#RightCut=62
LeftCutNose = len(rib_nose) -5#-4#-16
RightCutNose = len(rib_nose) -4#-3#-13
#HoeheCutGerade=-0.01 # Koordinate
#HoeheCutGesamt=0.0 # Koordinate
HoeheCutGerade = rib_nose[RightCutNose][1]+3*0.005* (c_0/0.5) # Koordinate ### Alejandro
HoeheCutGesamt = rib_nose[RightCutNose][1]+3*0.0125* (c_0/0.5) #0.015 #0.0 # Koordinate ### Alejandro
##parameters specifying the position of the cuts and the height of the cut


#Begin_Cut=points_nose[0:LeftCut]
#End_Cut=points_nose[RightCut:len(points_nose)]
BeginCutNose=rib_nose[LeftCutNose]
EndCutNose=rib_nose[RightCutNose]

WingModel.ConstrainedSketch(name='__profile__', sheetSize=2*c)
WingModel.sketches['__profile__'].Spline(points=rib_nose[:])
WingModel.sketches['__profile__'].Line(point1=(rib_nose[0][0], rib_nose[0][1]), point2=(rib_nose[-1][0], rib_nose[-1][1]))
##create the sketch of the nose
WingModel.Part(dimensionality=THREE_D, name='Rib-Nose-New', type=DEFORMABLE_BODY)
WingModel.parts['Rib-Nose-New'].BaseShell(sketch=WingModel.sketches['__profile__'])
##create the uncutted nose part
del WingModel.sketches['__profile__']
#WingModel.ConstrainedSketch(gridSpacing=0.006, name='__profile__', sheetSize=c, transform=WingModel.parts['Rib-Nose-New'].MakeSketchTransform(
#    sketchPlane=WingModel.parts['Rib-Nose-New'].faces.findAt((rib_nose[0][0], 0.0, 0.0), (0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, sketchUpEdge=WingModel.parts['Rib-Nose-New'].edges.findAt((
#    rib_nose[0][0], 0.0, 0.0), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
WingModel.ConstrainedSketch(gridSpacing=0.006, name='__profile__', sheetSize=c, transform=WingModel.parts['Rib-Nose-New'].MakeSketchTransform(
    sketchPlane=WingModel.parts['Rib-Nose-New'].faces.findAt((rib_nose[0][0], rib_nose[0][1], 0.0), (0.0, 0.0, 1.0)), sketchPlaneSide=SIDE1, sketchUpEdge=WingModel.parts['Rib-Nose-New'].edges.findAt((
#    rib_nose[0][0], (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    (rib_nose[0][0]+rib_nose[-1][0])/2, (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
##0.110603, 0.009319, 0.0



WingModel.sketches['__profile__'].Line(point1=(BeginCutNose[0], rib_nose[-1][1]-HoeheCutGesamt), point2=(BeginCutNose[0], HoeheCutGerade))
WingModel.sketches['__profile__'].Line(point1=(EndCutNose[0], rib_nose[-1][1]-HoeheCutGesamt), point2=(EndCutNose[0], HoeheCutGerade))
WingModel.sketches['__profile__'].Arc3Points(point1=(BeginCutNose[0], HoeheCutGerade), point2=(EndCutNose[0], HoeheCutGerade), point3=((EndCutNose[0]+BeginCutNose[0])/2, HoeheCutGesamt))
WingModel.sketches['__profile__'].Line(point1=(EndCutNose[0], rib_nose[-1][1]-HoeheCutGesamt), point2=(BeginCutNose[0], rib_nose[-1][1]-HoeheCutGesamt))

WingModel.parts['Rib-Nose-New'].CutExtrude(
    flipExtrudeDirection=OFF, sketch=
    WingModel.sketches['__profile__'], sketchOrientation=RIGHT,
    sketchPlane=WingModel.parts['Rib-Nose-New'].faces.findAt((
    points_nose[0][0], rib_nose[0][1], 0.0), ), sketchPlaneSide=SIDE1, 
    sketchUpEdge=WingModel.parts['Rib-Nose-New'].edges.findAt((
#    points_nose[0][0], (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0),))
    (rib_nose[0][0]+rib_nose[-1][0])/2, (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0),))
##create the cutted nose part
del WingModel.sketches['__profile__']

######################

#===============================================================================
# 4-PartDefinition
#===============================================================================

# part definition -all 

WingModel.ConstrainedSketch(name='__profile__', sheetSize=2*c)     
WingModel.sketches['__profile__'].Spline(points= (rib_trail_up[:]  ))     
WingModel.sketches['__profile__'].Spline(points= (rib_box_up[:]     ))
WingModel.sketches['__profile__'].Spline(points= (points_nose[:]    ))     
WingModel.sketches['__profile__'].Spline(points= (rib_box_low[:]    ))
WingModel.sketches['__profile__'].Spline(points= (rib_trail_low[:]  ))
WingModel.sketches['__profile__'].Line(point1 = rib_trail_low[-1], point2 = rib_trail_up[0]  )
          
# first spar

WingModel.sketches['__profile__'].Line(point1=points_nose[0], point2=points_nose[-1])

# 2nd spar

WingModel.sketches['__profile__'].Line(point1=points_box_up[0], point2=points_box_low[-1])

WingModel.Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
WingModel.parts['Part-1'].BaseShellExtrude(depth=width, sketch=
    WingModel.sketches['__profile__'])
del WingModel.sketches['__profile__']

### Tasche einfuegen - Small gaps close to the rib, WHY??

## Trail
frontplane=WingModel.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=XYPLANE)
##DatumPlaneByPrincipalPlane(), this method creates a Feature object and a Datumplane object through the origin along one of the principal planes
WingModel.ConstrainedSketch(gridSpacing=0.06, name='__profile__', 
    sheetSize=2*c, transform=
    WingModel.parts['Part-1'].MakeSketchTransform(
    sketchPlane=WingModel.parts['Part-1'].datums[frontplane.id], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=WingModel.parts['Part-1'].edges.findAt((rib_trail_up[0][0], (rib_trail_up[0][1]+rib_trail_low[-1][1])/2, 0.0), ), 
    sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
WingModel.parts['Part-1'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=WingModel.sketches['__profile__'])
##this method projects the vertices of specified edges, and datum points from the part onto the specified ConstrainedSketch object. The vertices and datum points appear on the sketch as reference geometry.
#WingModel.sketches['__profile__'].rectangle(point1=(
#    BeginCutTrail[0], 0), point2=(EndCutTrail[0], -1))
WingModel.sketches['__profile__'].rectangle(point1=(
    BeginCutTrail[0], BeginCutTrail[1]+(0.005* (c_0/0.5))), point2=(EndCutTrail[0], -1)) #Modified by Alejandro Valverde
WingModel.parts['Part-1'].CutExtrude(flipExtrudeDirection=ON, 
    sketch=WingModel.sketches['__profile__'], sketchOrientation=
    RIGHT, sketchPlane=WingModel.parts['Part-1'].datums[2], 
    sketchPlaneSide=SIDE1, sketchUpEdge=
    WingModel.parts['Part-1'].edges.findAt((rib_trail_up[0][0], (rib_trail_up[0][1]+rib_trail_low[-1][1])/2, 0.0), 
    ))

## Nose   
frontplane=WingModel.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=XYPLANE)
WingModel.ConstrainedSketch(gridSpacing=0.06, name='__profile__', 
    sheetSize=2*c, transform=
    WingModel.parts['Part-1'].MakeSketchTransform(
    sketchPlane=WingModel.parts['Part-1'].datums[frontplane.id], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=WingModel.parts['Part-1'].edges.findAt((
#    points_nose[0][0], (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    (rib_nose[0][0]+rib_nose[-1][0])/2, (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0), ), sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
WingModel.parts['Part-1'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=WingModel.sketches['__profile__'])
WingModel.sketches['__profile__'].rectangle(point1=(
    BeginCutNose[0], BeginCutNose[1]+(0.005* (c_0/0.5))), point2=(EndCutNose[0], -1)) #Modified by Alejandro Valverde
WingModel.parts['Part-1'].CutExtrude(flipExtrudeDirection=ON, 
    sketch=WingModel.sketches['__profile__'], sketchOrientation=
    RIGHT, sketchPlane=WingModel.parts['Part-1'].datums[2], 
    sketchPlaneSide=SIDE1, sketchUpEdge=
#    WingModel.parts['Part-1'].edges.findAt((points_nose[0][0], (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0), 
    WingModel.parts['Part-1'].edges.findAt(((rib_nose[0][0]+rib_nose[-1][0])/2, (rib_nose[0][1]+rib_nose[-1][1])/2, 0.0), 
    ))

#===============================================================================
# 5-Materials definition
#===============================================================================

winkelGlass=0

WingModel.Material(name='Aluminium_wingSkin')
WingModel.materials['Aluminium_wingSkin'].Density(table=((2.7e-6, ), )) #Modified, Alejandro Valverde, units in mm
WingModel.materials['Aluminium_wingSkin'].Elastic(table=((70000.0, 0.3), )) #Modified, Alejandro Valverde, units in mm


## http://www.gordoncomposites.com/products/TDS/GC-70-UL.pdf
#WingModel.Material(name='GFRP-UD-SparRear')
#WingModel.materials['GFRP-UD-SparRear'].Elastic(type=LAMINA, table=(
#    (41368543759.0, 10342135940.0, 0.29, 4481592241.0, 4481592241.0, 3378431074.0), ))
#WingModel.materials['GFRP-UD-SparRear'].elastic.FailStress(table=(
#    (1048003109.0, 765318060.0, 466884350.0, 146168855.0, 51021204.0, 0.0, ), ))
#WingModel.materials['GFRP-UD-SparRear'].elastic.FailStrain(table=((0.025, 0.019, ), ))

WingModel.rootAssembly.regenerate()

# WingModel.Material(name='CFRP-UD-Box-M40')
# WingModel.Material(name='CFRP-UD-Box-G')
# WingModel.Material(name='CFRP-UD-Box-G-Hoff')
# WingModel.Material(name='Schwalbe')

# https://www.acpsales.com/upload/Mechanical-Properties-of-Carbon-Fiber-Composite-Materials.pdf

# WingModel.Material(name='CFRP-UD-Box-M55')
# WingModel.materials['CFRP-UD-Box-M55'].Elastic(table=((300000000000.0, 12000000000.0, 0.3, 5000000000.0, 5000000000.0, 5000000000.0), ), type=LAMINA)
# WingModel.materials['CFRP-UD-Box-M55'].elastic.FailStress(table=((1600000000.0, 1300000000.0, 50000000.0, 250000000.0, 75000000.0, 0.0, ), ))

### Giulios Daten

WingModel.Material(name='CFRP-UD-Box')
WingModel.materials['CFRP-UD-Box'].Elastic(table=((135000000000.0, 10000000000.0, 0.3, 5000000000.0, 5000000000.0, 3846000000.0), ), type=LAMINA)
WingModel.materials['CFRP-UD-Box'].elastic.FailStress(table=((1500000000.0, 1200000000.0, 55000000.0, 250000000.0, 70000000.0, 0.0, ), ))
WingModel.materials['CFRP-UD-Box'].elastic.FailStrain(table=(( 0.0105, 0.0085, 0.005, 0.025, 0.014), ))

# WingModel.Material(name='Glass-Fabric')
# WingModel.materials['Glass-Fabric'].Elastic(table=((25000000000.0, 25000000000.0, 0.2, 4000000000.0, 4000000000.0, 4000000000.0), ), type=LAMINA)
# WingModel.materials['Glass-Fabric'].elastic.FailStress(table=((440000000.0, 425000000.0, 440000000.0, 425000000.0, 40000000.0, 0.0, ), ))
# WingModel.materials['Glass-Fabric'].elastic.FailStrain(table=((0.0175, 0.017, 0.0175, 0.017, 1.0), ))

# WingModel.Material(name='CFK-Fabric')
# WingModel.materials['CFK-Fabric'].Elastic(table=((70000000000.0, 70000000000.0, 0.1, 5000000000.0, 5000000000.0, 5000000000.0), ), type=LAMINA)
# WingModel.materials['CFK-Fabric'].elastic.FailStress(table=((600000000.0, 570000000.0, 600000000.0, 570000000.0, 90000000.0, 0.0, ), ))
# WingModel.materials['CFK-Fabric'].elastic.FailStrain(table=((0.0085, 0.008, 0.0085, 0.008, 0.018), ))

#===============================================================================
# 6-FacesAndSets
#===============================================================================

faces_skin2a=WingModel.parts['Part-1'].faces.findAt(( ( points_box_up[1][0],points_box_up[1][1],0.0   ) ,),)
faces_skin2b=WingModel.parts['Part-1'].faces.findAt(( ( points_box_low[1][0],points_box_low[1][1],0.0  ) ,),)

faces_skin1=WingModel.parts['Part-1'].faces.findAt(   ( (points_nose[1][0],points_nose[1][1],0.0 ) ,),                                                                
                                                                  ( (points_nose[-2][0],points_nose[-2][1],0.0 ) ,),)


faces_skin3=WingModel.parts['Part-1'].faces.findAt( ( ( points_trail1_low[1][0],points_trail1_low[1][1],0.0 ) ,),
                                                               ( ( points_trail3_low[1][0],points_trail3_low[1][1],0.0  ) ,),
                                                               ( ( points_trail3_up[1][0],points_trail3_up[1][1],0.0 ) ,),
                                                               ( ( points_trail4_low[1][0],points_trail4_low[1][1],0.0  ) ,),
                                                               ( ( points_trail4_up[1][0],points_trail4_up[1][1],0.0 ) ,),
                                                                  ( ( points_trail2_low[1][0],points_trail2_low[1][1],0.0  ) ,),
                                                               ( ( points_trail1_up[1][0],points_trail1_up[1][1],0.0 ) ,),
                                                               ( ( points_trail2_up[1][0],points_trail2_up[1][1],0.0  ) ,),)
                                                               
faces_skin4 = WingModel.parts['Part-1'].faces.getByBoundingBox(0.99*c, -0.1*c, -0.01*width, 1.01*c, 0.1*c, 1.01*width)

##setskin=WingModel.parts['Part-1'].Set(faces=faces_skin , name='skin_set')
WingModel.parts['Part-1'].Surface(name='Surf-Skin1', side2Faces=faces_skin1)
WingModel.parts['Part-1'].Surface(name='Surf-Skin2a', side1Faces=faces_skin2a)
WingModel.parts['Part-1'].Surface(name='Surf-Skin2b', side1Faces=faces_skin2b)
WingModel.parts['Part-1'].Surface(name='Surf-Skin3', side2Faces=faces_skin3)
WingModel.parts['Part-1'].Surface(name='Surf-Skin4', side2Faces=faces_skin4)

WingModel.parts['Part-1'].SurfaceByBoolean(name='Surf-Skin', surfaces=
    (WingModel.parts['Part-1'].surfaces['Surf-Skin1'], WingModel.parts['Part-1'].surfaces['Surf-Skin3'], 
    WingModel.parts['Part-1'].surfaces['Surf-Skin2a'], WingModel.parts['Part-1'].surfaces['Surf-Skin2b'],
    WingModel.parts['Part-1'].surfaces['Surf-Skin4']))

faces_box=WingModel.parts['Part-1'].faces.findAt( ( ( points_box_up[1][0],points_box_up[1][1],0.0  ) ,),
                                                               ( ( points_box_low[1][0],points_box_low[1][1],0.0 ) ,),( ( points_box_low[-2][0],points_box_low[-2][1],0.0  ) ,),)

faces_box_oben=WingModel.parts['Part-1'].faces.findAt( ( ( points_box_up[1][0],points_box_up[1][1],0.0  ) ,),)

faces_box_unten=WingModel.parts['Part-1'].faces.findAt( ( ( points_box_low[1][0],points_box_low[1][1],0.0  ) ,),( ( points_box_low[-2][0],points_box_low[-2][1],0.0  ) ,),)


faces_front=WingModel.parts['Part-1'].faces.findAt(
                                                                ( ( points_nose[1][0], points_nose[1][1], 0.0 ) ,),
                                                                ( ( points_nose[-2][0], points_nose[-2][1], 0.0 ) ,),)


faces_rear_unten1=WingModel.parts['Part-1'].faces.findAt( ( ( points_trail1_low[1][0], points_trail1_low[1][1], 0.0 ) ,))

faces_rear_unten2=WingModel.parts['Part-1'].faces.findAt( ( ( points_trail4_low[1][0], points_trail4_low[1][1], 0.0 ) ,),)

faces_rear_oben=WingModel.parts['Part-1'].faces.findAt( ( ( points_trail1_up[1][0],points_trail1_up[1][1],0.0 ) ,),
                                                               ( ( points_trail3_up[1][0],points_trail3_up[1][1],0.0 ) ,),
                                                               ( ( points_trail4_up[1][0],points_trail4_up[1][1],0.0 ) ,),
                                                               ( ( points_trail2_up[1][0],points_trail2_up[1][1],0.0) ,),)
                                                               
faces_rear_end = WingModel.parts['Part-1'].faces.getByBoundingBox(0.99*c, -0.1*c, -0.01*width, 1.01*c, 0.1*c, 1.01*width)                                                               

faces_frontrear_G=WingModel.parts['Part-1'].faces.findAt( ( ( points_trail3_low[1][0],points_trail3_low[1][1],0.0  ) ,),)
faces_corr_front=WingModel.parts['Part-1'].faces.findAt( ( ( points_nose[-2][0], points_nose[-2][1], 0.0 ) ,),)

##faces_rear_end=WingModel.parts['Part-1'].faces.findAt( ( ( c,0,width ) ,),)

#setspar1=WingModel.parts['Part-1'].faces.findAt(( ( points_nose[0][0]   ,0.0,0.001 ) , )  , )
setspar1=WingModel.parts['Part-1'].faces.findAt(( ( points_nose[0][0]   ,points_nose[0][1] ,0.001 ) , )  , )

setspar1_Up=WingModel.parts['Part-1'].faces.findAt(( ( points_nose[0][0]   ,points_nose[0][1]-0.0025 ,0.001 ) , )  , )
#setspar1_mid=WingModel.parts['Part-1'].faces.findAt(( ( points_nose[0][0]   ,0.0,0.001 ) , )  , )
setspar1_mid=WingModel.parts['Part-1'].faces.findAt(( ( points_nose[0][0]   ,points_nose[0][1],0.001 ) , )  , )
setspar1_Low=WingModel.parts['Part-1'].faces.findAt(( ( points_nose[0][0]   ,points_nose[-1][1]+0.0025,0.001 ) , )  , )

#setspar2=WingModel.parts['Part-1'].faces.findAt(( ( points_box_up[0][0]   ,0.0,0.001 ) , )  , )
setspar2=WingModel.parts['Part-1'].faces.findAt(( ( points_box_up[0][0]   ,points_box_up[0][1] ,0.001 ) , )  , )

setbox=WingModel.parts['Part-1'].Set(faces=faces_box , name='box_set')

setbox_oben=WingModel.parts['Part-1'].Set(faces=faces_box_oben , name='box_set_oben')
setbox_unten=WingModel.parts['Part-1'].Set(faces=faces_box_unten , name='box_set_unten')

setfront=WingModel.parts['Part-1'].Set(faces=faces_front , name='front_set')
setrearoben=WingModel.parts['Part-1'].Set(faces=faces_rear_oben , name='rear_set_oben')
setrearunten1=WingModel.parts['Part-1'].Set(faces=faces_rear_unten1 , name='rear_set_unten1')
setrearunten2=WingModel.parts['Part-1'].Set(faces=faces_rear_unten2 , name='rear_set_unten2')


##setfrontrear_G=WingModel.parts['Part-1'].Set(faces=faces_frontrear_G , name='frontrear_set_G')
##set_corr_front=WingModel.parts['Part-1'].Set(faces=faces_corr_front , name='corr_front')

setrear_end=WingModel.parts['Part-1'].Set(faces=faces_rear_end , name='rear_end')

setfull_Spar1=WingModel.parts['Part-1'].Set(faces=setspar1 , name='fullpart_set_Spar1')

setfull_Spar1_Up=WingModel.parts['Part-1'].Set(faces=setspar1_Up , name='fullpart_set_Spar1_Up')
setfull_Spar1_Mid=WingModel.parts['Part-1'].Set(faces=setspar1_mid , name='fullpart_set_Spar1_Mid')
setfull_Spar1_Low=WingModel.parts['Part-1'].Set(faces=setspar1_Low , name='fullpart_set_Spar1_Low')

setfull_Spar2=WingModel.parts['Part-1'].Set(faces=setspar2 , name='fullpart_set_Spar2')

meshSection=WingModel.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=points_box_up[-4][0], 
    principalPlane=YZPLANE)
meshSection1=WingModel.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=points_box_up[-1][0], 
    principalPlane=YZPLANE)

#===============================================================================
# 7-FacesAndSets
#===============================================================================

# Sets for airfoil sections at control points
#for i in range(N):
#    b_plane=WingModel.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=y_C[N+i], principalPlane=XYPLANE)
#    WingModel.parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=WingModel.parts['Part-1'].datums[b_plane.id], faces=WingModel.parts['Part-1'].faces.getByBoundingBox(-0.01*c, -c, -0.1*c, 1.01*c, c, width+0.1*c))
#    edges_section_control1 = WingModel.parts['Part-1'].edges.getByBoundingBox(-0.01*c, -0.001, y_C[N+i]-0.1*element_width, 1.01*c, c, y_C[N+i]+0.1*element_width)   
#    edges_section_control2 = WingModel.parts['Part-1'].edges.getByBoundingBox(-0.01*c, -c, y_C[N+i]-0.1*element_width, 1.01*c, 0.001, y_C[N+i]+0.1*element_width)   
#    edges_section_control3 = WingModel.parts['Part-1'].edges.findAt(((points_nose[1][0], points_nose[1][1], y_C[N+i]) ,),)
#    WingModel.parts['Part-1'].Set(name='Section'+str(i+1) , edges=edges_section_control1+edges_section_control2+edges_section_control3)
for i in range(N):
    b_plane=WingModel.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=y_C[N+i], principalPlane=XYPLANE)
    WingModel.parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=WingModel.parts['Part-1'].datums[b_plane.id], faces=WingModel.parts['Part-1'].faces.getByBoundingBox(-0.01*c, -c, -0.1*c, 1.01*c, c, width+0.1*c))
    edges_section_control1 = WingModel.parts['Part-1'].edges.findAt(((mod6_UP_r[0], mod6_UP_r[1], y_C[N+i]) ,),)
    edges_section_control2 = WingModel.parts['Part-1'].edges.findAt(((mod6_UP_m[0], mod6_UP_m[1], y_C[N+i]) ,),)
    edges_section_control3 = WingModel.parts['Part-1'].edges.findAt(((mod6_UP_f[0], mod6_UP_f[1], y_C[N+i]) ,),)
    edges_section_control4 = WingModel.parts['Part-1'].edges.findAt(((mod6_DN_m[0], mod6_DN_m[1], y_C[N+i]) ,),)
    edges_section_control5 = WingModel.parts['Part-1'].edges.findAt(((rib_nose[LeftCutNose-1][0], rib_nose[LeftCutNose-1][1], y_C[N+i]) ,),)
    edges_section_control6 = WingModel.parts['Part-1'].edges.findAt(((rib_nose[RightCutNose+1][0], rib_nose[RightCutNose+1][1], y_C[N+i]) ,),)
    edges_section_control7 = WingModel.parts['Part-1'].edges.findAt(((rib_trail_low[LeftCutTrail-1][0], rib_trail_low[LeftCutTrail-1][1], y_C[N+i]) ,),)
    edges_section_control8 = WingModel.parts['Part-1'].edges.findAt(((rib_trail_low[RightCutTrail+1][0], rib_trail_low[RightCutTrail+1][1], y_C[N+i]) ,),)
    edges_section_control9 = WingModel.parts['Part-1'].edges.findAt(((rib_trail_up[0][0], (rib_trail_up[0][1]+rib_trail_low[-1][1])/2, y_C[N+i]) ,),)
    edges_section_control10 = WingModel.parts['Part-1'].edges.findAt(((points_nose[1][0], points_nose[1][1], y_C[N+i]) ,),)
    WingModel.parts['Part-1'].Set(name='Section'+str(i+1) , edges=edges_section_control1+edges_section_control2+edges_section_control3+edges_section_control4+edges_section_control5+edges_section_control6+edges_section_control7+edges_section_control8+edges_section_control9+edges_section_control10)	


#x_SC = 0.125
#a_plane=WingModel.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=x_SC, principalPlane=YZPLANE)
#WingModel.parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=WingModel.parts['Part-1'].datums[a_plane.id], 
#    faces=WingModel.parts['Part-1'].faces.getByBoundingBox(-0.01*c, -c, y_C[N]-1E-4, 1.01*c, c, 1.01*width))
    
anzDivision=round(Anz_Unterteilung/anzRibs)

WingModel.rootAssembly.Instance(dependent=OFF, name='Part-1-1', part=WingModel.parts['Part-1'])


for i in xrange(1,anzRibs+1):
    WingModel.rootAssembly.Instance(dependent=OFF, name='Part-rib_nose-'+str(i), part=WingModel.parts['Rib-Nose-New'])
    WingModel.rootAssembly.Instance(dependent=OFF, name='Part-rib_trail_low-'+str(i), part=WingModel.parts['Part-Trail'])
    #WingModel.rootAssembly.Instance(dependent=OFF, name='Part-rib_box-'+str(i), part=WingModel.parts['Part-Box-'+str(BucklingSpar)])
    WingModel.rootAssembly.translate(instanceList=('Part-rib_nose-'+str(i), ) , vector=(0.0, 0.0, width-(i-1)*anzDivision*width/Anz_Unterteilung))
    WingModel.rootAssembly.translate(instanceList=('Part-rib_trail_low-'+str(i), ) , vector=(0.0, 0.0, width-(i-1)*anzDivision*width/Anz_Unterteilung))
    #WingModel.rootAssembly.translate(instanceList=('Part-rib_box-'+str(i), ) , vector=(0.0, 0.0, width-(i-1)*anzDivision*width/Anz_Unterteilung))

#Manual addition, Alejandro Valverde
WingModel.rootAssembly.Instance(dependent=OFF, name='Part-rib_nose-'+str(i+1), part=WingModel.parts['Rib-Nose-New'])
WingModel.rootAssembly.Instance(dependent=OFF, name='Part-rib_trail_low-'+str(i+1), part=WingModel.parts['Part-Trail'])
#WingModel.rootAssembly.Instance(dependent=OFF, name='Part-rib_box-'+str(i+1), part=WingModel.parts['Part-Box-'+str(BucklingSpar)])
# WingModel.rootAssembly.translate(instanceList=('Part-rib_nose-'+str(i+1), ) , vector=(0.0, 0.0, width-(i-1)*anzDivision*width/Anz_Unterteilung)) #This has been changed, Alejandro Valverde
# WingModel.rootAssembly.translate(instanceList=('Part-rib_trail_low-'+str(i+1), ) , vector=(0.0, 0.0, width-(i-1)*anzDivision*width/Anz_Unterteilung))

ribs1=[]
ribs2=[]
ribs3=[]
ribs4=[]
ribs5=[]
ribs6=[]
ribs7=[]
ribs8=[]
ribs9=[]

for i in xrange(1,anzRibs+1):
    ribs1.append(WingModel.rootAssembly.instances['Part-rib_trail_low-'+str(i)])
    #ribs8.append(WingModel.rootAssembly.instances['Part-rib_box-'+str(i)])
    ribs9.append(WingModel.rootAssembly.instances['Part-rib_nose-'+str(i)])

#Manual addition, Alejandro Valverde
ribs1.append(WingModel.rootAssembly.instances['Part-rib_trail_low-'+str(i+1)])
ribs9.append(WingModel.rootAssembly.instances['Part-rib_nose-'+str(i+1)])

ribs1.append(WingModel.rootAssembly.instances['Part-1-1'])

ribs1=tuple(ribs1)
ribs8=tuple(ribs8)
ribs9=tuple(ribs9)

rib=ribs1+ribs8+ribs9

    
WingModel.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(         rib         ), name='Wing',
    originalInstances=SUPPRESS)

for i in range(N):
    if i == N-1 or i == 0 or i == 12 or i == 15 or i == 19:
        edges_section_SC = WingModel.parts['Wing'].edges.getByBoundingBox(-0.01*c, -c, y_C[N+i]-0.1*element_width, 1.01*c, c, y_C[N+i]+0.1*element_width)
        WingModel.parts['Wing'].Set(name='Section'+str(i+1)+'_SC_edges', edges=edges_section_SC)
    else:
        pass

## hier werden die Kanten vorne und hinten erzeugt, um spaeter den AOA auszulesen... uebernommen von Mesh.py - Obtain edges that would be used during the meshing operations
front_nose_plane=WingModel.parts['Wing'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)
front_nose_planeOFF=WingModel.parts['Wing'].DatumPlaneByPrincipalPlane(offset=(mod7_UP[1]+mod7_DN[1])/2, principalPlane=XZPLANE)
faces_nodalpoints=WingModel.parts['Wing'].faces.findAt( ( ( points_nose[25][0], points_nose[25][1], 0.001), ),)

#for i in xrange(1,anzRibs+1): #anzRibs+1
#    faces_nodalpoints=faces_nodalpoints+WingModel.parts['Wing'].faces.findAt( ( ( points_trail3_up[0][0], points_trail4_up[0][1]-0.00001, width-(i-1)*anzDivision*width/Anz_Unterteilung) ,),)
#
#for i in xrange(1,anzRibs+1): #anzRibs+1
#    faces_nodalpoints=faces_nodalpoints+WingModel.parts['Wing'].faces.findAt(((
#    points_box_up[0][0], points_box_up[0][1]-0.00001, width-(i-1)*anzDivision*width/Anz_Unterteilung), ),)+ \
#    WingModel.parts['Wing'].faces.findAt(((points_box_up[1][0], points_box_up[1][1]-0.00001, width-(i-1)*anzDivision*width/Anz_Unterteilung), ),)
#                                         
for i in range(Anz_Unterteilung): #anzRibs+1
    faces_nodalpoints=faces_nodalpoints+WingModel.parts['Wing'].faces.findAt(((
    points_nose[25][0], points_nose[25][1], y_C[N+i]+0.001), ),)+ \
    WingModel.parts['Wing'].faces.findAt(((points_nose[25][0], points_nose[25][1], y_C[N+i]-0.001), ),)
#
#faces_nodalpoints = faces_nodalpoints+WingModel.parts['Wing'].faces.getByBoundingBox(
#    points_trail4_low[0][0], points_trail4_low[0][1], -0.01*width, 1.01*c, points_trail4_up[-1][1], 1.01*width)
faces_nodalpoints = faces_nodalpoints+WingModel.parts['Wing'].faces.getByBoundingBox(
    points_box_up[-1][0], -c, y_C[N]-1E-4, 1.01*c, c, 1.01*width)
    
#WingModel.parts['Wing'].PartitionFaceByDatumPlane(datumPlane=
#    WingModel.parts['Wing'].datums[front_nose_plane.id], faces=faces_nodalpoints)
#WingModel.parts['Wing'].PartitionFaceByDatumPlane(datumPlane=
#    WingModel.parts['Wing'].datums[front_nose_plane.id], faces=WingModel.parts['Wing'].faces.getByBoundingBox(
#    -0.01, -c, -0.01, x_spar1+0.01, c, 1.01*width) + WingModel.parts['Wing'].faces.getByBoundingBox(
#    x_spar2+0.01, -c, -0.01, c+0.01, c, 1.01*width))
WingModel.parts['Wing'].PartitionFaceByDatumPlane(datumPlane=
    WingModel.parts['Wing'].datums[front_nose_planeOFF.id], faces=WingModel.parts['Wing'].faces.getByBoundingBox(
    x_spar1-0.01, -c, -0.01, x_spar2+0.01, c, 1.01*width))

#partition face by sketch
dP = WingModel.parts['Wing'].DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
#dP = WingModel.parts['Wing'].datum.keys()[-1]
dA = WingModel.parts['Wing'].DatumAxisByPrincipalAxis(principalAxis=YAXIS)
#dA = WingModel.parts['Wing'].datum.keys()[-1]

s = WingModel.ConstrainedSketch(gridSpacing=0.006, name='__profile__', sheetSize=c, transform=WingModel.parts['Wing'].MakeSketchTransform(
    sketchPlane=WingModel.parts['Wing'].datums[dP.id], sketchPlaneSide=SIDE2, 
    sketchUpEdge=WingModel.parts['Wing'].datums[dA.id], 
    origin=(0.0, 0.0, 0.0)))
s.Line(point1=(-0.01, 0.0), point2=(0.01, 0.0))	
WingModel.parts['Wing'].PartitionFaceBySketchThruAll(sketchPlane=WingModel.parts['Wing'].datums[dP.id], 
    sketchUpEdge=WingModel.parts['Wing'].datums[dA.id], 
    faces=WingModel.parts['Wing'].faces.getByBoundingBox(-0.01, -c, -0.01, c+0.01, c, 1.01*width), 
    sketchPlaneSide=SIDE2, sketch=s)
del WingModel.sketches['__profile__']

s = WingModel.ConstrainedSketch(gridSpacing=0.006, name='__profile__', sheetSize=c, transform=WingModel.parts['Wing'].MakeSketchTransform(
    sketchPlane=WingModel.parts['Wing'].datums[dP.id], sketchPlaneSide=SIDE2, 
    sketchUpEdge=WingModel.parts['Wing'].datums[dA.id], 
    origin=(0.0, 0.0, 0.0)))
s.Line(point1=(-c+0.01, 0.0), point2=(-c-0.01, 0.0))	
#s.Line(point1=(c-0.1, 0.0), point2=(c+0.1, 0.0))	
WingModel.parts['Wing'].PartitionFaceBySketchThruAll(sketchPlane=WingModel.parts['Wing'].datums[dP.id], 
    sketchUpEdge=WingModel.parts['Wing'].datums[dA.id], 
    faces=WingModel.parts['Wing'].faces.getByBoundingBox(c-0.01, -c, -0.01, c+0.01, c, 1.01*width), 
    sketchPlaneSide=SIDE2, sketch=s)
del WingModel.sketches['__profile__']

	
for i in xrange(1, anzRibs):
    c_plane=WingModel.parts['Wing'].DatumPlaneByPrincipalPlane(offset=width-i*anzDivision*element_width, principalPlane=XYPLANE)
    WingModel.parts['Wing'].PartitionFaceByDatumPlane(datumPlane=WingModel.parts['Wing'].datums[c_plane.id], 
    faces=WingModel.parts['Wing'].faces.getByBoundingBox(-0.01*c, -c, -0.1*c, 1.01*c, c, width+0.1*c))

## Sets for cross sections at wing tip for shear center determination
#vertex1_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_nose[0][0], 0, width), 1E-05) # 0
#vertex2_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_box_up[0][0], 0, width), 1E-05) # points_trail3_up[0][0]
#vertex3_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_nose[0][0], points_nose[0][1], width), 1E-05)
#vertex4_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_box_up[0][0], points_box_up[0][1], width), 1E-05)
#WingModel.parts['Wing'].Set(name='MeasuringPoint'+str(anzRibs)+'_Front' , vertices=vertex1_section_SC)
#WingModel.parts['Wing'].Set(name='MeasuringPoint'+str(anzRibs)+'_Rear' , vertices=vertex2_section_SC)
#WingModel.parts['Wing'].Set(name='Rib'+str(anzRibs)+'_SC_vertices' , vertices=vertex1_section_SC+vertex2_section_SC)
#WingModel.parts['Wing'].Set(name='LoadingPoint'+str(anzRibs)+'_Front' , vertices=vertex3_section_SC)
#WingModel.parts['Wing'].Set(name='LoadingPoint'+str(anzRibs)+'_Rear' , vertices=vertex4_section_SC)
#
## Sets for cross sections at control points for reading out angles
#for i in range(N):
#    vertex1_section = WingModel.parts['Wing'].vertices.getByBoundingSphere((rib_box_up[-1][0], 0, y_C[N+i]), 1E-04)     # ((0, 0, y_C[N+i]), 1E-04)
#    vertex2_section = WingModel.parts['Wing'].vertices.getByBoundingSphere((rib_box_up[0][0], 0, y_C[N+i]), 1E-04)    # ((points_trail4_up[0][0], 0, y_C[N+i]), 1E-04)
#    WingModel.parts['Wing'].Set(name='MeasuringPointSection'+str(i+1)+'_Front' , vertices=vertex1_section)
#    WingModel.parts['Wing'].Set(name='MeasuringPointSection'+str(i+1)+'_Rear' , vertices=vertex2_section)
#    
## Sets for cross sections at ribs for reading out angles
#for i in range(anzRibs):
#    vertex1_rib = WingModel.parts['Wing'].vertices.getByBoundingSphere((rib_box_up[-1][0], 0, width-i*element_width*anzDivision), 1E-04)    # ((0, 0, width-i*element_width*anzDivision), 1E-04)
#    vertex2_rib = WingModel.parts['Wing'].vertices.getByBoundingSphere((rib_box_up[0][0], 0, width-i*element_width*anzDivision), 1E-04)   # ((points_trail4_up[0][0], 0, width-i*element_width*anzDivision), 1E-04)
#    WingModel.parts['Wing'].Set(name='MeasuringPointRib'+str(anzRibs-i)+'_Front' , vertices=vertex1_rib)
#    WingModel.parts['Wing'].Set(name='MeasuringPointRib'+str(anzRibs-i)+'_Rear' , vertices=vertex2_rib)
    

WingModel.parts['Wing'].sets['fullpart_set_Spar2']
regions = WingModel.parts['Wing'].sets['fullpart_set_Spar2']
WingModel.parts['Wing'].flipNormal(regions=regions)
WingModel.parts['Wing'].sets['front_set']
regions = WingModel.parts['Wing'].sets['front_set']
WingModel.parts['Wing'].flipNormal(regions=regions)
WingModel.parts['Wing'].sets['rear_end']
regions = WingModel.parts['Wing'].sets['rear_end']
WingModel.parts['Wing'].flipNormal(regions=regions)
WingModel.parts['Wing'].sets['rear_set_oben']
regions = WingModel.parts['Wing'].sets['rear_set_oben']
WingModel.parts['Wing'].flipNormal(regions=regions)
WingModel.parts['Wing'].sets['rear_set_unten1']
regions = WingModel.parts['Wing'].sets['rear_set_unten1']
WingModel.parts['Wing'].flipNormal(regions=regions)
WingModel.parts['Wing'].sets['rear_set_unten2']
regions = WingModel.parts['Wing'].sets['rear_set_unten2']
WingModel.parts['Wing'].flipNormal(regions=regions)

#===============================================================================
# 8-Geometrical sets
#===============================================================================

# point on trail rib
# point1=[rib_trail_up[-1][0]+airfoilDIST,rib_trail_up[-1][1]-airfoilTH/fairfoilTH], point2=[rib_trail_low[0][0]+airfoilDIST,rib_trail_low[0][1]+airfoilTH/fairfoilTH]
#xPoTE = ((rib_trail_up[-1][0]+airfoilDIST)+(rib_trail_low[0][0]+airfoilDIST))/2+0.001
#yPoTE = ((rib_trail_up[-1][1]-airfoilTH/fairfoilTH)+(rib_trail_low[0][1]+airfoilTH/fairfoilTH))/2


faces_setrib=WingModel.parts['Wing'].faces.findAt( 
    ((points_nose[1][0]-0.001, points_nose[1][1]-0.001, width) ,), 
    ((points_trail1_up[0][0]-0.001, points_trail1_up[0][1], width) ,),
    ((points_trail2_up[0][0]-0.001, points_trail2_up[0][1], width) ,),
    ((points_trail2_low[0][0]+0.001, points_trail2_low[0][1]+0.001, width) ,), 
    ((points_trail1_up[0][0]-0.001, points_trail1_up[0][1], width) ,),
    ((points_trail3_up[0][0]-0.001, points_trail3_up[0][1], width) ,),
    ((points_trail1_low[0][0]+0.001, points_trail1_low[0][1]+0.001, width) ,),
    ((points_trail3_low[0][0]+0.001, points_trail3_low[0][1]+0.001, width) ,), 
    ((points_trail1_low[0][0]+0.001, points_trail1_low[0][1]+0.001, width) ,),)
#    ((points_box_up[0][0]-0.001, points_box_up[0][1], width) ,),     #sim1,sim3                                                                     
#    ((points_box_up[0][0]-0.001, points_box_up[0][1]-0.001, width) ,),  #sim2
#    ((points_box_low[0][0]+0.001, points_box_low[0][1]+0.001, width) ,),)
                                                                                          
                  
                                                                                                  

for i in xrange(2,anzRibs+1): #anzRibs+1
    faces_setrib=faces_setrib+WingModel.parts['Wing'].faces.findAt(                   ( ( points_nose[1][0]-0.001, points_nose[1][1]-0.001,         width-(i-1)*anzDivision*width/Anz_Unterteilung) , ),
                                                                                                  ( ( points_trail1_up[0][0]-0.001, points_trail1_up[0][1], width-(i-1)*anzDivision*width/Anz_Unterteilung         ) , ),
                                                                                                  ( ( points_trail2_up[0][0]-0.001, points_trail2_up[0][1],       width-(i-1)*anzDivision*width/Anz_Unterteilung        ) , ),
                                                                                                  ( ( points_trail2_low[0][0]+0.001, points_trail2_low[0][1]+0.001,  width-(i-1)*anzDivision*width/Anz_Unterteilung             ) , ),
                                                                                                  ( (points_trail3_up[0][0]-0.001, points_trail3_up[0][1],     width-(i-1)*anzDivision*width/Anz_Unterteilung         ) , ),
                                                                                                  ( (points_trail1_low[0][0]+0.001, points_trail1_low[0][1]+0.001,                 width-(i-1)*anzDivision*width/Anz_Unterteilung    ) , ),
                                                                                                  ( (points_trail3_low[0][0]+0.001, points_trail3_low[0][1]+0.001, width-(i-1)*anzDivision*width/Anz_Unterteilung      ),  ),)
#                                                                                                  ( (points_box_up[0][0]-0.001, points_box_up[0][1],   width-(i-1)*anzDivision*width/Anz_Unterteilung      ),  ),  #sim1,sim3
#                                                                                                  ( (points_box_up[0][0]-0.001, points_box_up[0][1]-0.001,   width-(i-1)*anzDivision*width/Anz_Unterteilung      ),  ),  #sim2
#                                                                                                  ( (points_box_low[0][0]+0.001, points_box_low[0][1]+0.001, width-(i-1)*anzDivision*width/Anz_Unterteilung  ),  ),)
#Manual addition - Alejandro Valverde
faces_setrib=faces_setrib+WingModel.parts['Wing'].faces.findAt(                   ( ( points_nose[1][0]-0.001, points_nose[1][1]-0.001,         0.0) , ),
                                                                                              ( ( points_trail1_up[0][0]-0.001, points_trail1_up[0][1], 0.0         ) , ),
                                                                                              ( ( points_trail2_up[0][0]-0.001, points_trail2_up[0][1],       0.0        ) , ),
                                                                                              ( ( points_trail2_low[0][0]+0.001, points_trail2_low[0][1]+0.001,  0.0             ) , ),
                                                                                              ( (points_trail3_up[0][0]-0.001, points_trail3_up[0][1],     0.0         ) , ),
                                                                                              ( (points_trail1_low[0][0]+0.001, points_trail1_low[0][1]+0.001,                 0.0    ) , ),
                                                                                              ( (points_trail3_low[0][0]+0.001, points_trail3_low[0][1]+0.001, 0.0      ),  ),)


setrib=WingModel.parts['Wing'].Set(faces=faces_setrib, name='Set-Rib')

## cut wing (from 07_Assembly...)
#WingModel.parts['Wing'].PartitionFaceByDatumPlane(datumPlane=
#    WingModel.parts['Wing'].datums[front_nose_plane.id], faces=WingModel.parts['Wing'].faces.getByBoundingBox(
#    -0.01, -c, -0.01, x_spar1+0.01, c, 1.01*width) + WingModel.parts['Wing'].faces.getByBoundingBox(
#    x_spar2+0.01, -c, -0.01, c+0.01, c, 1.01*width))

WingModel.rootAssembly.regenerate()

WingModel.parts['Wing'].Set( name='wingtip' , edges=WingModel.parts['Wing'].edges.getByBoundingBox(-0.01*c, -c, width-0.01*c, 1.01*c, c, width+0.01*c))
WingModel.rootAssembly.regenerate()


    
anzRibs_rest = float(anzRibs-3)
anzDiv_sec3 = round(anzRibs_rest/2)
anzDiv_sec2 = anzRibs_rest-anzDiv_sec3
elements_rest = N-anzRibs
elements_per_division = (elements_rest-(elements_rest % anzRibs))/anzRibs
elements_at_root = elements_per_division+elements_rest % anzRibs


#Point on front spar
xPOFS = (mod7_UP[0]+mod7_DN[0])/2
yPOFS = (mod7_UP[1]+mod7_DN[1])/2
xPOFS1 = (mod7_UP[0]+3*mod7_DN[0])/4
yPOFS1 = (mod7_UP[1]+3*mod7_DN[1])/4
xPOFS3 = (3*mod7_UP[0]+mod7_DN[0])/4
yPOFS3 = (3*mod7_UP[1]+mod7_DN[1])/4

# Commented out, Alejandro Valverde
# for i in range(anzRibs-1):
#     faces_setspar1_ribs = WingModel.parts['Wing'].faces.findAt( ( ( xPOFS1, yPOFS1, width-(i+1)*element_width*anzDivision+0.001 ) ,), (( xPOFS1, yPOFS1, width-(i+1)*element_width*anzDivision-0.001 ), ),
# 		( ( xPOFS3, yPOFS3, width-(i+1)*element_width*anzDivision+0.001 ) ,), (( xPOFS3, yPOFS3, width-(i+1)*element_width*anzDivision-0.001 ), )) 
#     WingModel.parts['Wing'].Set(faces=faces_setspar1_ribs , name='Set_Spar1_Ribs'+str(anzRibs-1-i))	

# faces_setspar1_root = WingModel.parts['Wing'].faces.findAt( ( ( xPOFS1, yPOFS1, y_C[N]-1E-4 ) ,) , (( xPOFS3, yPOFS3, y_C[N]-1E-4 ) ,),)
# WingModel.parts['Wing'].Set(faces=faces_setspar1_root , name='Set_Spar1_root')
# faces_setspar1_tip = WingModel.parts['Wing'].faces.findAt( ( ( xPOFS1, yPOFS1, y_C[2*N-1]+1E-4 ) ,) , ( ( xPOFS3, yPOFS3, y_C[2*N-1]+1E-4 ) ,),)
# WingModel.parts['Wing'].Set(faces=faces_setspar1_tip , name='Set_Spar1_tip')



# faces_setspar1_sec4 = faces_setspar1_tip+WingModel.parts['Wing'].sets['Set_Spar1_Ribs'+str(anzRibs-1)].faces

# faces_setspar1_sec3 = WingModel.parts['Wing'].sets['Set_Spar1_Ribs'+str(anzRibs-2)].faces
# for i in range(int(anzDiv_sec3)-1):
#     faces_setspar1_sec3 = faces_setspar1_sec3+WingModel.parts['Wing'].sets['Set_Spar1_Ribs'+str(anzRibs-(2+i+1))].faces

# faces_setspar1_sec2 = WingModel.parts['Wing'].sets['Set_Spar1_Ribs'+str(2)].faces
# for i in range(int(anzDiv_sec2)-1):
#     faces_setspar1_sec2 = faces_setspar1_sec2+WingModel.parts['Wing'].sets['Set_Spar1_Ribs'+str(2+i+1)].faces

# faces_setspar1_sec1 = faces_setspar1_root+WingModel.parts['Wing'].sets['Set_Spar1_Ribs'+str(1)].faces



#set_Spar1_sec1=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec1 , name='Set_sec1')
#set_Spar1_sec2=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec2 , name='Set_sec2')
#set_Spar1_sec3=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec3 , name='Set_sec3')
#set_Spar1_sec4=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec4 , name='Set_sec4')
# faces_setspar1_sec1 = WingModel.parts['Wing'].faces.getByBoundingBox(x_spar1-0.001, -c, -0.01, x_spar1m+0.001, c, 1.01*width) #Commented out, Alejandro Valverde
#faces_setspar1_sec1 = WingModel.parts['Wing'].faces.getByBoundingBox(x_spar1-0.001, -c, -0.01, x_spar1+0.001, c, 1.01*width)
# set_Spar1_sec1=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec1 , name='Set_sec1') #Commented out, Alejandro Valverde
# set_Spar1_sec2=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec1 , name='Set_sec2') #Commented out, Alejandro Valverde
# set_Spar1_sec3=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec1 , name='Set_sec3') #Commented out, Alejandro Valverde
# set_Spar1_sec4=WingModel.parts['Wing'].Set(faces=faces_setspar1_sec1 , name='Set_sec4') #Commented out, Alejandro Valverde

faces_nodalpoints = WingModel.parts['Wing'].faces.getByBoundingBox(
    -0.01*c, -c, y_C[N]-1E-4, points_nose[0][0]+1E-04, c, 1.01*width)
    
#WingModel.parts['Wing'].PartitionFaceByDatumPlane(datumPlane=
#    WingModel.parts['Wing'].datums[front_nose_plane.id], faces=faces_nodalpoints)


#Point on rear spar
#xPORS = (mod8_UP[0]+mod8_DN[0])/2
#yPORS = (mod8_UP[1]+mod8_DN[1])/2
yPORS = yPOFS
xPORS = numpy.interp(yPOFS, [mod8_DN[1],mod8_UP[1]], [mod8_DN[0],mod8_UP[0]])


# Sets for cross sections at wing tip for shear center determination
#vertex1_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_nose[0][0], 0, width), 1E-05) # 0
#vertex2_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_box_up[0][0], 0, width), 1E-05) # points_trail3_up[0][0]
#vertex1_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_nose[0][0], (mod7_UP[1]+mod7_DN[1])/2, width), 1E-05) # 0
#vertex2_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_box_up[0][0], (mod7_UP[1]+mod7_DN[1])/2, width), 1E-05) # points_trail3_up[0][0]
vertex1_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((xPOFS, yPOFS, width), 1E-05) # 0
vertex2_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((xPORS, yPORS, width), 1E-05) # points_trail3_up[0][0]
vertex3_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_nose[0][0], points_nose[0][1], width), 1E-05)
vertex4_section_SC = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_box_up[0][0], points_box_up[0][1], width), 1E-05)
WingModel.parts['Wing'].Set(name='MeasuringPoint'+str(anzRibs)+'_Front' , vertices=vertex1_section_SC)
WingModel.parts['Wing'].Set(name='MeasuringPoint'+str(anzRibs)+'_Rear' , vertices=vertex2_section_SC)
WingModel.parts['Wing'].Set(name='Rib'+str(anzRibs)+'_SC_vertices' , vertices=vertex1_section_SC+vertex2_section_SC)
WingModel.parts['Wing'].Set(name='LoadingPoint'+str(anzRibs)+'_Front' , vertices=vertex3_section_SC)
WingModel.parts['Wing'].Set(name='LoadingPoint'+str(anzRibs)+'_Rear' , vertices=vertex4_section_SC)

# Sets for cross sections at control points for reading out angles
for i in range(N):
    vertex1_section = WingModel.parts['Wing'].vertices.getByBoundingSphere((0, 0, y_C[N+i]), 1E-04)     #((rib_box_up[-1][0], 0, y_C[N+i]), 1E-04)
    vertex2_section = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_trail4_up[0][0], 0, y_C[N+i]), 1E-04)    #((rib_box_up[0][0], 0, y_C[N+i]), 1E-04)
    WingModel.parts['Wing'].Set(name='MeasuringPointSection'+str(i+1)+'_Front' , vertices=vertex1_section)
    WingModel.parts['Wing'].Set(name='MeasuringPointSection'+str(i+1)+'_Rear' , vertices=vertex2_section)
    
# Sets for cross sections at ribs for reading out angles
for i in range(anzRibs):
    vertex1_rib = WingModel.parts['Wing'].vertices.getByBoundingSphere((0, 0, width-i*element_width*anzDivision), 1E-04)    # ((rib_box_up[-1][0], 0, width-i*element_width*anzDivision), 1E-04)
    vertex2_rib = WingModel.parts['Wing'].vertices.getByBoundingSphere((points_trail4_up[0][0], 0, width-i*element_width*anzDivision), 1E-04)   # ((rib_box_up[0][0], 0, width-i*element_width*anzDivision), 1E-04)
    WingModel.parts['Wing'].Set(name='MeasuringPointRib'+str(anzRibs-i)+'_Front' , vertices=vertex1_rib)
    WingModel.parts['Wing'].Set(name='MeasuringPointRib'+str(anzRibs-i)+'_Rear' , vertices=vertex2_rib)


##############################################################
#===============================================================================
# INCLUSION OF WING-BOX
#===============================================================================

del mdb.models['Model-1'] #Delete old model

#Find faces of rear spar
# zRange = np.linspace(0.01*width, 0.99*width, 100)
# yRange = np.linspace(0.05*mod8_DN[1], 0.95*mod8_UP[1], 4)
# faces_list=[]
# for zLoop in zRange: #xrange is the same as range
#     for yLoop in yRange:
#         face_found=WingModel.parts['Wing'].faces.findAt(( ( rib_box_low[-1][0]   ,yLoop,zLoop ) , )  , )
        
#         if face_found not in faces_list: #if the face hasn't already been found
#             faces_list.append(face_found)
#         else:
#             print('Face found is being neglected, it was already found')

# faces_tuple = tuple(faces_list)

# WingModel.parts['Wing'].Set(faces=faces_tuple , name='setToDeleteWingBox')

# sets to remove
setsToRemove = ('fullpart_set_Spar1', 'fullpart_set_Spar1_Low', 'fullpart_set_Spar1_Up')
WingModel.parts['Wing'].RemoveFaces(deleteCells=False, 
    faceList=list(WingModel.parts['Wing'].sets['fullpart_set_Spar1_Low'].faces))
WingModel.parts['Wing'].RemoveFaces(deleteCells=False, 
    faceList=list(WingModel.parts['Wing'].sets['fullpart_set_Spar2'].faces))

for setToRemove in setsToRemove:
    del WingModel.parts['Wing'].sets[setToRemove]

if sys.version_info.major == 2:
    execfile('mainBuildAndExecuteWingBox.py') #Load parametric study values
elif sys.version_info.major == 3:
    exec(open("./"+'mainBuildAndExecuteWingBox.py').read())

#Remove upper and lower flanges of wing-box
# faces_list=[]
# faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( design.cutWingTip - (0.1*design.cutWingTip),yInSkinUP2(design, design.C3/2), design.C3/2 ) , )  , )[0]) #[0] access the Face object
# faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( design.cutWingTip - (0.1*design.cutWingTip),yInSkinDN2(design, design.C3/2), design.C3/2 ) , )  , )[0]) #[0] access the Face object
# faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( design.cutWingTip - (0.1*design.cutWingTip),yInSkinUP2(design, design.B/2), design.B/2 ) , )  , )[0]) #[0] access the Face object
# faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( design.cutWingTip - (0.1*design.cutWingTip),yInSkinDN2(design, design.B/2), design.B/2 ) , )  , )[0]) #[0] access the Face object

# if (design.innerRibs_gap + design.B) == mesh.d:
#     for xPosInner in design.innerRibsXpos:
#         if xPosInner == design.innerRibsXpos[0]: #For the rib close to the root

#             positionInX = xPosInner - ((xPosInner - design.cutWingRoot)/2)

#         else:

#             positionInX = xPosInner - ((xPosInner - design.innerRibsXpos[index-1])/2)
#         faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( xPosInner - (0.1*design.cutWingTip),yInSkinUP2(design, design.C3/2), design.C3/2 ) , )  , )[0]) #[0] access the Face object
#         faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( xPosInner - (0.1*design.cutWingTip),yInSkinDN2(design, design.C3/2), design.C3/2 ) , )  , )[0]) #[0] access the Face object
#         # faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( xPosInner - (0.1*design.cutWingTip),yInSkinUP2(design, design.B/2), design.B/2 ) , )  , )[0]) #[0] access the Face object
#         # faces_list.append(WingModel.parts['RibBox'].faces.findAt(( ( xPosInner - (0.1*design.cutWingTip),yInSkinDN2(design, design.B/2), design.B/2 ) , )  , )[0])

#         index += 1
# WingModel.parts['RibBox'].RemoveFaces(deleteCells=False, faceList=faces_list)

model.rootAssembly.rotate(angle=-90.0, axisDirection=
    (0.0, 1.0, 0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Lattice-1', ))
model.rootAssembly.translate(instanceList=('Lattice-1', ), vector=(0.0, 0.0, -design.cutWingRoot))
model.rootAssembly.translate(instanceList=('Lattice-1', ), vector=(c_0*PositionRearSpar, 0.0, 0.0))
model.rootAssembly.translate(instanceList=('Lattice-1', ), vector=(0.0, mod8_DN[1]-design.cutDown, 0.0))

#Partition in wing skin for creating nodes at the desired position
datumPlane_partitionSkin = WingModel.parts['Wing'].DatumPlaneByPrincipalPlane(offset=c_0*PositionRearSpar - (design.B/2), principalPlane=YZPLANE)
WingModel.parts['Wing'].PartitionFaceByDatumPlane(datumPlane=WingModel.parts['Wing'].datums[datumPlane_partitionSkin.id], faces=WingModel.parts['Wing'].sets['box_set'].faces)

WingModel.rootAssembly.InstanceFromBooleanMerge(domain=
    GEOMETRY, instances=(WingModel.rootAssembly.instances['Wing-1'], 
    # WingModel.rootAssembly.instances['Lattice-1'],
    WingModel.rootAssembly.instances['RibBox-1'])
    , name='Wing-noLattice', originalInstances=SUPPRESS)

WingModel.rootAssembly.InstanceFromBooleanMerge(domain=
    GEOMETRY, instances=(WingModel.rootAssembly.instances['Wing-noLattice-1'], 
    WingModel.rootAssembly.instances['Lattice-1'])
    , name='Wing-final', originalInstances=SUPPRESS)

finalPart = WingModel.parts['Wing-final']
#################################################################

#===============================================================================
# 9-Section definition - A lot of commented code hasn't been included
#===============================================================================

WingModel.HomogeneousShellSection(idealization=NO_IDEALIZATION,     # Alu
    integrationRule=SIMPSON, material='Aluminium_wingSkin', name='stringer_section', 
    numIntPts=5, poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=
    GRADIENT, thickness=t_stringer, thicknessField='', thicknessModulus=None, 
    thicknessType=UNIFORM, useDensity=OFF)

WingModel.CompositeShellSection(idealization=NO_IDEALIZATION,         # UD
    integrationRule=SIMPSON, layup=(SectionLayer(thickness=0.001, orientAngle=0,
    material='CFRP-UD-Box'),), name='stringer_section_UD', poissonDefinition=DEFAULT, 
    preIntegrate=OFF, symmetric=False, temperature=GRADIENT, thicknessModulus=
    None, thicknessType=UNIFORM, useDensity=OFF)

## Spar-Rear-Layup thicker    
WingModel.CompositeShellSection(idealization=NO_IDEALIZATION,    # is updated before job submission
    integrationRule=SIMPSON, layup=(SectionLayer(thickness=1, orientAngle=0, material='Glass-Fabric', plyName='ply1111'),), 
    name='Spar-Rear_GLASS_ONLY', poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=
    GRADIENT, thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

## Spar-Rear-Layup thicker    
WingModel.CompositeShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, layup=(SectionLayer(thickness=0.00002, orientAngle=FaserWinkel, material='Glass-Fabric', plyName='ply1111'), 
    SectionLayer(thickness=thickness_bucklingSpar, orientAngle=FaserWinkel, material='CFRP-UD-Box', plyName='ply2222'),
    SectionLayer(thickness=0.00002, orientAngle=FaserWinkel, material='Glass-Fabric', plyName='ply3333'),), 
    name='Spar-Rear-Layup_thick', poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=
    GRADIENT, thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

fakt=1.0                                                                                                     # old layup for box
WingModel.CompositeShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, layup=(SectionLayer(thickness=fakt*0.000125, 
    material='CFRP-UD-Box'), SectionLayer(thickness=fakt*0.000125, 
    material='CFRP-UD-Box'), SectionLayer(thickness=fakt*0.000125, orientAngle=90.0, material='CFRP-UD-Box'), 
    SectionLayer(thickness=fakt*0.000125, orientAngle=90.0, material='CFRP-UD-Box'), 
    SectionLayer(thickness=fakt*0.000125, orientAngle=90.0, 
    material='CFRP-UD-Box'), SectionLayer(thickness=fakt*0.000125, 
    orientAngle=90.0, material='CFRP-UD-Box'), SectionLayer(
    thickness=fakt*0.000125, material='CFRP-UD-Box'), SectionLayer(
    thickness=fakt*0.000125, material='CFRP-UD-Box'),), name='Box-LayUp', 
    poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=
    GRADIENT, thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

# cfk box
# finalPart.SectionAssignment(offset=0.0, 
#     offsetField='', offsetType=MIDDLE_SURFACE, region=
#     finalPart.sets['fullpart_set_Spar2'], sectionName=
#     'Spar-Rear-Layup_thick', thicknessAssignment=FROM_SECTION)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['Set-Rib'], sectionName=
    'stringer_section', thicknessAssignment=FROM_SECTION)

#Buckling Spar
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['box_set_oben'], sectionName=
    'Box-LayUp', thicknessAssignment=FROM_SECTION)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['box_set_unten'], sectionName=
    'Box-LayUp', thicknessAssignment=FROM_SECTION)
# finalPart.SectionAssignment(offset=0.0, 
#     offsetField='', offsetType=MIDDLE_SURFACE, region=
#     finalPart.sets['Set_sec1'], sectionName=
#     'Box-LayUp', thicknessAssignment=FROM_SECTION)

finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['front_set'], sectionName=
    'Box-LayUp', thicknessAssignment=FROM_SECTION)


finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['rear_set_unten1'], sectionName=
    'Box-LayUp', thicknessAssignment=FROM_SECTION)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['rear_set_unten2'], sectionName=
    'Box-LayUp', thicknessAssignment=FROM_SECTION)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['rear_set_oben'], sectionName=
    'Box-LayUp', thicknessAssignment=FROM_SECTION)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    finalPart.sets['rear_end'], sectionName=
    'Box-LayUp', thicknessAssignment=FROM_SECTION)

#===============================================================================
# 10-Meshing
#===============================================================================

mesh_wing = structtype()
mesh_wing.courseSize=float(paraRead_wing.sizeCourse)
mesh_wing.fineSize=float(paraRead_wing.sizeFine)
mesh_wing.ElemType = mesh.ElemType
mesh_wing.d = mesh.d
#sizeInner=0.005#0.01#0.02
#sizeAll=0.01
# finalPart.seedPart(deviationFactor=0.1,minSizeFactor=0.1, size=fineSize)

finalPart.seedPart(deviationFactor=0.1,minSizeFactor=0.1, size=mesh_wing.courseSize)

# partToMesh.seedEdgeBySize(constraint = FINER, deviationFactor=0.1, edges=setOfEdges, size=mesh_wing.courseSize)

# setMeshControls                jau dat klappt

try:
    for i in range(Anz_Unterteilung):
      # WingModel.parts['Wing'].setMeshControls(regions=WingModel.parts['Wing'].faces.findAt(   ( ( points_box_up[0][0] ,0.0, y_C[N+i] ) , ) , )  ,technique=STRUCTURED) # spar
      #  WingModel.parts['Wing'].setMeshControls(regions=WingModel.parts['Wing'].faces.findAt(   ( ( points_nose[0][0] ,0.0, y_C[N+i] ) , ) , )  ,technique=STRUCTURED) # spar
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_box_up[0][0] ,points_box_up[0][1], y_C[N+i] ) , ) , )  ,technique=STRUCTURED) # spar
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_nose[0][0] ,points_nose[0][1], y_C[N+i] ) , ) , )  ,technique=STRUCTURED) # spar
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_nose[10][0]   ,points_nose[10][1], y_C[N+i] ) , ) , )  ,technique=STRUCTURED) # nose
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_nose[-10][0]   ,points_nose[-10][1], y_C[N+i] ) , ) , )  ,technique=STRUCTURED) # nose
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_trail1_up[1][0]   ,points_trail1_up[1][1], y_C[N+i]  ) , ) , )  ,technique=STRUCTURED)
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_trail1_low[1][0]   ,points_trail1_low[1][1], y_C[N+i] ) , ) , )  ,technique=STRUCTURED)
    ########    finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( c  ,0.0001,y*width/Anz_Unterteilung+0.001) , ) , )  ,technique=STRUCTURED)
    ########    finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( c  ,-0.0001,y*width/Anz_Unterteilung+0.001) , ) , )  ,technique=STRUCTURED)
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_box_up[1][0]   ,points_box_up[1][1], y_C[N+i] ) , ) , )  ,technique=STRUCTURED) # box
        finalPart.setMeshControls(regions=finalPart.faces.findAt(   ( ( points_box_low[-2][0]   ,points_box_low[-2][1], y_C[N+i] ) , ) , )  ,technique=STRUCTURED)

##finalPart.setMeshControls(regions=faces_setrib,technique=STRUCTURED) # rippe
except Exception as e:
    print('Error: Error found when defining the mesh controls')

# setElementType             jau dat klappt

# finalPart.setElementType(elemTypes=(ElemType(elemCode=S4R, elemLibrary=STANDARD), ElemType(elemCode=S3, elemLibrary=STANDARD),), regions=finalPart.sets['fullpart_set_Spar1'])
# finalPart.setElementType(elemTypes=(ElemType(elemCode=S4R, elemLibrary=STANDARD), ElemType(elemCode=S3, elemLibrary=STANDARD),), regions=finalPart.sets['fullpart_set_Spar2'])
# finalPart.setElementType(elemTypes=(ElemType(elemCode=S4R, elemLibrary=STANDARD), ElemType(elemCode=S3, elemLibrary=STANDARD),), regions=finalPart.sets['C-box'])
# finalPart.setElementType(elemTypes=(ElemType(elemCode=S4R, elemLibrary=STANDARD), ElemType(elemCode=S3, elemLibrary=STANDARD),), regions=setrib)    # rippe

## Mesh lattice
edges_lattice = finalPart.edges.getByBoundingBox(c_0*PositionRearSpar - 1.05*design.B,
                                                mod8_DN[1] + (0.75*design.cutGap_y),
                                                0.0,
                                                c_0*PositionRearSpar + 0.05*design.B,
                                                mod8_UP[1] - (0.75*design.cutGap_y),
                                                width)

finalPart.Set(name='edges_lattice_fine_mesh' , edges=edges_lattice)

# for edge in edges_lattice:
finalPart.seedEdgeBySize(constraint = FINER, deviationFactor=0.1, edges=edges_lattice, size=mesh.fineSize)

# #Partition on Skin
# datumPlane_partitionSkin = finalPart.DatumPlaneByPrincipalPlane(offset=c_0*PositionRearSpar - 1.05*design.B, principalPlane=YZPLANE)
# finalPart.PartitionFaceByDatumPlane(datumPlane=finalPart.datums[datumPlane_partitionSkin.id], faces=finalPart.sets['box_set'].faces)

if False: #Additional option to fine mesh area on the wing skin on the top of the lattice
    try:
        zRange = np.linspace(0.01*width, 0.95*width, 20)
        edges_tuple_fineMesh_Skin=()
        for zLoop in zRange: #xrange is the same as range
            # edge_found1=finalPart.edges.findAt(( ( c_0*PositionRearSpar - 1.05*design.B,ySkinForXUP(pointsSkin, c_0*PositionRearSpar - 1.05*design.B),zLoop ) , )  , )
            # edge_found2=finalPart.edges.findAt(( ( c_0*PositionRearSpar,ySkinForXUP(pointsSkin, c_0*PositionRearSpar),zLoop ) , )  , )
            # edge_found3=finalPart.edges.findAt(( ( c_0*PositionRearSpar - 1.05*design.B,ySkinForXDN(pointsSkin, c_0*PositionRearSpar - 1.05*design.B),zLoop ) , )  , )
            # edge_found4=finalPart.edges.findAt(( ( c_0*PositionRearSpar,ySkinForXDN(pointsSkin, c_0*PositionRearSpar),zLoop ) , )  , )
            edge_found1=finalPart.edges.getByBoundingSphere( (c_0*PositionRearSpar - 1.05*design.B,ySkinForXUP(pointsSkin, c_0*PositionRearSpar - 1.05*design.B),zLoop), 0.001 )
            edge_found2=finalPart.edges.getByBoundingSphere( (c_0*PositionRearSpar,ySkinForXUP(pointsSkin, c_0*PositionRearSpar),zLoop), 0.001 )
            edge_found3=finalPart.edges.getByBoundingSphere( (c_0*PositionRearSpar - 1.05*design.B,ySkinForXDN(pointsSkin, c_0*PositionRearSpar - 1.05*design.B),zLoop), 0.001 )
            edge_found4=finalPart.edges.getByBoundingSphere( (c_0*PositionRearSpar,ySkinForXDN(pointsSkin, c_0*PositionRearSpar),zLoop), 0.001 )
            
            if edge_found1 and edge_found1 not in edges_tuple_fineMesh_Skin: #if the edge hasn't already been found
                edges_tuple_fineMesh_Skin += (edge_found1, )
            else:
                print('Edge found is being neglected, it was already found')
            if edge_found2 and edge_found2 not in edges_tuple_fineMesh_Skin: #if the edge hasn't already been found
                edges_tuple_fineMesh_Skin += (edge_found2, )
            else:
                print('Edge found is being neglected, it was already found')
            if edge_found3 and edge_found3 not in edges_tuple_fineMesh_Skin: #if the edge hasn't already been found
                edges_tuple_fineMesh_Skin += (edge_found3, )
            else:
                print('Edge found is being neglected, it was already found')
            if edge_found4 and edge_found4 not in edges_tuple_fineMesh_Skin: #if the edge hasn't already been found
                edges_tuple_fineMesh_Skin += (edge_found4, )
            else:
                print('Edge found is being neglected, it was already found')

        finalPart.Set(name='edges_skin_fine_mesh' , edges=edges_tuple_fineMesh_Skin)
        for j in range(len(edges_tuple_fineMesh_Skin)):
            finalPart.seedEdgeBySize(constraint = FINER, deviationFactor=0.1, edges=edges_tuple_fineMesh_Skin[j], size=mesh.fineSize)
    except Exception as e:
        print('Warning: Set for finer mesh on skin could not be created')

## generate mesh

finalPart.generateMesh()
WingModel.rootAssembly.regenerate()

#===============================================================================
# 11-Material orientations - Not included because it is only necessary for laminated with different lay-ups
#===============================================================================

#===============================================================================
# 12-Sets for output
#===============================================================================

for i in range(Anz_Unterteilung):   # REIHENFOLGE : OBEN HINTEN ; OBEN VORNE ; UNTEN VORNE ; UNTEN HINTEN
    WingModel.rootAssembly.Set(name='Point-OL-'+str(i),  vertices=WingModel.rootAssembly.instances['Wing-1'].vertices.getByBoundingSphere((points_nose[0][0], points_nose[0][1], y_C[N+i]), 1E-03))
    WingModel.rootAssembly.Set(name='Point-OR-'+str(i),  vertices=WingModel.rootAssembly.instances['Wing-1'].vertices.getByBoundingSphere((points_box_up[0][0], points_box_up[0][1], y_C[N+i]), 1E-03))
    WingModel.rootAssembly.Set(name='Point-UR-'+str(i),  vertices=WingModel.rootAssembly.instances['Wing-1'].vertices.getByBoundingSphere((points_box_low[-1][0], points_box_low[-1][1], y_C[N+i]), 1E-03))
    WingModel.rootAssembly.Set(name='Point-UL-'+str(i),  vertices=WingModel.rootAssembly.instances['Wing-1'].vertices.getByBoundingSphere((points_nose[-1][0], points_nose[-1][1], y_C[N+i]), 1E-03))
##
#    node_P_OR = WingModel.rootAssembly.sets['Point-OR-'+str(i)].nodes
#    WingModel.rootAssembly.Set(name='NODES-P_OR-'+str(i), nodes=node_P_OR)
#    node_P_UR = WingModel.rootAssembly.sets['Point-UR-'+str(i)].nodes
#    WingModel.rootAssembly.Set(name='NODES-P_UR-'+str(i), nodes=node_P_UR)
#    node_P_OL = WingModel.rootAssembly.sets['Point-OL-'+str(i)].nodes
#    WingModel.rootAssembly.Set(name='NODES-P_OL-'+str(i), nodes=node_P_OL)
#    node_P_UL = WingModel.rootAssembly.sets['Point-UL-'+str(i)].nodes
#    WingModel.rootAssembly.Set(name='NODES-P_UL-'+str(i), nodes=node_P_UL)



for y in xrange(1,anzRibs+1):   
    WingModel.rootAssembly.Set(name='LEAD-'+str(int(anzRibs)-y),  vertices=WingModel.rootAssembly.instances['Wing-1'].vertices.getByBoundingSphere((0, 0, width-(y-1)*anzDivision*width/Anz_Unterteilung), 1E-03))     
    WingModel.rootAssembly.Set(name='REAR-'+str(int(anzRibs)-y),  vertices=WingModel.rootAssembly.instances['Wing-1'].vertices.getByBoundingSphere((points_trail4_up[0][0], 0, width-(y-1)*anzDivision*width/Anz_Unterteilung), 1E-03))
#    node_LEAD = WingModel.rootAssembly.sets['LEAD-'+str(int(anzRibs)-y)].nodes
#    WingModel.rootAssembly.Set(name='NODES-LEAD-'+str(int(anzRibs)-y), nodes=node_LEAD)
#    node_REAR = WingModel.rootAssembly.sets['REAR-'+str(int(anzRibs)-y)].nodes
#    WingModel.rootAssembly.Set(name='NODES-REAR-'+str(int(anzRibs)-y), nodes=node_REAR)
    
## Point on which the testload acts; check with sections defined in 07_Assembly_Dom
#WingModel.rootAssembly.Set(name='point_SC',  vertices=WingModel.rootAssembly.instances['Wing-1'].vertices.getByBoundingBox(x_SC-1E-4, 0.2*0.12*c, width-1E-4, x_SC+1E-4, 0.12*c, width+1E-4))

WingModel.rootAssembly.regenerate()