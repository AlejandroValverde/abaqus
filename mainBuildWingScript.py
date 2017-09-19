from scipy.optimize import fsolve
import math
import pdb #pdb.set_trace()
import os
import numpy as np
import time

class structtype():
	
	def __init__(self):

		pass

	def loadParameters(self, fileName):

		attributes = ()

		file = open(fileName, 'r')

		lines = file.readlines()

		for i in range(int((len(lines))/2)):

			nameParater = lines[(i*2)]
			valueParater = lines[(2*i)+1]

			valueParater = valueParater.replace('\n','')
			nameParater = nameParater.replace('\n','')

			valueParater = valueParater.replace('\r\n','')
			nameParater = nameParater.replace('\r\n','')

			setattr(self, nameParater, valueParater)

			attributes += (nameParater, )

		file.close()

		self.attrs = attributes

	def createInputMatlab(self, fileName, design_wing, design):

		# imp
		# design_wing.C2r, design_wing.s, design.M, design.N, design.cutGap_y

		file = open(fileName, 'w')

		for atr in self.attrs:

			file.write(atr + ';')

		file.write('C2r' + ';')
		file.write('s' + ';')
		file.write('L0' + ';')
		file.write('r0' + ';')

		file.write('\n')
		
		for atr in self.attrs:

			file.write(getattr(self, atr) + ';')

		file.write(str(design_wing.C2r) + ';')
		file.write(str(design_wing.s) + ';')
		file.write(str(design.L0) + ';')
		file.write(str(design.r0) + ';')

		file.close()

def parametersFromWing(design_wing):

	span = design_wing.span
	c_0 = design_wing.c_0
	N = design_wing.N
	PositionFrontSpar = design_wing.PositionFrontSpar
	PositionRearSpar = design_wing.PositionRearSpar

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

	Anz_Unterteilung = N                                                            # Check for consisteny with N in parameters
	c = c_0                                                                         # chord length [m], Check for consisteny with c_0 in parameters
	##rectangular wing, the chord lengths are the same along the spannwise direction
	width = span/2                                                                     # width of the wing (half of the spanwidth) [m], check for consisteny with b in parameters
	x_spar1 = PositionFrontSpar*c     #0.25*c     #0.16*c     # 0.21*c                                   # position of spar on chord
	x_spar2 = PositionRearSpar*c#(PositionFrontSpar+PositionRearSpar+AngleRearSpar)*c     #0.35*c     #0.26*c     # 0.31*c                          # position of spar on chord
	x_spar1m = PositionFrontSpar*c     #(PositionFrontSpar+AngleFrontSpar)*c                               # position of spar on chord
	x_spar2m = PositionRearSpar*c#(PositionFrontSpar+PositionRearSpar)*c                            # position of spar on chord

	anzRibs = 6
	t_stringer = 0.001
	element_width = width/N
	x_str1 = x_spar2+ (c-x_spar2)/4	                # position of spar on chord
	x_str2 = x_spar2+ (c-x_spar2)/2           # position of spar on chord1
	x_str3 = c-0.02*c #x_spar2+ 3*(c-x_spar2)/4

	#===============================================================================
	# 2-PointsAndSurfaces
	#===============================================================================

	#Import airfoil points
	Cpxy = np.genfromtxt('NACA0012_199points.dat', delimiter=None, skip_header=1)

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

	skin_points_upX = [[0] for i in range(int(n_skinpts/2))]
	skin_points_upY = [[0] for i in range(int(n_skinpts/2))]
	for i in range(int(n_skinpts/2)):
		skin_points_upX[i] = x[i]*c_0
		skin_points_upY[i] = y[i]*c_0

	#skin_points_upX.sort(reverse=True)
	#skin_points_upY.sort(reverse=True)
	#skin_points_upX.sort()
	#skin_points_upY.sort()	
	skin_points_upX = np.flipud(skin_points_upX) 
	skin_points_upY = np.flipud(skin_points_upY) 
		
	skin_points_dnX = [[0] for i in range(int(n_skinpts/2))]
	skin_points_dnY = [[0] for i in range(int(n_skinpts/2))]
	for i in range(int(n_skinpts/2)):
		j = i+int(n_skinpts/2)
		skin_points_dnX[i] = x[j]*c_0
		skin_points_dnY[i] = y[j]*c_0

	##define nodes to find edges in 06_modelling...
	mod6_UP_r = [0,0];
	mod6_UP_m = [0,0];
	mod6_UP_f = [0,0];
	distSPARsectUP = x_spar2m - x_spar1
	distSPARsectDN = x_spar2 - x_spar1m
	for i in range(int(n_skinpts/2)):	
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
	for i in range(int(n_skinpts/2)+1,n_skinpts):	
		if x[i]*c_0 > x_spar1m + distSPARsectDN/2 and mod6_DN_m[0] == 0:
			mod6_DN_m[0] = x[i]*c_0;
			mod6_DN_m[1] = y[i]*c_0;


	mod7_UP = [x_spar1, np.interp(x_spar1, skin_points_upX, skin_points_upY)]
	mod7_DN = [x_spar1m, np.interp(x_spar1m, skin_points_dnX, skin_points_dnY)]

	mod8_UP = [x_spar2m, np.interp(x_spar2m, skin_points_upX, skin_points_upY)]
	mod8_DN = [x_spar2, np.interp(x_spar2, skin_points_dnX, skin_points_dnY)]

	# out

	design_wing.C2r = mod8_UP[1] - mod8_DN[1]
	design_wing.C2f = mod7_UP[1] - mod7_DN[1]
	design_wing.C2diff = (design_wing.C2f - design_wing.C2r) / 2
	design_wing.s = width
	design_wing.C3 = (PositionRearSpar*c_0) - (PositionFrontSpar*c_0)

	return design_wing

def readSearchLandrForWing(fileName):

	dataDict = {}

	file = open(fileName, 'r')

	lines = file.readlines()

	for i in range(int(len(lines)/2)):

		nameParater = lines[(i*2)]
		valueParater = lines[(2*i)+1]

		valueParater = valueParater.replace('\r\n','')
		nameParater = nameParater.replace('\r\n','')

		valueParater = valueParater.replace('\n','')
		nameParater = nameParater.replace('\n','')

		dataDict[nameParater] = float(valueParater)

	file.close()

	return dataDict

def f(x, *dimFix):

	L, r = x

	C2r, s, M, N, cutGap_y = dimFix

	eq1 = C2r - (2*(M-1)*math.sqrt( (3/4)* ( math.pow((2*math.pow(r,2)),2) + math.pow((2*math.pow(L,2)),2) ) )) - (2 * (math.pow(r,2) + cutGap_y) )

	eq2 = s - ((N - 1)*math.sqrt( math.pow((2*math.pow(r,2)),2) + math.pow((2*math.pow(L,2)),2)  ))

	return (eq1, eq2)

#####################
print('\n'+'### Build wing model')
# Obtain parameters

######## Chiral
paraRead = structtype()
paraRead_wing = structtype()
design = structtype()
design_wing = structtype()

#Read parameters
paraRead.loadParameters('inputAbaqus.txt')

# pdb.set_trace()
design.M = int(paraRead.M) #Number of unit cells in transversal direction
design.N = int(paraRead.N) #20 - Number of unit cells in spanwise direction
design.L0 = float(paraRead.L) #half length, initial values
design.r0 = float(paraRead.r)  #Node radius, initial values
design.cutGap_y = float(paraRead.cutGap_y) #Gap between the lattice and the skin in the y direction

######## Wing
paraRead_wing.loadParameters('inputAbaqus_wing.txt')

design_wing.span = float(paraRead_wing.span)
design_wing.c_0 = float(paraRead_wing.c_0)
design_wing.N = int(paraRead_wing.N)
design_wing.PositionFrontSpar = float(paraRead_wing.PositionFrontSpar)
design_wing.PositionRearSpar = float(paraRead_wing.PositionRearSpar)

design_wing = parametersFromWing(design_wing) #Obtain rest of parameters from wing

dimFix = (design_wing.C2r, design_wing.s, design.M, design.N, design.cutGap_y)

# r0 = 5
# L0 = 13

[sol,infodict,ier,msg] = fsolve(f, (math.sqrt(design.L0), math.sqrt(design.r0)), args=dimFix, full_output=1)

if ier == 1:

	print('-> Solution found')

	print('--- > L: '+str(math.pow(sol[0], 2)))
	print('--- > r: '+str(math.pow(sol[1], 2)))

	file = open('findrAndL_found.txt', 'w')

	values = (math.pow(sol[0], 2), math.pow(sol[1], 2)) 
	labels = ('L', 'r')

	for label, value in zip(labels, values):
		
		file.write(label + '\n')

		file.write(str(value) + '\n')

	file.close()

else:
	print('-> Not solution found for r and L, message: '+msg)
	print('-> Using Matlab now')

	paraRead.createInputMatlab('inputAbaqus_Matlab.txt', design_wing, design)

	os.system('matlab -nodesktop -nosplash -r "searchParameters_Matlab;quit;"')

	#Wait for matlab to finish its execution
	i = 1
	flagMatlabExecution = False
	while True:
		time.sleep(20)

		print('-> Waiting for Matlab to finish, '+str(i*20)+'seconds elapsed')
		i += 1
		
		if i > 7:
			raise ValueError('Error: Aborting Matlab execution, it is taking so much time')

		for f in os.listdir(os.getcwd()):
		    if f.startswith('Matlab_execution.txt'):
		        
		        file = open('Matlab_execution.txt', 'r')

		        lines = file.readlines()

		        for i in range(int((len(lines))/2)):

		        	execIndicator = lines[(2*i)+1]

		        	execIndicator = execIndicator.replace('\n','')
		        	execIndicator = execIndicator.replace('\r\n','')

		        	file.close()
		        	# os.remove(f)

		        	if execIndicator == 'OK':
		        		flagMatlabExecution = True

		        	elif execIndicator == 'FAIL':
		        		raise ValueError('Error: Matlab execution failed')

		        	else:
		        		raise ValueError('Not possible to read Matlab exe file')

		if flagMatlabExecution:
			print('-> Matlab successfully executed')
			break #Exit while loop

out = os.system('abaqus cae script=codeForWing_temp.py')