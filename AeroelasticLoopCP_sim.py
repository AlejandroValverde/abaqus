## Previous Alejandro Valverde
from moduleBuildWingBoxAndPostProc import loadParameters
# class structtype():
#     pass

# paraRead_wing = structtype()
# paraRead_wing = loadParameters(paraRead_wing, 'inputAbaqus_wing.txt')

t_step = 0.1 #float(paraRead_wing.t_step)

## AEROELASTIC LOOP

print('BEGINNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN')

D_lines = [line.rstrip('\n') for line in open(folder_path+'/Damping'+str(iSIM)+'_'+str(elements_inWeb)+'.txt')]
DD_lines=map(int, D_lines[0].split())
Damp_lines=DD_lines[0]
##os.remove(folder_path+'/Damping.txt')

t_tot = 1
# t_step = 0.2#1.0#1#0.15 #FibrOrientation # 0.25
t_step_0 = t_step #t_step defined in codeForWing_sim

afterBUCKLING = 0

velocity = V
Ma=velocity/340.29
aeroQ = 0.5*rho*(velocity**2) # dynamic pressure [Pa]
chord = c_0

alpha = alpha_init

displTE_old = 0.0
displTE = [0.0]*N
displLE = [0.0]*N
twist = [0.0]*N

# step parameter
dampCoeff_Step1 = Damp_lines
dampMagnit=Damp_lines
initialIncrement = 0.005
varsFieldOut = ('U', 'UR', 'RF', 'P','EE','S','E','DAMAGET', 'CFAILURE', 'COORD')
varsHistoryOut = ('ALLIE','ALLSD')

###
model = mdb.models['Model-SparAngle-'+str(int(jobNumber))]
assembly = mdb.models['Model-SparAngle-'+str(int(jobNumber))].rootAssembly
instance = mdb.models['Model-SparAngle-'+str(int(jobNumber))].rootAssembly.instances['Wing-final-1']

# fail stress CFRP
##model.materials['CFRP-UD-Box'].elastic.FailStress(table=((
##	1450000000.0, 1400000000.0, 55000000.0, 170000000.0, 90000000.0, 0.0, ), ))
WingModel.materials['CFRP-UD-Box'].elastic.FailStress(table=((1417500000.0, 1147500000.0, 50000000.0, 250000000.0, 70000000.0, 0.0, ), ))
WingModel.materials['CFRP-UD-Box'].elastic.FailStrain(table=(( 0.0105, 0.0085, 0.005, 0.025, 0.014), ))

# path_xfoilProgram = os.path.abspath('C:/Temp/user_element/xfoil6.96/bin/xfoil.exe')
path_xfoilProgram = os.path.abspath(xfoil_path)
path_xfoilFiles = os.path.join(folder_path,'XFOIL_files')
path_save = os.path.join(folder_path+'/Results-'+str(iSIM)+'_'+str(elements_inWeb) ,'Model-SparAngle-'+str(int(jobNumber)))

#f = open(folder_path+'/AeroCoeff/coeffs'+'-'+str(jobNumber)+'-'+str(int(1000*PositionFrontSpar)) + '-' + str(1000*int(AngleFrontSpar)) + '-' + str(int(1000*PositionRearSpar)) + '-' + str(int(1000*AngleRearSpar))+'-'+str(int(1000*thickness_bucklingSpar[digge]))+'-'+str(int(FaserWinkel[faserw]))+'.txt','w')
#f.close()
#f = open(folder_path+'/AeroCoeff/coeffs_Out_'+'-'+str(jobNumber)+'-'+str(int(1000*PositionFrontSpar)) + '-' + str(1000*int(AngleFrontSpar)) + '-' + str(int(1000*PositionRearSpar)) + '-' + str(int(1000*AngleRearSpar))+'-'+str(int(1000*thickness_bucklingSpar[digge]))+'-'+str(int(FaserWinkel[faserw]))+'.txt','w')
#f.close()  
#f = open(folder_path+'/AeroCoeff/p'+'-'+str(jobNumber)+'-'+str(int(1000*PositionFrontSpar)) + '-' + str(1000*int(AngleFrontSpar)) + '-' + str(int(1000*PositionRearSpar)) + '-' + str(int(1000*AngleRearSpar))+'-'+str(int(1000*thickness_bucklingSpar[digge]))+'-'+str(int(FaserWinkel[faserw]))+'.txt','w')
#f.close()
f = open(folder_path+'/resultsOUT/r'+'-'+str(iSIM)+'_'+str(elements_inWeb) +'.txt','w')
f.close()
 
print('BEGINNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN 1')


#######################################################################
# FLIP NORMALS ON WINGBOX
# Commented out, Alejandro Valverde
setWINGBOXoben = model.parts['Wing-final'].sets['box_set_oben']
setWINGBOXunten = model.parts['Wing-final'].sets['box_set_unten']

for item in setWINGBOXoben.faces:
	model.parts['Wing-final'].flipNormal(regions=Region(
		faces=model.parts['Wing-final'].faces.findAt(coordinates=(item.pointOn[0],))
		))

for item in setWINGBOXunten.faces:
	model.parts['Wing-final'].flipNormal(regions=Region(
		faces=model.parts['Wing-final'].faces.findAt(coordinates=(item.pointOn[0],))
		))

#######################################################################
## B.C.


#######################################################################
## B.C. # Already defined 

# Encastre
# model.rootAssembly.Set(edges=model.rootAssembly.instances['Wing-final-1'].edges.findAt(
#     ((points_nose[0][0], 0.0001, 0.0), ),((points_box_up[-1][0], points_box_up[-1][1],  0.0), ), ((points_box_up[1][0], points_box_up[1][1],  0.0), ), ((points_nose[0][0], -0.0001, 0.0), ), ((rib_trail_low[0][0], 0.0001, 0.0), ), 
#     ((points_box_low[1][0], points_box_low[1][1],  0.0), ), ((rib_trail_low[0][0], -0.0001, 0.0), ), ), name='All_BC')
# ##model.EncastreBC(createStepName='Initial', localCsys=None, name='BC-1', region=model.rootAssembly.sets['Set-265'])
# model.EncastreBC(createStepName='Initial', localCsys=None, name='All_BC', region=model.rootAssembly.sets['All_BC'])

# # Boundary condition with degree of freedom in spanwise direction (sliding allowed)
# edges_BC1 = model.rootAssembly.instances['Wing-final-1'].edges.findAt(((rib_trail_up[-2][0], rib_trail_up[-2][1], 0.0) ,),)
# edges_BC2 = model.rootAssembly.instances['Wing-final-1'].edges.findAt(((rib_trail_low[1][0], rib_trail_low[1][1], 0.0) ,),)
# edges_BC4 = model.rootAssembly.instances['Wing-final-1'].edges.findAt(((points_nose[-2][0], points_nose[-2][1], 0.0) ,),)
# edges_BC5 = model.rootAssembly.instances['Wing-final-1'].edges.findAt(((points_nose[1][0], points_nose[1][1], 0.0) ,),)
# edges_BC6 = model.rootAssembly.instances['Wing-final-1'].edges.findAt(((points_nose[len(points_nose)/2+1][0], points_nose[len(points_nose)/2+1][1], 0.0) ,),)
# edges_BC3 = model.rootAssembly.instances['Wing-final-1'].edges.findAt(((points_trail3_low[-1][0], points_trail3_low[-1][1], 0.0) ,),)
# edges_BC7 = model.rootAssembly.instances['Wing-final-1'].edges.findAt(((points_trail4_low[-1][0], points_trail4_low[-1][1], 0.0) ,),)

# model.rootAssembly.Set(name='Edges_BC', edges=edges_BC1+edges_BC2+edges_BC3+edges_BC4+edges_BC5+edges_BC6+edges_BC7)

##model.DisplacementBC(name='BC-2', createStepName='Initial', 
##    region=model.rootAssembly.sets['Edges_BC'], u1=SET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=SET, 
##    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
##
##mdb.models['Model-SparAngle-1'].rootAssembly.SetByBoolean(name='ALL_BC', sets=(
##    mdb.models['Model-SparAngle-1'].rootAssembly.sets['Set-265'], 
##    mdb.models['Model-SparAngle-1'].rootAssembly.sets['Edges_BC']))
##
########
##
##mdb.models['Model-SparAngle-1'].rootAssembly.regenerate()
##mdb.models['Model-SparAngle-1'].rootAssembly.ReferencePoint(point=
##    mdb.models['Model-SparAngle-1'].rootAssembly.instances['Wing-final-1'].vertices[270])
##mdb.models['Model-SparAngle-1'].rootAssembly.Set(name='AAAAA', referencePoints=
##    (mdb.models['Model-SparAngle-1'].rootAssembly.referencePoints[144], ))
##mdb.models['Model-SparAngle-1'].Coupling(controlPoint=
##    mdb.models['Model-SparAngle-1'].rootAssembly.sets['AAAAA'], couplingType=
##    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
##    'Constraint-1', surface=
##    mdb.models['Model-SparAngle-1'].rootAssembly.sets['ALL_BC'], u1=ON, u2=ON, 
##    u3=ON, ur1=ON, ur2=ON, ur3=ON)
##mdb.models['Model-SparAngle-1'].boundaryConditions['BC-2'].suppress()
##mdb.models['Model-SparAngle-1'].boundaryConditions['BC-1'].suppress()
##mdb.models['Model-SparAngle-1'].DisplacementBC(amplitude=UNSET, createStepName=
##    'Initial', distributionType=UNIFORM, fieldName='', localCsys=None, name=
##    'BC-3', region=mdb.models['Model-SparAngle-1'].rootAssembly.sets['AAAAA'], 
##    u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
##mdb.models['Model-SparAngle-1'].rootAssembly.Set(edges=
##    mdb.models['Model-SparAngle-1'].rootAssembly.instances['Wing-final-1'].edges.getSequenceFromMask(
##    ('[#0:10 #40 #0 #80000000 :14 #18400000 #100 ]', ), ), name='Set-117')
##
##
#########################

#######################################################################
# aerodynamic profile: points from file
CS = []
index = 1
#with open(start_dir+'/NACA0015.dat') as f:
#with open(start_dir+'/fteroFOIL150200.dat') as f:
with open(start_dir+'/NACA0012_199points.dat') as f:
    next(f)
    for line in f:
        CS.append([float(x) for x in line.split()])

f.close()

CS = [[-chord*a, chord*b] for [a,b] in CS]
CS = [[-a, b] for [a,b] in CS]


print('BEGINNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN 3')

##-----------------------------------------------------
## coordinates for aerodynamic pressure distribution for undeformed
## normalize to unit chord
#CS_norm = [map(lambda x: x/chord, p) for p in CS]
## write *.dat file with normalized coordinates of undeformed aerofoil for XFOIL input
#coord = file(os.path.join(path_xfoilFiles,'coords_'+str(jobNumber)+'-'+str(index-1)+'-'+str(1)+'.dat'),'w')
#coord.write(jobname+'     beginning of '+'Step-1'+'\n')
#for line in CS_norm:
#    coord.write('     '.join(str(round(x,7)) for x in line) + '\n')
#
#coord.close()


print('BEGINNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN 4')

## coordinates of aerofoil nodes for Step-1/undeformed
#aerofoilCoords = []
#aerofoilCoordsDef = []
#for i in range(N):
#	aerofoilCoords.append([])
#	aerofoilCoordsDef.append([])
#	aerofoilCoords[i] = sectionCoordsInitial(instance)

# coordinates of cambered aerofoil nodes for Step-1/undeformed
#camberline = np.genfromtxt(start_dir+'/fteroCL.dat', delimiter=None, skip_header=1)
#xCL = camberline[:, 0]
#yCL = camberline[:, 1]
xCL = np.linspace(0,chord,100)
yCL = [0]*100
aerofoilCoords = []
aerofoilCoordsDef = []
for i in range(N):
	aerofoilCoords.append([])
	aerofoilCoordsDef.append([])
	aerofoilCoords[i] = sectionCoordsInitialCAMBERED(instance,xCL,yCL)

#airfoilX = []
#airfoilY = []	
#for i in range(len(aerofoilCoords[0])):
#	airfoilX.append(aerofoilCoords[0][i][0])
#	airfoilY.append(aerofoilCoords[0][i][1])


# aeroelastic loop
loopNoConvergence = True
index = 1
	
# for index in range(1,int(t_tot/t_step)+1):
# for index in range(1,3):

print('haalllooooooooooooooooooooooooooooooooooooooooooooooooooooooooo')

print(os.getcwd())
twistVORHER=0

#alphaPOLAR
#clPOLAR
#cdPOLAR
#alphaPOLARcp
#topCpPRE
#bottomCpPRE
#execfile(start_dir+'\LoadPREDEFINEDaero.py')
#execfile(start_dir+'/LoadPREDEFINEDaeroREFINED.py')
#execfile(start_dir+'\LoadPREDEFINEDaeroNACA0012.py')
execfile(start_dir+'\loadPREDEFINEDaeroNACA0012_2.py')

#if index == 3:

displTEold=0

while loopNoConvergence:
	#f = open(folder_path+'/testW'+'-'+str(index) +'.txt','w')
	#f.close()
        os.remove(folder_path+'/Geklappt'+str(iSIM)+'_'+str(elements_inWeb)+'.txt')
        geklappt = open(folder_path+'/Geklappt'+str(iSIM)+'_'+str(elements_inWeb)+'.txt','a')
        geklappt.write('No'+'\n')
        geklappt.close()
        
	print(loopNoConvergence)
	primaryJobName = jobname + '_' + str(jobNumber) + '-' + str(index)
        
	if 	index == 1:# and alpha == 0):
		dpos = [1]*N
		lldisc = 1
		alphaPOLARt = [alphaPOLAR]*N
		alphaPOLARcpt = [alphaPOLARcp]*N
	else:
		dpos = [ x+1 for x in range(N)]
		lldisc = N
		#N = 5
		alphaPOLARt = []
		alphaPOLARcpt = []
		for i in range(N):
			alphaPOLARt.append([])
			alphaPOLARt[i][:] = [x + twist[i] for x in alphaPOLAR]
			alphaPOLARcpt.append([])
			alphaPOLARcpt[i][:] = [x + twist[i] for x in alphaPOLARcp]			
	print(1)	
	#------------------------------------------------------------------
	# WEISSINGER	
	liftslope = [0]*N
	alpha0 = [0]*N
	aeroCoeffsLL1 = [0]*N
	aeroCoeffsLL2 = [0]*N 
	#------------------------------------------------------------------
	# calculate linear polar of each section
	for dp in range(lldisc):
		# calculate parameters for linear lifting line: coeff 2
		#alphaLL2 = alpha + 0.99 #0.99
		#flagLL = True
		#while flagLL:
		#	alphaLL2 = alphaLL2 + 0.01
		#	flagLL = getPressureXfoilSpan(path_xfoilFiles,jobNumber,index,path_xfoilProgram,Re,Ma,alphaLL2,dpos[dp])
		#aeroCoeffsLL2[dp] = extractCoeffs(path_xfoilFiles, jobNumber, index, dpos[dp])
		#alphaLL2 = 2.0#alpha + 1.0
		alphaLL2 = alpha + 1.0		
		aeroCoeffsLL2[dp] = [numpy.interp(alphaLL2, alphaPOLARt[dp], clPOLAR), numpy.interp(alphaLL2, alphaPOLARt[dp], cdPOLAR)]
                
		# calculate parameters for linear lifting line: coeff 1
		#alphaLL1 = alpha - 0.01 # 0.01
		#flagNoConvXfoil = True
		#while flagNoConvXfoil:
		#	alphaLL1 = alphaLL1 + 0.01
		#	flagNoConvXfoil = getPressureXfoilSpan(path_xfoilFiles,jobNumber,index,path_xfoilProgram,Re,Ma,alphaLL1,dpos[dp])
		#aeroCoeffsLL1[dp] = extractCoeffs(path_xfoilFiles, jobNumber, index, dpos[dp])
		#alphaLL1 = 0.0#alpha # alpha0[0]*180/pi
		alphaLL1 = alpha # alpha0[0]*180/pi		
		aeroCoeffsLL1[dp] = [numpy.interp(alphaLL1, alphaPOLARt[dp], clPOLAR), numpy.interp(alphaLL1, alphaPOLARt[dp], cdPOLAR)]
		
		# calculate parameters for lifting line: lift slope and alpha0
		liftslope[dp] = (aeroCoeffsLL2[dp][0] - aeroCoeffsLL1[dp][0])/((alphaLL2-alphaLL1)*pi/180) 
		alpha0[dp] = (aeroCoeffsLL2[dp][0]*(alphaLL1*pi/180)-aeroCoeffsLL1[dp][0]*(alphaLL2*pi/180))/(aeroCoeffsLL2[dp][0] - aeroCoeffsLL1[dp][0])
	# calculate lifting line
        
	if index == 1:# and alpha == 0):
		for i in range(N-1):
			alpha0[i+1] = alpha0[0]
			liftslope[i+1] = liftslope[0]
			aeroCoeffsLL1[i+1] = aeroCoeffsLL1[0]
			aeroCoeffsLL2[i+1] = aeroCoeffsLL2[0]
	#------------------------------------------------------------------
	# linear weissinger
	alphaefflin, w_ij_trailing, gamma = weissinger(chord,span,alpha,alpha0,liftslope,2*N)
	print(2)
        
	#------------------------------------------------------------------
	# nonlinear weissinger
	#index = 10
	alphatol = 0.001
	alphaNL = 0 # alpha nonlinear
	# if (index != 1 and alpha >= alphaNL):  #min(alpha0) != 0
	if alpha >= alphaNL:  #min(alpha0) != 0
		#alphamin = int(min(alpha0))-1
		alphamin = int(min(alpha0)*180/pi)-5
		alphamax = alpha + 1
		alphastep = 0.1#0.5
		alen = int((alphamax-alphamin)/alphastep + 1)
		clTable = [[0 for x in range(alen)] for x in range(N)]
		alphaTable = [[0 for x in range(alen)] for x in range(N)]
		#clTable = [[0 for x in range(len(clPOLAR))] for x in range(N)]
		#alphaTable = [[0 for x in range(len(clPOLAR))] for x in range(N)]
		# Call XFOIL to generate Polar
		for j in range(N):
			## Generate Polar from alpha = alphamin
			#flagLL = getPolarXfoil(path_xfoilFiles,jobNumber,index,path_xfoilProgram,Re,Ma,dpos[j],alphamin,alphamax,alphastep)
			#aeroCoeffs = extractCoeffsPolar(path_xfoilFiles, jobNumber, index, dpos[j])
			#for i in range(len(alphaTable[j])):
			#	if i < len(aeroCoeffs):
			#		alphaTable[j][i] = aeroCoeffs[i][0]
			#		clTable[j][i] = aeroCoeffs[i][1]
			#	else:
			#		alphaTable[j][i] = aeroCoeffs[len(aeroCoeffs)-1][0]
			#		clTable[j][i] = aeroCoeffs[len(aeroCoeffs)-1][1]
			for i in range(len(alphaTable[j])):
				alphaTable[j][i] = alphamin+i*alphastep
				clTable[j][i] = numpy.interp(alphamin+i*alphastep,alphaPOLARt[j],clPOLAR)		
			#alphaTable[j] = alphaPOLARt[j]
			#clTable[j] = clPOLAR
##			
		# Calculate Non Linear Weissinger   for x in range(len(alphaefflin)): alphaefflin[x]*180/pi
		#f = open(folder_path+'/testW'+'-'+str(index) +'.txt','a')
		#f.write(str(alphaefflin)+'\n'+str(w_ij_trailing)+'\n'+str(gamma)+'\n'+str(alphaTable)+'\n'+str(clTable)+'\n'+str(chord)+'\n'+str(alpha)+'\n'+str(alpha0)+'\n'+str(2*N)+'\n')
		#f.close()
#		alphaeff = nlweissinger(alphaefflin, w_ij_trailing, gamma, alphaTable, clTable, chord, alpha, alpha0, 2*N)
		alphaeff = nlweissingerOPT(alphaefflin, w_ij_trailing, gamma, alphaTable, clTable, chord, alpha, alpha0, 2*N)
		print(3)
		
	else:
		alphaeff = [0]*(N)
		for i in range(N):
			alphaeff[i] = alphaefflin[N+i]
			
	# get pressure from xfoil with alpha effective from weissinger
	topCp = []
	bottomCp = []
	aeroCoeffs = [0]*N
	clWing = [0]*N
	cdWing = [0]*N
	clNew = [0]*N
	cdNew = [0]*N
	alphatol = 0.001
        
	for i in range(len(alphaeff)):
	    topCp.append([])
	    bottomCp.append([])
	    #alphall = alphaeff[i]*180/pi - alphatol
	    #flagLL = True
	    #while flagLL:
	    #	alphall = alphall + alphatol
	    #	flagLL = getPressureXfoilSpan(path_xfoilFiles,jobNumber,index,path_xfoilProgram,Re,Ma,alphall,dpos[i])
        #
	    #topCp[i], bottomCp[i] = dataAeroCpRoot(path_xfoilFiles,jobNumber,index,chord, dpos[i])
        #
	    #aeroCoeffs[i] = extractCoeffs(path_xfoilFiles, jobNumber, index, dpos[i])
	    
	    testALPHA = 0
	    posALPHA = 0
	    alphall = alphaeff[i]*180/pi
	    print(alphall)
	    while testALPHA == 0:
                			
			#if alphaPOLAR[posALPHA] < alphall:
			if alphaPOLARcpt[i][posALPHA] < alphall:#alphaPOLAR[posALPHA] < alphall:
				posALPHA = posALPHA+1
			else:
				testALPHA = 1
	    
	    #numpy.interp(alphall, [alphaPOLAR[posALPHA-1],alphaPOLAR[posALPHA]], [clPOLAR[posALPHA-1],clPOLAR[posALPHA]])
            	    
	    topCpA = []
	    bottomCpA = []
	    for j in range(len(topCpPRE[0])):
#	    	topCpA.append([topCpPRE[posALPHA-1][j][0], numpy.interp(alphall, [alphaPOLAR[posALPHA-1],alphaPOLAR[posALPHA]], [topCpPRE[posALPHA-1][j][1],topCpPRE[posALPHA][j][1]])])
#	    	bottomCpA.append([bottomCpPRE[posALPHA-1][j][0], numpy.interp(alphall, [alphaPOLAR[posALPHA-1],alphaPOLAR[posALPHA]], [bottomCpPRE[posALPHA-1][j][1],bottomCpPRE[posALPHA][j][1]])])
	    	topCpA.append([topCpPRE[posALPHA-1][j][0], numpy.interp(alphall, [alphaPOLARcpt[i][posALPHA-1],alphaPOLARcpt[i][posALPHA]], [topCpPRE[posALPHA-1][j][1],topCpPRE[posALPHA][j][1]])])
	    	bottomCpA.append([bottomCpPRE[posALPHA-1][j][0], numpy.interp(alphall, [alphaPOLARcpt[i][posALPHA-1],alphaPOLARcpt[i][posALPHA]], [bottomCpPRE[posALPHA-1][j][1],bottomCpPRE[posALPHA][j][1]])])
                 
	    topCp[i] = topCpA
	    bottomCp[i] = bottomCpA
                
	    aeroCoeffs[i] = [numpy.interp(alphall, alphaPOLARt[i], clPOLAR), numpy.interp(alphall, alphaPOLARt[i], cdPOLAR)]
        
	    
	    clWing[i] = aeroCoeffs[i][0]
	    cdWing[i] = aeroCoeffs[i][1]
	    clNew[i] = aeroCoeffs[i][0]*cos(alpha/180.0*pi-alphaeff[i])-aeroCoeffs[i][1]*sin(alpha/180.0*pi-alphaeff[i])
	    cdNew[i] = aeroCoeffs[i][1]*cos(alpha/180.0*pi-alphaeff[i])+aeroCoeffs[i][0]*sin(alpha/180.0*pi-alphaeff[i])
	    
	    ## f = open(folder_path+'/AeroCoeff/coeffs.txt','a')  
	    ## f.write('c_l = '+ str(aeroCoeffs[i][0])+ '\t')
	    ## f.write('c_d = '+ str(aeroCoeffs[i][1])+ '\t')
	    ## f.write('alphaeff = '+ str(alphaeff[i]*180.0/pi)+ '\t')
	    ## f.write('alpha = '+ str(alpha)+ '\n')		
	    ## f.close()		
	    #f = open(folder_path+'/AeroCoeff/coeffs'+'-'+str(jobNumber)+'-'+str(int(1000*PositionFrontSpar)) + '-' + str(1000*int(AngleFrontSpar)) + '-' + str(int(1000*PositionRearSpar)) + '-' + str(int(1000*AngleRearSpar))+'-'+str(int(1000*thickness_bucklingSpar[digge]))+'-'+str(int(FaserWinkel[faserw]))+'.txt','a')  
	    #f.write('c_l = '+ str(clNew[i])+ '\t')
	    #f.write('c_d = '+ str(cdNew[i])+ '\t')
	    #f.write('alphaeff = '+ str(alphaeff[i]*180.0/pi)+ '\t')
	    #f.write('alpha = '+ str(alpha)+ '\n')		
	    #f.close()
	
##	f1 = open(folder_path+'/AeroCoeff/coeffs_Out_'+'-'+str(jobNumber)+'-'+str(int(1000*PositionFrontSpar)) + '-' + str(1000*int(AngleFrontSpar)) + '-' + str(int(1000*PositionRearSpar)) + '-' + str(int(1000*AngleRearSpar))+'-'+str(int(1000*thickness_bucklingSpar[digge]))+'-'+str(int(FaserWinkel[faserw]))+'.txt','a')  
	f1 = open(folder_path+'/AeroCoeff/coeffs_Out_.txt','a')  
	f1.write('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\n' + 'Cl = '+ str(clNew)+'\n' + 'Cd = '+ str(cdNew) +'\n')
	f1.close()
	    
	clWingMean = numpy.mean(clWing)
	cdWingMean = numpy.mean(cdWing)
	clcdWing = clWingMean/cdWingMean
	clNewMean = numpy.mean(clNew)
	cdNewMean = numpy.mean(cdNew)
	clcdNew = clNewMean/cdNewMean
                
        hebelarm=y_C[len(y_C)/2.0:len(y_C)]
        cm = y_C
        Mb = [0.0]*N
        Area_Hebelarm = [0.06*0.5]*N
        for i in range(len(hebelarm)):
                Mb[i]=hebelarm[i]*clNew[i]*Area_Hebelarm[i]*0.5*rho*(sqrt(t_step)*V)*(sqrt(t_step)*V)

        Mb_Ges=sum(Mb)
	f = open(folder_path+'/AeroCoeff/coeffWing'+'-'+str(iSIM) +'_'+str(elements_inWeb) + '-' +str(jobNumber)+'-'+ str(index) +'.txt','a')
	# f.wrte('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\t' + 'cl_wing = '+ str(clWingMean)+'\t' + 'cd_wing = '+ str(cdWingMean) +'\t' + 'cl/cd = '+ str(clcdWing)+'\n')
##	f.write('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\t' + 'cl_wing = '+ str(clNewMean)+'\t' + 'cd_wing = '+ str(cdNewMean) +'\t' + 'cmd = '+ str(Mb_Ges)+'\n')
	f.write(str(jobNumber)+ '\t' + str(index) + '\t' + str(clNewMean)+'\t' +  str(cdNewMean) +'\t' + str(Mb_Ges)+'\t' + str(t_step)+'\n')
	f.close()
        
        print('Nach Weissinger:')
        print(os.getcwd())
	# calculate imperfection
	if index == 1:
		model.BuckleStep(name='Step-EV', previous='Initial', numEigen=10)
		model.steps['Step-EV'].setValues(vectors=100, maxIterations=1000)
		model.steps['Step-EV'].setValues(blockSize=DEFAULT, eigensolver=LANCZOS, maxBlocks=DEFAULT, minEigen=0.0)
                
#		pressureApplyEV(topCp,bottomCp,aerofoilCoords,span,index,aeroQ, model, assembly, instance,N,chord,folder_path,t_step)
#		pressureApplyEV(topCp,bottomCp,aerofoilCoords,span/2,index,aeroQ, model, assembly, instance,N,chord,folder_path,t_step)
		pressureApplyEVcamber(topCp,bottomCp,aerofoilCoords,span/2,index,aeroQ, model, assembly, instance,N,chord,folder_path,t_step,xCL,yCL)
		
		mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
			explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
			memory=90, memoryUnits=PERCENTAGE, model='Model-SparAngle-'+str(int(jobNumber)), modelPrint=OFF, 
			multiprocessingMode=DEFAULT, name='Wing-EV-'+str(int(jobNumber)), nodalOutputPrecision=SINGLE, 
			numCpus=1, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
			userSubroutine='', waitHours=0, waitMinutes=0)
			
		model.keywordBlock.synchVersions(storeNodesAndElements=False)
		n_kw=len(model.keywordBlock.sieBlocks)
		model.keywordBlock.insert(n_kw-2, '*NODE FILE, LAST MODE=10, GLOBAL=YES\nU')
		mdb.jobs['Wing-EV-'+str(int(jobNumber))].setValues(numCpus=1, numDomains=1)
		mdb.jobs['Wing-EV-'+str(int(jobNumber))].submit(consistencyChecking=OFF)
		mdb.jobs['Wing-EV-'+str(int(jobNumber))].waitForCompletion()
		
		model.keywordBlock.synchVersions(storeNodesAndElements=False)
		n_replace = GetBlockPosition('Model-SparAngle-'+str(int(jobNumber)),'*NODE FILE, LAST MODE=')
		model.keywordBlock.replace(n_replace, """ """)
		model.keywordBlock.replace(n_replace, """ """)
		
		del model.steps['Step-EV']
		del assembly.surfaces['pSurf-0']
		
		print('Nach EV:')
		execfile(folder_path+'/rootBendingMoment/evODB_Falk.py') 
		# define new step
		evFile = open(folder_path+'/resultsOUT/r'+'-'+str(iSIM)+'_'+str(elements_inWeb) +'.txt','a')
                evFile.write(str(minEV) +'\n')   
                evFile.close()
                print('VOOOOOOOOOOOOOORRRRRR NLLLL:')

	aeroLoopStepsNL(jobname,jobNumber,index,model,dampCoeff_Step1,initialIncrement,varsFieldOut,varsHistoryOut,dampMagnit)
	#------------------------------------------------------------------
	# apply pressure by mapped field
	if index == 1:
		roin = []
		pNorm = []
#	roin, pNorm = pressureApplyFalk(topCp,bottomCp,aerofoilCoords,span,index,aeroQ, model, assembly, instance,N,chord,roin,pNorm,folder_path,t_step)
#	roin, pNorm = pressureApplyFalk(topCp,bottomCp,aerofoilCoords,span/2,index,aeroQ, model, assembly, instance,N,chord,roin,pNorm,folder_path,t_step)

##	if index ==1:
##		t_step=0.1#minEV # fuer lineare extrapolation
	roin, pNorm = pressureApplyFalkcamber(topCp,bottomCp,aerofoilCoords,span/2,index,aeroQ, model, assembly, instance,N,chord,roin,pNorm,folder_path,t_step,xCL,yCL)
	
	# job
	if index ==1:
			jobType = ANALYSIS
	else:
			jobType = RESTART

        # jobLoop = mdb.Job(name=primaryJobName, 
				  # model='Model-1', 
				  # userSubroutine='',
				  # type=jobType, 
				  # explicitPrecision=SINGLE,
				  # nodalOutputPrecision=SINGLE, description='',
				  # parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT,
				  # numDomains=6, numCpus=6, scratch='',
				  # echoPrint=OFF, modelPrint=ON, contactPrint=OFF, historyPrint=OFF)
        print('Vor Job-Definition:')
        print(os.getcwd())
        jobLoop = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
                                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
                                memory=90, memoryUnits=PERCENTAGE, model='Model-SparAngle-'+str(int(jobNumber)), modelPrint=OFF, 
                                multiprocessingMode=DEFAULT, name=primaryJobName, nodalOutputPrecision=SINGLE, 
                                numCpus=4, numDomains=4, numGPUs=0, queue=None, scratch='', type=jobType, 
                                userSubroutine='', waitHours=0, waitMinutes=0)				  

        if index == 1:
                model.keywordBlock.synchVersions(storeNodesAndElements=False)
                pos_step=0
                for block in model.keywordBlock.sieBlocks:
                        if string.lower(block[0:len('*Step')])==string.lower('*Step'):
                                break
                        else:
                                pos_step=pos_step+1		
                
                model.keywordBlock.insert(pos_step-1,'*IMPERFECTION, FILE=Wing-EV-'+str(int(jobNumber))+', STEP=1\n1,5e-7\n2,3e-7\n3,1e-7\n4,1e-8\n5,1e-8\n6,1e-7\n7,1e-8')
                
                
                model.keywordBlock.synchVersions(storeNodesAndElements=False)
                delConflict = True
                while delConflict:
                        try:
                                n_conflict = GetBlockPosition('Model-SparAngle-'+str(int(jobNumber)),'*Conflict')
                                model.keywordBlock.replace(n_conflict, """ """)
                        except:
                                delConflict = False				
                                

        jobLoop.writeInput(consistencyChecking=OFF)
        
	#-----------------------------------------------------
	# save completed model
	mdb.saveAs(pathName=path_save)
        
	#-----------------------------------------------------
	# run job
	# print('Exiting script') #Manual exit, Alejandro Valverde
	# exit()
	jobLoop.submit()
	jobLoop.waitForCompletion()
        
	os.remove(folder_path+'/Geklappt'+str(iSIM)+'_'+str(elements_inWeb)+'.txt')
        geklappt = open(folder_path+'/Geklappt'+str(iSIM)+'_'+str(elements_inWeb)+'.txt','a')
        geklappt.write('Yes'+'\n')
        geklappt.close()
        
        for i in range(N):
                #twist[i], displLE[i], displTE[i], aerofoilCoordsDef[i] = coordsDeformedAerofoil(path_xfoilFiles,primaryJobName,jobNumber,index,chord,i+1)
                twist[i], displLE[i], displTE[i], aerofoilCoordsDef[i] = coordsDeformedAerofoilCAMBER(path_xfoilFiles,primaryJobName,jobNumber,index,chord,i+1,xCL,yCL)
        
        ## loop convergence criterion
        #if index > 1:
        #        deltaDisplTE = abs(displTEold - displTE[len(displTE)-1])
        #        if (index*t_step >= 1 and abs(abs(twist[-1])-abs(twistVORHER)) < 0.01):    # 0.001
        #                loopNoConvergence = False
        # loop convergence criterion
        #maxSTRESSout = maxSTRESS(primaryJobName,jobNumber)
        #print(maxSTRESSout)
        #print('maaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaxstreeeeeeeeeeeeeeeeeeeeeeeessssss')
        t_step_readout=t_step
        if index >= 150:
                loopNoConvergence = False
        
        if index >= 1:
                deltaDisplTE = abs(displTEold - displTE[len(displTE)-1])
                if t_step == 1.0 and abs(abs(twist[-1])-abs(twistVORHER)) < 0.01:
                        loopNoConvergence = False
                        maxSTRESSout = maxSTRESS(primaryJobName,index)
                        f = open(folder_path+'/resultsOUT/r'+'-'+str(iSIM)+'_'+str(elements_inWeb) +'.txt','a')
                        f.write(str(twist[-1])+'\n'+str(maxSTRESSout)+'\n')
                        f.close()
                
                if t_step<1.0 and abs(abs(twist[-1])-abs(twistVORHER)) < 0.01:
                        t_step = t_step+t_step_0 #1.0
        
##                else:
##                        t_step = 1.0
##						if t_step<1.0:
##                                                        t_step = t_step+0.1#1.0
##                                                
##                                                afterBUCKLING = 1
##						f = open(folder_path+'/resultsOUT/r'+'-'+str(iSIM) +'.txt','a')
##						f.write(str(minEV)+'\n'+str(twist[-1])+'\n')
##						f.close()
                
##						
						
        TwistDifference=abs(abs(twist[-1])-abs(twistVORHER))
        displTEold = displTE[len(displTE)-1]
        twistVORHER=twist[-1]
        
        # save twist to text file 
        f = open(folder_path+'/AeroCoeff/twistWing'+'-'+str(iSIM)+'_'+str(elements_inWeb)+'.txt','a')
        f.write('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\n' + 'twist = '+ str(twist)+'\n' + 'displLE = '+ str(displLE) +'\n' + 'displTE = '+ str(displTE)+'\n')
        f.close()
        
        
        # save pNorm for post processing
        #if index == 1:
        #        f = open(folder_path+'/AeroCoeff/pNorm.txt','a')  
         #       f.write(str(pNorm))
          #      f.close()
          
        # execfile(folder_path+'/rootBendingMoment/pressureODB_Falk_OnlyLastFrame.py')
        execfile(folder_path+'/rootBendingMoment/pressureODB_Falk_OnlyLastFrame_sim.py')
        if loopNoConvergence == False:
                execfile(folder_path+'/rootBendingMoment/pressureODB_AverageStress.py')

        ##        execfile(folder_path+'/rootBendingMoment/pressureODB_Falk_OnlyTwist.py')        
        
        
        print('GUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDNGUDN')

        if index>1:
                try:
                        primaryJobName = 'Wing-stat_'+str(jobNumber)+'-' + str(index-1)
                        sttName = primaryJobName + '.stt'
                        odbName = primaryJobName + '.odb'
                        mdlName = primaryJobName + '.mdl'
                        os.remove(folder_path+'/Results-'+str(iSIM)+'_'+str(elements_inWeb)+'/'+sttName)
                        # os.remove(folder_path+'/Results-'+str(iSIM)+'_'+str(elements_inWeb)+'/'+odbName) #commented out, Alejandro Valverde
                        os.remove(folder_path+'/Results-'+str(iSIM)+'_'+str(elements_inWeb)+'/'+mdlName)
                except:
                        fehlerbeimloeschen=1
        
        index = index + 1

print(span)


##for fileName in range(1,index-1):
##        primaryJobName = 'Wing-stat_'+str(jobNumber)+'-' + str(fileName)
##        sttName = primaryJobName + '.stt'
##        odbName = primaryJobName + '.odb'
##        mdlName = primaryJobName + '.mdl'
##        os.remove(folder_path+'/Results-'+str(iSIM)+'/'+sttName)
##        os.remove(folder_path+'/Results-'+str(iSIM)+'/'+odbName)
##        os.remove(folder_path+'/Results-'+str(iSIM)+'/'+mdlName)


# # remove xfoil files
# fileList = os.listdir(path_xfoilFiles)
# for fileName in fileList:
	# os.remove(path_xfoilFiles+"/"+fileName)			
