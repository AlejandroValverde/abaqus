#===============================================================================
# Functions Urban
#===============================================================================


from abaqus import *
from abaqusConstants import *
import visualization
import os
import re
import sys
import shutil
import subprocess as sp
import numpy
import math
import regionToolset
import copy
# import numpy as np
# from scipy import interpolate
from parameters_Dom_sim import *

def coordsDeformedAerofoilFalk(odb,pathXfoilData,jobName,jobNumber,jobIndex,chordLength,dpos, iterat):
    # open output database
##    odbName = jobName + '.odb'
##    odb = visualization.openOdb(odbName)
    
    stepName = 'Step-%s' %jobIndex
        
    setName = 'SECTION%s' %dpos
    
    aerofoilY0 = odb.rootAssembly.instances['WING-FINAL-1'].nodeSets[setName]
        
    aerofoilCoordsX = []
    aerofoilCoordsZ = []
    
    for k in aerofoilY0.nodes:
        aerofoilCoordsX.append(k.coordinates[0])
        aerofoilCoordsZ.append(k.coordinates[1])
    
    displAerofoilX = []
    displAerofoilZ = []
    
    dispField = odb.steps[stepName].frames[iterat].fieldOutputs['U']
    dispSubField = dispField.getSubset(region=aerofoilY0)
    displValues = dispSubField.values
    
    for node in displValues:
            displAerofoilX.append(node.data[0])
            displAerofoilZ.append(node.data[1])
    
    ##odb.close()
    #-----------------------------------------------------
        
    aerofoilData = zip(aerofoilCoordsX,displAerofoilX,aerofoilCoordsZ,displAerofoilZ)
    
    aerofoilData_top = []
    aerofoilData_bottom = []
    
    for value in aerofoilData:
                if (not(1e-07 > value[2] > -1e-07 and value[0] > 0.2)):
                        if value[2] >= 0.0:
                                aerofoilData_top.append(value)
                        else:
                                aerofoilData_bottom.append(value)		
    
    
    aerofoilData_top.sort(key=lambda x: x[0], reverse=True)
    aerofoilData_bottom.sort(key=lambda x: x[0])
    
    #-----------------------------------------------------
    # find trailing edge displacement of bottom node of TE gap
    displTE = aerofoilData_bottom[-1][3]
    displTEx = aerofoilData_bottom[-1][1]
    
    # find leading edge displacement of bottom node of LE
    displLE = aerofoilData_top[-1][3]
    displLEx = aerofoilData_top[-1][1]
    
    # twist
    twist = numpy.arctan((displTE-displLE)/(chordLength-displLEx+displTEx))*180/pi
    
    return twist

	
def coordsDeformedAerofoilFalkCAMBER(odb,pathXfoilData,jobName,jobNumber,jobIndex,chordLength,dpos, iterat,xCL, yCL):
    # open output database
##    odbName = jobName + '.odb'
##    odb = visualization.openOdb(odbName)
    
    stepName = 'Step-%s' %jobIndex
        
    setName = 'SECTION%s' %dpos
    
    aerofoilY0 = odb.rootAssembly.instances['WING-FINAL-1'].nodeSets[setName]
        
    aerofoilCoordsX = []
    aerofoilCoordsZ = []
    
    for k in aerofoilY0.nodes:
        aerofoilCoordsX.append(k.coordinates[0])
        aerofoilCoordsZ.append(k.coordinates[1])
    
    displAerofoilX = []
    displAerofoilZ = []
    
    dispField = odb.steps[stepName].frames[iterat].fieldOutputs['U']
    dispSubField = dispField.getSubset(region=aerofoilY0)
    displValues = dispSubField.values
    
    for node in displValues:
            displAerofoilX.append(node.data[0])
            displAerofoilZ.append(node.data[1])
    
    ##odb.close()
    #-----------------------------------------------------
        
    aerofoilData = zip(aerofoilCoordsX,displAerofoilX,aerofoilCoordsZ,displAerofoilZ)
    
    aerofoilData_top = []
    aerofoilData_bottom = []
    
    for value in aerofoilData:
                if (not(1e-07 > value[2] > -1e-07 and value[0] > 0.2)):
                        if numpy.interp(value[0], xCL, yCL) <= value[2]:
                                aerofoilData_top.append(value)
                        else:
                                aerofoilData_bottom.append(value)		
    
    
    aerofoilData_top.sort(key=lambda x: x[0], reverse=True)
    aerofoilData_bottom.sort(key=lambda x: x[0])
    
    #-----------------------------------------------------
    # find trailing edge displacement of bottom node of TE gap
    displTE = aerofoilData_bottom[-1][3]
    displTEx = aerofoilData_bottom[-1][1]
    
    # find leading edge displacement of bottom node of LE
    displLE = aerofoilData_top[-1][3]
    displLEx = aerofoilData_top[-1][1]
    
    # twist
    twist = numpy.arctan((displTE-displLE)/(chordLength-displLEx+displTEx))*180/pi
    
    return twist	
	
#######################################################################
# Assign twist angles

def twistangles_for_output(curTwist, position, breiteSection):

    count=int(position/breiteSection)                # array nummer des linken linearen interpolationsstuetzstelle
    first_twist_position=count*breiteSection
    second_twist_position=first_twist_position+breiteSection
    first_twist=curTwist[count]
    second_twist=curTwist[count+1]
    twistReturn = (second_twist-first_twist)/(second_twist_position-first_twist_position) * (position-first_twist_position) + first_twist
    return twistReturn

	
#######################################################################

#######################################################################
# OBTAIN AEROFOIL COORDINATES FOR STEP-1/UNDEFORMED

def sectionCoordsInitial(instance):

    aerofoilCoordsX = []
    aerofoilCoordsY = []
    
    # for k in instance.sets['Wing-1.Section1'].nodes:
    for k in instance.sets['Section1'].nodes:
        aerofoilCoordsX.append(k.coordinates[0])
        aerofoilCoordsY.append(k.coordinates[1])
    
    aerofoilData = zip(aerofoilCoordsX,aerofoilCoordsY)
    
    aerofoilData_top = []
    aerofoilData_bottom = []
    
    # sort data from top TE, around LE and bottom TE
    for value in aerofoilData:
        if value[1] >= 0.0:
            aerofoilData_top.append(value)
        else:
            aerofoilData_bottom.append(value)
    
    aerofoilData_top.sort(key=lambda x: x[0], reverse=True)
    aerofoilData_bottom.sort(key=lambda x: x[0])
    
    aerofoilDataFull = []
    for item in aerofoilData_top:
        aerofoilDataFull.append(item)
    for item in aerofoilData_bottom:
        aerofoilDataFull.append(item)
            
    return aerofoilDataFull

	
#######################################################################
# OBTAIN CAMBERED AEROFOIL COORDINATES FOR STEP-1/UNDEFORMED

def sectionCoordsInitialCAMBERED(instance,xCL,yCL):

    aerofoilCoordsX = []
    aerofoilCoordsY = []
    
    # for k in instance.sets['Wing-1.Section1'].nodes:
    for k in instance.sets['Section1'].nodes:
    	aerofoilCoordsX.append(k.coordinates[0])
    	aerofoilCoordsY.append(k.coordinates[1])
    
    aerofoilData = zip(aerofoilCoordsX,aerofoilCoordsY)
    
    aerofoilData_top = []
    aerofoilData_bottom = []
    
    # sort data from top TE, around LE and bottom TE by using camberline
    #numpy.interp(aerofoilData[0][0], xCL, yCL), aerofoilData[0][1]
    for value in aerofoilData:
    	if numpy.interp(value[0], xCL, yCL) <= value[1]:
    		aerofoilData_top.append(value)
    	else:
    		aerofoilData_bottom.append(value)
    
    aerofoilData_top.sort(key=lambda x: x[0], reverse=True)
    aerofoilData_bottom.sort(key=lambda x: x[0])
    
    aerofoilDataFull = []
    for item in aerofoilData_top:
    	aerofoilDataFull.append(item)
    for item in aerofoilData_bottom:
    	aerofoilDataFull.append(item)
            
    return aerofoilDataFull	
	
#######################################################################
# XFOIL ANALYSIS

def getPressureXfoilSpan(pathXfoilFiles,jobNumber,jobIndex,pathXfoilProgram,Re,Ma,alpha,dpos):

    #-----------------------------------------------------
    def issueCmd(cmd,echo=True):
        ps.stdin.write(cmd+'\n')
        if echo:
            print cmd
    #-----------------------------------------------------
    
    iterNo = 1000
    filePressure=os.path.join(pathXfoilFiles,'pressure_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat') #added dpos
    fileCoords=os.path.join(pathXfoilFiles,'coords_'+str(jobNumber)+'-'+str(jobIndex-1)+'-'+str(dpos)+'.dat') #added dpos
    fileXfoilCoords=os.path.join(pathXfoilFiles,'coordsXfoil_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat')
    fileCoeffs = os.path.join(pathXfoilFiles,'coeffs_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat') #added dpos
    fileOut = os.path.join(pathXfoilFiles,'xfoilOut_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.txt')
    filePlot = os.path.join(pathXfoilFiles,'plot_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.ps')

    flagNoConvXfoil = False
    
    #-----------------------------------------------------

    # remove if files already exist
    listFile2Remove=[filePressure,fileXfoilCoords,fileCoeffs,fileOut,filePlot]
    try:
        for j in listFile2Remove:
            if os.path.exists(j)==True:
                os.remove(j)
    except OSError:
        pass

    #-----------------------------------------------------

    log = open(fileOut,'w+')

    ps = sp.Popen([pathXfoilProgram],
                  stdin=sp.PIPE,
                  stdout=log,
                  stderr=None)
    issueCmd('PLOP')
    issueCmd('G')
    issueCmd('')
    issueCmd('LOAD')         
    issueCmd(fileCoords)
    issueCmd('PANE') # regenerate paneling
    issueCmd('SAVE')
    issueCmd(fileXfoilCoords)
    issueCmd('OPER')
    issueCmd('ITER')
    issueCmd(str(iterNo)) 
    issueCmd('VISC')
    issueCmd(str(Re))
    issueCmd('MACH')
    issueCmd(str(Ma))
    issueCmd('PACC')
    issueCmd(fileCoeffs)
    issueCmd('\n')      
    issueCmd('OPER')
    issueCmd('ALFA')
    issueCmd(str(alpha))
    issueCmd('HARD')
    issueCmd('CPWR')
    issueCmd(filePressure)
    issueCmd('PACC')
    issueCmd('\n')
    issueCmd('QUIT')
    
    ps.wait()
        
    log.seek(0)
    if 'VISCAL:  Convergence failed' in log.read():
        flagNoConvXfoil = True
        os.remove(filePressure)
        os.remove(fileCoeffs)
    log.close()
	
    return flagNoConvXfoil	
	
#######################################################################
# XFOIL POLAR

def getPolarXfoil(pathXfoilFiles,jobNumber,jobIndex,pathXfoilProgram,Re,Ma,dpos,amin,amax,astep):

    #-----------------------------------------------------
    def issueCmd(cmd,echo=True):
        ps.stdin.write(cmd+'\n')
        if echo:
            print cmd
    #-----------------------------------------------------
    
    iterNo = 200 #1000
    filePressure=os.path.join(pathXfoilFiles,'pressure_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat') #added dpos
    fileCoords=os.path.join(pathXfoilFiles,'coords_'+str(jobNumber)+'-'+str(jobIndex-1)+'-'+str(dpos)+'.dat') #added dpos
    fileXfoilCoords=os.path.join(pathXfoilFiles,'coordsXfoil_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat')
    fileCoeffs = os.path.join(pathXfoilFiles,'coeffs_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat') #added dpos
    fileOut = os.path.join(pathXfoilFiles,'xfoilOut_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.txt')
    filePlot = os.path.join(pathXfoilFiles,'plot_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.ps')

    flagNoConvXfoil = False
    
    #-----------------------------------------------------

    # remove if files already exist
    listFile2Remove=[filePressure,fileXfoilCoords,fileCoeffs,fileOut,filePlot]
	
    try:
        for j in listFile2Remove:
            if os.path.exists(j)==True:
                os.remove(j)
    except OSError:
        pass

    #-----------------------------------------------------

    log = open(fileOut,'w+')

    ps = sp.Popen([pathXfoilProgram],
                  stdin=sp.PIPE,
                  stdout=log,
                  stderr=None)

    issueCmd('PLOP')
    issueCmd('G')
    issueCmd('')
    issueCmd('LOAD')         
    issueCmd(fileCoords)
    issueCmd('PANE') # regenerate paneling
    issueCmd('SAVE')
    issueCmd(fileXfoilCoords)
    issueCmd('OPER')
    issueCmd('ITER')
    issueCmd(str(iterNo)) 
    issueCmd('VISC')
    issueCmd(str(Re))
    issueCmd('MACH')
    issueCmd(str(Ma))
    issueCmd('PACC')
    issueCmd(fileCoeffs)
    issueCmd('\n')      
    issueCmd('OPER')
    issueCmd('ASEQ' + ' ' + str(amin)+' '+str(amax)+' '+str(astep))
    issueCmd('PACC')
    issueCmd('\n')
    issueCmd('QUIT')
    
    ps.wait()
        
    log.seek(0)
    if 'VISCAL:  Convergence failed' in log.read():
        flagNoConvXfoil = True
        # os.remove(filePressure)
        # os.remove(fileCoeffs)
    log.close()
	
    return flagNoConvXfoil	
	
	
#######################################################################
# EXTRACT AERODYNAMIC COEFFICIENTS FROM FILE

def extractCoeffs(pathXfoilFiles, jobNumber, jobIndex, dpos):

    fileCoeffs = os.path.join(pathXfoilFiles,'coeffs_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat')

    with open(fileCoeffs) as infile:
        for line in infile:
            pass
        last = (line.strip()).split()
    
    # [CL, CD]
    aeroCoeffs = [float(last[1]), float(last[2])]

    return aeroCoeffs

#######################################################################
# EXTRACT XFOIL COEFFS POLAR
	
def extractCoeffsPolar(pathXfoilFiles, jobNumber, jobIndex, dpos):

    fileCoeffs = os.path.join(pathXfoilFiles,'coeffs_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat') # added dpos
	
    with open(fileCoeffs) as infile:
        coeffsLine = []
        aeroCoeffs = []
        i = 0
        for line in infile:
            coeffsLine.append(line)
            coeffsLine = (line.strip()).split()	
            if i > 11:
				aeroCoeffs.append(0)
				aeroCoeffs[i-12] = [float(coeffsLine[0]), float(coeffsLine[1])]
            i = i + 1

    return aeroCoeffs		
	
#######################################################################
# EXTRACT XFOIL COEFFS POLAR
	
def extractCoeffsPolarPRE(pathXfoilFiles, jobNumber, jobIndex, dpos):

    fileCoeffs = os.path.join(pathXfoilFiles,'coeffs_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat') # added dpos
	
    with open(fileCoeffs) as infile:
        coeffsLine = []
        aeroCoeffs = []
        i = 0
        for line in infile:
            coeffsLine.append(line)
            coeffsLine = (line.strip()).split()	
            if i > 11:
				aeroCoeffs.append(0)
				#aeroCoeffs[i-12] = [float(coeffsLine[0]), float(coeffsLine[1])]
				aeroCoeffs[i-12] = [float(coeffsLine[0]), float(coeffsLine[1]), float(coeffsLine[2])]
            i = i + 1

    return aeroCoeffs	
	
#######################################################################
# Linear Weissinger

def weissinger(chord,depthWing,alpha,al0,ls,n):
	
    #span = 1.0*depthWing/1000
    span = 1.0*depthWing
    step = span/n
    
    alpha0half = al0[::-1]
    alpha0 = [0]*n
    for i in range(n/2):
		alpha0[i] = alpha0half[i]
		alpha0[i+n/2] = al0[i]
	
    alphatot = [0]*(n)
    for i in range(n):
		alphatot[i] = (alpha*pi/180) - alpha0[i]
		
    liftslopehalf = ls[::-1]
    liftslope = [0]*n
    for i in range(n/2):
		liftslope[i] = liftslopehalf[i]
		liftslope[i+n/2] = ls[i]
	
	
    if max(alphatot) == 0:
		cl = [0]*(n/2)
		cl_tot = [0]*n
		alphaeff = [0]*(n/2)
		alphaeff_tot = [0]*n
		w_ij_trailing = numpy.zeros((n,n))
		gamma = [0]*n
    else:
		# vortex coordinates: 1/4 chord
		x_vortex = [0.0]*(n+1)
		# y_vortex = [-(span/2)+i/step for i in range(n+1) ]
		y_vortex = [-(span/2)+i*step for i in range(n+1) ]
		z_vortex = [0.0]*(n+1)
		
		# control point coordinates: 3/4 chord
		#x_control = [chord/1000/2.0]*(n)
		x_control = [chord/2.0]*(n)
		y_control = [0.0]*(n)
		for i in range(n):
			y_control[i] = 0.5*(y_vortex[i]+y_vortex[i+1])
		z_control = [0.0]*(n)

		# normal on control surface
		n_x = [sin(alphatot[x]) for x in range(len(alphatot))]
		n_y = [0.0]*n
		n_z = [cos(alphatot[x]) for x in range(len(alphatot))]
		
		# freestream velocity normalized
		u_inf = 1.0
		v_inf = 0.0
		w_inf = 0.0
		
		# BC: Sum of the normal velocity component induced by the wing's bound vortices w_b, by the wake w_i
		#     and by the free-stream velocity Q_inf will be zero
		#     w_b + w_i + Q_inf*alpha = 0
		RHS = [0]*n   # right hand side: -Q_inf*alpha 
		
		# Induced Velocity by Vortex i at Control Point j
		w_ij = numpy.zeros((n,n))
		w_ij_trailing = numpy.zeros((n,n))

		for j in range(n):
			RHS[j] = -(n_x[j]*u_inf+n_y[j]*v_inf+n_z[j]*w_inf)
			for i in range(n):
				# Induced Velocity by Horseshoe vortice for Lift Calculation
				w_ij_HS = horseshoe(x_control[j],y_control[j],z_control[j],x_vortex[i],y_vortex[i],z_vortex[i],x_vortex[i+1],y_vortex[i+1],z_vortex[i+1])
				w_ij[j,i] = w_ij_HS[0]*n_x[j] + w_ij_HS[1]*n_y[j] + w_ij_HS[2]*n_z[j]
				# Induced Velocity by Horseshoe Trailing vortice for Drag Calculation
				w_ij_HST = horseshoetrailing(x_control[j],y_control[j],z_control[j],x_vortex[i],y_vortex[i],z_vortex[i],x_vortex[i+1],y_vortex[i+1],z_vortex[i+1])
				w_ij_trailing[j,i] = w_ij_HST[0]*n_x[j] + w_ij_HST[1]*n_y[j] + w_ij_HST[2]*n_z[j]

		# Circulation gamma
		gamma = numpy.linalg.solve(w_ij,RHS)
		
		# Lift Coefficient and Effective AoA
		cl = [0]*(n/2)
		alphaeff = [0]*(n/2)
		for i in range(n/2):
			#cl[i] = 2/(chord/1000)*gamma[n/2+i]
			cl[i] = 2/(chord)*gamma[n/2+i]
			alphaeff[i] = cl[i]/liftslope[n/2+i] + alpha0[n/2+i]
        
		cl_tot = [0]*n
		alphaeff_tot = [0]*n
		for i in range(n):
			#cl_tot[i] = 2/(chord/1000)*gamma[i]
			cl_tot[i] = 2/(chord)*gamma[i]
			alphaeff_tot[i] = cl_tot[i]/liftslope[i] + alpha0[i]

		
    # # f = open('C:\Temp\script_opt1_ML\scripts\ll.txt','a')
    # f = open(folder_path+'/AeroCoeff/ll.txt','a')
    # for i in range(n):    
		# f.write('cl_y = '+ str(cl_tot[i])+ '\n')
    # f.write('\n')
    # for i in range(n):    
		# f.write('alphaeff = '+ str(alphaeff_tot[i])+ '\n')
    # f.write('\n')
    # f.close()
		
    return alphaeff_tot, w_ij_trailing, gamma

	
#######################################################################
# Non Linear Weissinger

def nlweissinger(alphaeff, w_ij_trailing, gamma, alphaTable, clTable, chord, alpha, alpha0, n):
    
    # alphatot = alpha*pi/180 - alpha0
    
    D = 0.05 # Damping Factor
    k = 1
    k_TOL = 10**-5
    runI = 0
		
    gammaNL = [0]*n
    gammaNLD = [0]*n
    cl = [0]*n
    alphai = [0]*n
    kt = [0]*n
	
    while k > k_TOL and runI < 50000:
		runI = runI+1
		if runI == 50000:
			alphaeffOLD = alphaeff
			print('alphaeffOLD = '+str(alphaeffOLD))
		for i in range(n):
			if i < n/2: 
				cl[i] = numpy.interp(alphaeff[i]*180/pi,alphaTable[n/2-1-i],clTable[n/2-1-i])
			else:
				cl[i] = numpy.interp(alphaeff[i]*180/pi,alphaTable[i-n/2],clTable[i-n/2])
			# Kutta Joukowsi: 
			# dL = rho*V*Gamma(y)dy = cl*rho/2*V^2*dy*c 
			# -> Gamma = (c*V)/2*cl (V=1) 
			#gammaNL[i] = cl[i]*chord/1000/2
			gammaNL[i] = cl[i]*chord/2
			gammaNLD[i] = (1-D)*gamma[i] + D*gammaNL[i]
			kt[i] = abs(gamma[i] - gammaNL[i])
	
		maxkt = max(kt)
		# # f = open('C:\Temp\script_opt1_ML\scripts\gamma.txt','a')
		# f = open(folder_path+'/AeroCoeff/gamma.txt','a')	
		# f.write('delta gamma max = '+ str(maxkt)+ '\n')
		# f.close()
	
		k = max(kt)
		# k = k - 0.1
		alphaizw = -numpy.dot(w_ij_trailing,gammaNLD)
		
		for i in range(n):
			alphai[i] = atan(alphaizw[i])
			alphaeff[i] = alpha*pi/180 - alphai[i] # coupled with xfoil alpha_a instead of alpha_tot
			gamma[i] = gammaNLD[i]
		if runI > 50000-50:
			print(gamma)
		if runI == 50000:
			print('alphaeffNEW = '+str(alphaeff))
			f = open('D:/faselu/ftero_f_testW/testWa.txt','w')
			f.write(str(alphaeff)+'\n'+str(alphaeffOLD)+'\n'+str(gammaNL)+'\n'+str(gamma)+'\n'+str(kt)+'\n')
			f.close()
    
	# # Lift Coefficient and Effective AoA
    # f = open('C:\Temp\script_opt1_ML\scripts\gamma.txt','a')
    # f.write('nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn' + '\n')
    # f.close()
	
    cl_nl = [0]*(n/2)
    alphaeff_nl = [0]*(n/2)
    for i in range(n/2):
		cl_nl[i] = cl[n/2+i]
		alphaeff_nl[i] = alphaeff[n/2+i]
	
    # f = open('C:\Temp\script_opt1_ML\scripts\ll.txt','a')
    # for i in range(len(cl)):    
		# f.write('cl_y_nl = '+ str(cl[i])+ '\n')
    # f.write('\n')
    # for i in range(len(cl)):    
        # f.write('alphaeff_nl = '+ str(alphaeff[i])+ '\n')
    # f.write('\n')
    # f.close()
	
    # clWing = numpy.mean(cl)

    return alphaeff_nl


#######################################################################
# Non Linear Weissinger FOR OPTIMIZATION	
	
def nlweissingerOPT(alphaeff, w_ij_trailing, gamma, alphaTable, clTable, chord, alpha, alpha0, n):
    
    # alphatot = alpha*pi/180 - alpha0
    
    D = 0.05 # Damping Factor
    k = 1
    k_TOL = 10**-5
    runI = 0
    nMAX = 5000#20000
		
    gammaNL = [0]*n
    gammaNLD = [0]*n
    cl = [0]*n
    alphai = [0]*n
    kt = [0]*n
    alphaeffOLD = [0]*n
    clOLD = [0]*n
    alphaeffMEAN = [0]*n
    clMEAN = [0]*n
	
    while k > k_TOL and runI < nMAX:
		runI = runI+1
		if runI == nMAX:
			for i in range(n):
				alphaeffOLD[i] = alphaeff[i]
				clOLD[i] = cl[i]
			#print('alphaeffOLDup = '+str(alphaeffOLD))
		for i in range(n):
			if i < n/2: 
				cl[i] = numpy.interp(alphaeff[i]*180/pi,alphaTable[n/2-1-i],clTable[n/2-1-i])
			else:
				cl[i] = numpy.interp(alphaeff[i]*180/pi,alphaTable[i-n/2],clTable[i-n/2])
			# Kutta Joukowsi: 
			# dL = rho*V*Gamma(y)dy = cl*rho/2*V^2*dy*c 
			# -> Gamma = (c*V)/2*cl (V=1) 
			#gammaNL[i] = cl[i]*chord/1000/2
			gammaNL[i] = cl[i]*chord/2
			gammaNLD[i] = (1-D)*gamma[i] + D*gammaNL[i]
			kt[i] = abs(gamma[i] - gammaNL[i])
	
		maxkt = max(kt)
		# # f = open('C:\Temp\script_opt1_ML\scripts\gamma.txt','a')
		# f = open(folder_path+'/AeroCoeff/gamma.txt','a')	
		# f.write('delta gamma max = '+ str(maxkt)+ '\n')
		# f.close()
	
		k = max(kt)
		# k = k - 0.1
		alphaizw = -numpy.dot(w_ij_trailing,gammaNLD)
		
		for i in range(n):
			alphai[i] = atan(alphaizw[i])
			alphaeff[i] = alpha*pi/180 - alphai[i] # coupled with xfoil alpha_a instead of alpha_tot
			gamma[i] = gammaNLD[i]
		#if runI > nMAX-50:
		#	print(gamma[0])
		if runI == nMAX:
			for i in range(n):
				alphaeffMEAN[i] = (alphaeff[i]+alphaeffOLD[i])/2
				clMEAN[i] = (cl[i]+clOLD[i])/2
			#print('alphaeffNEW = '+str(alphaeff))
			#print('alphaeffOLD = '+str(alphaeffOLD))
			#f = open('D:/faselu/ftero_f_testW/testWa.txt','w')
			#f.write(str(alphaeff)+'\n'+str(alphaeffOLD)+'\n'+str(gammaNL)+'\n'+str(gamma)+'\n'+str(kt)+'\n')
			#f.close()
    
	# # Lift Coefficient and Effective AoA
    # f = open('C:\Temp\script_opt1_ML\scripts\gamma.txt','a')
    # f.write('nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn' + '\n')
    # f.close()
	
    cl_nl = [0]*(n/2)
    alphaeff_nl = [0]*(n/2)
    if runI < nMAX:
    	for i in range(n/2):
    		cl_nl[i] = cl[n/2+i]
    		alphaeff_nl[i] = alphaeff[n/2+i]
    else:
    	for i in range(n/2):
    		cl_nl[i] = clMEAN[n/2+i]
    		alphaeff_nl[i] = alphaeffMEAN[n/2+i]
			
    # f = open('C:\Temp\script_opt1_ML\scripts\ll.txt','a')
    # for i in range(len(cl)):    
		# f.write('cl_y_nl = '+ str(cl[i])+ '\n')
    # f.write('\n')
    # for i in range(len(cl)):    
        # f.write('alphaeff_nl = '+ str(alphaeff[i])+ '\n')
    # f.write('\n')
    # f.close()
	
    # clWing = numpy.mean(cl)

    return alphaeff_nl	
	
#######################################################################
# CALCULATE INDUCED VELOCITY OF A HORSESHOE VORTICE

def horseshoe(X,Y,Z,X1,Y1,Z1,X2,Y2,Z2):

    infinity = 10.0**8

    V1 = vortex(X,Y,Z,infinity,Y1,Z1,X1,Y1,Z1)
    V2 = vortex(X,Y,Z,X1,Y1,Z1,X2,Y2,Z2)
    V3 = vortex(X,Y,Z,X2,Y2,Z2,infinity,Y2,Z2)

    vel = [V1[0]+V2[0]+V3[0], V1[1]+V2[1]+V3[1],  V1[2]+V2[2]+V3[2]]

    return vel

#######################################################################
# CALCULATE INDUCED VELOCITY OF THE TRAILING HORSESHOE VORTICES

def horseshoetrailing(X,Y,Z,X1,Y1,Z1,X2,Y2,Z2):

    infinity = 10.0**8

    V1 = vortex(X,Y,Z,infinity,Y1,Z1,X1,Y1,Z1)
    V3 = vortex(X,Y,Z,X2,Y2,Z2,infinity,Y2,Z2)

    vel = [V1[0]+V3[0], V1[1]+V3[1], V1[2]+V3[2]]

    return vel
   
#######################################################################
# CALCULATE INDUCED VELOCITY OF A VORTEX

def vortex(X,Y,Z,X1,Y1,Z1,X2,Y2,Z2):
	
    TOL=10.0**-10;
    # CALCULATION OF R1 X R2 
    R1R2X = (Y-Y1)*(Z-Z2) - (Z-Z1)*(Y-Y2)
    R1R2Y = (X-X1)*(Z-Z2) - (Z-Z1)*(X-X2)
    R1R2Z = (X-X1)*(Y-Y2) - (Y-Y1)*(X-X2)

    # CALCULATION OF (Ri X R2 )**2 
    SQUARE = R1R2X*R1R2X + R1R2Y*R1R2Y + R1R2Z*R1R2Z

    # CALCULATION OF RO(R1/R(R1)-R2/R(R2))
    # R1 = sqrt((X-X1)*(X-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1))
    # R2 = sqrt((X-X2)*(X-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2))

    R1s = (X-X1)*(X-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1)
    R2s = (X-X2)*(X-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2)
    R1 = sqrt(R1s)
    R2 = sqrt(R2s)
	
    if (R1 < TOL or (R2 < TOL or SQUARE < TOL)):
        vel = [0, 0, 0]
    else:
        ROR1=(X2-X1)*(X-X1) + (Y2-Y1)*(Y-Y1) + (Z2-Z1)*(Z-Z1)
        ROR2=(X2-X1)*(X-X2) + (Y2-Y1)*(Y-Y2) + (Z2-Z1)*(Z-Z2) 
    
        COEFF = 1/(4*pi*SQUARE)*(ROR1/R1-ROR2/R2)
        vel = [R1R2X*COEFF, R1R2Y*COEFF, R1R2Z*COEFF]

    return vel


#######################################################################
# XFOIL PRESSURE COEFFICIENT DATA

def dataAeroCpRootold(pathXfoilFiles,jobNumber,jobIndex,chordLength, dpos):

    # read pressure coefficients from xfoil file
    filePressure = os.path.join(pathXfoilFiles,'pressure_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat')

    xfoilCp = []
      
	# read xfoil pressure from file
    with open(filePressure) as f:
        next(f)
        for line in f:
            xfoilCp.append([float(x) for x in line.split()])
    f.close()

    # scale xfoil x coordinates to real length
    xfoilCp = [[aa*chordLength,bb] for [aa, bb] in xfoilCp]
    
    # divide into top/bottom
    xfoilCp_top = xfoilCp[:(len(xfoilCp)/2)]
    xfoilCp_bottom = xfoilCp[(len(xfoilCp)/2):]

    # reverse direction of top surface vector (LE to TE)
    xfoilCp_top.sort(key=lambda x: x[0])  

    return xfoilCp_top, xfoilCp_bottom
 
 #######################################################################
# XFOIL PRESSURE COEFFICIENT DATA ROOT CLAMPING -> PRESSURE ON D-SPAR

def dataAeroCpRoot(pathXfoilFiles,jobNumber,jobIndex,chordLength, dpos):
# xDspar = [-coords_top_web1[0],-coords_bottom_web1[0]]

    # read pressure coefficients from xfoil file
    filePressure = os.path.join(pathXfoilFiles,'pressure_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat')

    # tolerance for finding the closest node to the beginning of D-spar
    # tolDspar = 8.0

    xfoilCp = []
    # xfoilCpNoDspar_top = []
    # xfoilCpNoDspar_bottom = []
      
# read xfoil pressure from file
    with open(filePressure) as f:
        next(f)
        for line in f:
            # xfoilCp.append([float(x) for x in line.split()])
            try:
				xfoilCp.append([float(x) for x in line.split()])
            except:
				xfoilCp.append([float(x) for x in line.split('-')])
				xfoilCp[len(xfoilCp)-1][1] = -xfoilCp[len(xfoilCp)-1][1]			
    f.close()

    # scale xfoil x coordinates to real length
    xfoilCp = [[aa*chordLength,bb] for [aa, bb] in xfoilCp]
    
    # divide into top/bottom
    xfoilCp_top = xfoilCp[:(len(xfoilCp)/2)]
    xfoilCp_bottom = xfoilCp[(len(xfoilCp)/2):]
    
    # select only pressure data outside D-spar
    # for value in xfoilCp_top:
        # if value[0] >= xDspar[0]-tolDspar:
            # xfoilCpNoDspar_top.append(value)
    # for value in xfoilCp_bottom:
        # if value[0] >= xDspar[1]-tolDspar:
            # xfoilCpNoDspar_bottom.append(value)

    # reverse direction of top surface vector (LE to TE)
    xfoilCp_top.sort(key=lambda x: x[0])  

    return xfoilCp_top, xfoilCp_bottom
 
 #######################################################################
# PRESSURE INTERPOLATION AND APPLICATION ROOT CLAMPING -> PRESSURE ON D-SPAR

def pressureApplyFalk(topCp,bottomCp,aerofoilCoords,depthProfile,stepNo,dynPressure,aeroModel,aeroAssembly,aeroInstance,lldis,chord,roin,pNorm,folder_path,t_step):

    dtolSurf = 0.0001

    currentStep = 'Step-' + str(stepNo)

    # slice through lists to get vectors: xfoil
    xfoilX_top = []
    xfoilCp_top = []
    xfoilX_bott = []
    xfoilCp_bott = []
	
	# append empty lists in first lldis indexes
    for i in range(lldis):
		xfoilX_top.append([])
		xfoilCp_top.append([])
		xfoilX_bott.append([])
		xfoilCp_bott.append([])

	# append x and cp from xfoil data
    for j in range(lldis):
		for item in topCp[j]:
			xfoilX_top[j].append(item[0])
			xfoilCp_top[j].append(item[1])
		for item in bottomCp[j]:
			xfoilX_bott[j].append(item[0])
			xfoilCp_bott[j].append(item[1])
 
    # slice through lists to get vectors: fe mesh
    # aerofoil data
    feX_top = []
    feX_bott = []
    feZ_top = []
    feZ_bott = []
	
    for i in range(lldis):
		feX_top.append([])
		feX_bott.append([])
		feZ_top.append([])
		feZ_bott.append([])
		# s_top = 0
		# s_bott = 0
		for j in range(len(aerofoilCoords[i])):
			if aerofoilCoords[i][j][1] > 0:
				feX_top[i].append(aerofoilCoords[i][j][0])
				feZ_top[i].append(aerofoilCoords[i][j][1])
				# s_top += 1
			else:
				feX_bott[i].append(aerofoilCoords[i][j][0])
				feZ_bott[i].append(aerofoilCoords[i][j][1])
				# s_bott += 1

			
	# interpolate cp
    feCp_top = []
    feCp_bott = []
	
    fePress_top = []
    fePress_bott =[]
	
    for j in range(lldis):
		feCp_top.append([])
		feCp_bott.append([])
		feCp_top[j] = numpy.interp(feX_top[j], xfoilX_top[j], xfoilCp_top[j], left=None, right=None)
		feCp_bott[j] = numpy.interp(feX_bott[j], xfoilX_bott[j], xfoilCp_bott[j], left=None, right=None)

		# change Cp to pressure
		fePress_top.append([])
		fePress_bott.append([])
		fePress_top[j] = [x * dynPressure for x in feCp_top[j]]
		fePress_bott[j] = [x * dynPressure for x in feCp_bott[j]]

    # define one list with top and bottom data
    fePress = []
    aerofoilDataX = []
    aerofoilDataZ = []
	
    for j in range(lldis):
		fePress.append([])
		aerofoilDataX.append([])
		aerofoilDataZ.append([])
		for item in fePress_top[j]:
			fePress[j].append(item)
		for item in fePress_bott[j]:
			fePress[j].append(item)
		for item in feX_top[j]:
			aerofoilDataX[j].append(item)
		for item in feX_bott[j]:
			aerofoilDataX[j].append(item)	
		for item in feZ_top[j]:
			aerofoilDataZ[j].append(item)
		for item in feZ_bott[j]:
			aerofoilDataZ[j].append(item)	
    
	#-----------------------------------------------------
    # apply interpolated pressure to pressure surfaces
    n_skinpts = len(fePress[0])   

    read_out = [[0 for x in range(4)] for x in range(n_skinpts*lldis)]

	# pressure in cae is applied on undeformed geometry, therefore the geometry data only needs to be calculated in step 1
    if stepNo == 1:
		for y in range(lldis):
			for x in range(1,n_skinpts+1):
				read_out[x-1+y*n_skinpts][0] = float(aerofoilDataX[y][x-1])
				read_out[x-1+y*n_skinpts][2] = float(depthProfile/lldis*(y + 0.5))
				read_out[x-1+y*n_skinpts][1] = float(aerofoilDataZ[y][x-1])
				read_out[x-1+y*n_skinpts][3] = float(fePress[y][x-1])
				
				f = open(folder_path+'/AeroCoeff/p.txt','a')  
				f.write('x = '+ str(read_out[x-1+y*n_skinpts][0])+ '\t')
				f.write('y = '+ str(read_out[x-1+y*n_skinpts][1])+ '\t')
				f.write('z = '+ str(read_out[x-1+y*n_skinpts][2])+ '\t')
				f.write('p = '+ str(read_out[x-1+y*n_skinpts][3])+ '\n')
				f.close()	
    else:
		for y in range(lldis):
			for x in range(1,n_skinpts+1):
				read_out[x-1+y*n_skinpts][0] = roin[x-1+y*n_skinpts][0]
				read_out[x-1+y*n_skinpts][1] = roin[x-1+y*n_skinpts][1]
				read_out[x-1+y*n_skinpts][2] = roin[x-1+y*n_skinpts][2]
				read_out[x-1+y*n_skinpts][3] = float(fePress[y][x-1])
				
				f = open(folder_path+'/AeroCoeff/p.txt','a')  
				f.write('x = '+ str(read_out[x-1+y*n_skinpts][0])+ '\t')
				f.write('y = '+ str(read_out[x-1+y*n_skinpts][1])+ '\t')
				f.write('z = '+ str(read_out[x-1+y*n_skinpts][2])+ '\t')
				f.write('p = '+ str(read_out[x-1+y*n_skinpts][3])+ '\n')
				f.close()	
	

    read_out= tuple(read_out)

    # aeroModel.MappedField(name='AD_field3', xyzPointData=read_out )
    aeroModel.MappedField(name='AD_field-'+str(stepNo), xyzPointData=read_out )
    
    factor_pressure = 1.0*(stepNo*t_step)
    if factor_pressure > 1:
		factor_pressure = 1
	
    if stepNo == 1:
		pNorm = [0]*(len(aerofoilCoords[0])-1)
		for i in range(len(aerofoilCoords[0])-1):
			xmina = aerofoilCoords[0][i][0]
			xmaxa = aerofoilCoords[0][i+1][0]
			zmina = aerofoilCoords[0][i][1]
			zmaxa = aerofoilCoords[0][i+1][1]
			ymin = -dtolSurf
			ymax = depthProfile+dtolSurf
			xmin = min(xmina, xmaxa)
			xmax = max(xmina, xmaxa)
			zmin = min(zmina, zmaxa)
			zmax = max(zmina, zmaxa)
			aeroAssembly.Surface( name='pSurf-'+str(i), 
				side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )
			
			# pSset for ODB PDLOAD calculation
			aeroAssembly.Surface( name='pSset-'+str(i), 
				side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )	

			if ((zmin >= 0.0) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((zmin < 0.0) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [pi-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((xmaxa-xmina) == 0):
				pNorm[i] = [-pi/2,sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
					
			if i > 0:
				aeroAssembly.SurfaceByBoolean(name='pSurf-0', surfaces=(aeroAssembly.surfaces['pSurf-0'], 
					aeroAssembly.surfaces['pSurf-'+str(i)], ))
				del aeroAssembly.surfaces['pSurf-'+str(i)]
			# ### OPTIONAL
			# # remove double selection with flange
			# aeroAssembly.Surface(name='pSurf-0', 
				# side2Elements=aeroInstance.elements.sequenceFromLabels(
				# labels=range(833,841), ))   
			# f = open('C:\Temp\script_opt1_ML2\scripts\pODB.txt','a')
			# f.write(str(pNorm[i])+ '\n')
			# f.close()
    else:	
		aeroModel.loads['press'+str(stepNo-1)].deactivate(currentStep)
		
    # aeroModel.Pressure(name='press'+str(stepNo), 
		# createStepName=currentStep, distributionType=FIELD, field='AD_field3', magnitude=factor_pressure,
		# region=aeroAssembly.surfaces['pSurf-0'])
    aeroModel.Pressure(name='press'+str(stepNo), 
		createStepName=currentStep, distributionType=FIELD, field='AD_field-'+str(stepNo), magnitude=factor_pressure,
		region=aeroAssembly.instances['Wing-final-1'].surfaces['Surf-Skin'])#['pSurf-0']) #Alejandro Valverde

    return read_out, pNorm

	
 #######################################################################
# PRESSURE INTERPOLATION AND APPLICATION ROOT CLAMPING -> PRESSURE ON D-SPAR

def pressureApplyFalkcamber(topCp,bottomCp,aerofoilCoords,depthProfile,stepNo,dynPressure,aeroModel,aeroAssembly,aeroInstance,lldis,chord,roin,pNorm,folder_path,t_step,xCL,yCL):

    dtolSurf = 0.0001

    currentStep = 'Step-' + str(stepNo)

    # slice through lists to get vectors: xfoil
    xfoilX_top = []
    xfoilCp_top = []
    xfoilX_bott = []
    xfoilCp_bott = []
	
	# append empty lists in first lldis indexes
    for i in range(lldis):
		xfoilX_top.append([])
		xfoilCp_top.append([])
		xfoilX_bott.append([])
		xfoilCp_bott.append([])

	# append x and cp from xfoil data
    for j in range(lldis):
		for item in topCp[j]:
			xfoilX_top[j].append(item[0])
			xfoilCp_top[j].append(item[1])
		for item in bottomCp[j]:
			xfoilX_bott[j].append(item[0])
			xfoilCp_bott[j].append(item[1])
 
    # slice through lists to get vectors: fe mesh
    # aerofoil data
    feX_top = []
    feX_bott = []
    feZ_top = []
    feZ_bott = []
	
    for i in range(lldis):
		feX_top.append([])
		feX_bott.append([])
		feZ_top.append([])
		feZ_bott.append([])
		# s_top = 0
		# s_bott = 0
		for j in range(len(aerofoilCoords[i])):
			#if aerofoilCoords[i][j][1] > 0:
			#	feX_top[i].append(aerofoilCoords[i][j][0])
			#	feZ_top[i].append(aerofoilCoords[i][j][1])
			#	# s_top += 1
			#else:
			#	feX_bott[i].append(aerofoilCoords[i][j][0])
			#	feZ_bott[i].append(aerofoilCoords[i][j][1])
			#	# s_bott += 1
			if numpy.interp(aerofoilCoords[i][j][0], xCL, yCL) <= aerofoilCoords[i][j][1]:
				feX_top[i].append(aerofoilCoords[i][j][0])
				feZ_top[i].append(aerofoilCoords[i][j][1])
				# s_top += 1
			else:
				feX_bott[i].append(aerofoilCoords[i][j][0])
				feZ_bott[i].append(aerofoilCoords[i][j][1])
				# s_bott += 1	
			
	# interpolate cp
    feCp_top = []
    feCp_bott = []
	
    fePress_top = []
    fePress_bott =[]
	
    for j in range(lldis):
		feCp_top.append([])
		feCp_bott.append([])
		feCp_top[j] = numpy.interp(feX_top[j], xfoilX_top[j], xfoilCp_top[j], left=None, right=None)
		feCp_bott[j] = numpy.interp(feX_bott[j], xfoilX_bott[j], xfoilCp_bott[j], left=None, right=None)

		# change Cp to pressure
		fePress_top.append([])
		fePress_bott.append([])
		fePress_top[j] = [x * dynPressure for x in feCp_top[j]]
		fePress_bott[j] = [x * dynPressure for x in feCp_bott[j]]

    # define one list with top and bottom data
    fePress = []
    aerofoilDataX = []
    aerofoilDataZ = []
	
    for j in range(lldis):
		fePress.append([])
		aerofoilDataX.append([])
		aerofoilDataZ.append([])
		for item in fePress_top[j]:
			fePress[j].append(item)
		for item in fePress_bott[j]:
			fePress[j].append(item)
		for item in feX_top[j]:
			aerofoilDataX[j].append(item)
		for item in feX_bott[j]:
			aerofoilDataX[j].append(item)	
		for item in feZ_top[j]:
			aerofoilDataZ[j].append(item)
		for item in feZ_bott[j]:
			aerofoilDataZ[j].append(item)	
    
	#-----------------------------------------------------
    # apply interpolated pressure to pressure surfaces
    n_skinpts = len(fePress[0])   

    read_out = [[0 for x in range(4)] for x in range(n_skinpts*lldis)]

	# pressure in cae is applied on undeformed geometry, therefore the geometry data only needs to be calculated in step 1
    if stepNo == 1:
		for y in range(lldis):
			for x in range(1,n_skinpts+1):
				read_out[x-1+y*n_skinpts][0] = float(aerofoilDataX[y][x-1])
				read_out[x-1+y*n_skinpts][2] = float(depthProfile/lldis*(y + 0.5))
				read_out[x-1+y*n_skinpts][1] = float(aerofoilDataZ[y][x-1])
				read_out[x-1+y*n_skinpts][3] = float(fePress[y][x-1])
				
				f = open(folder_path+'/AeroCoeff/p.txt','a')  
				f.write('x = '+ str(read_out[x-1+y*n_skinpts][0])+ '\t')
				f.write('y = '+ str(read_out[x-1+y*n_skinpts][1])+ '\t')
				f.write('z = '+ str(read_out[x-1+y*n_skinpts][2])+ '\t')
				f.write('p = '+ str(read_out[x-1+y*n_skinpts][3])+ '\n')
				f.close()	
    else:
		for y in range(lldis):
			for x in range(1,n_skinpts+1):
				read_out[x-1+y*n_skinpts][0] = roin[x-1+y*n_skinpts][0]
				read_out[x-1+y*n_skinpts][1] = roin[x-1+y*n_skinpts][1]
				read_out[x-1+y*n_skinpts][2] = roin[x-1+y*n_skinpts][2]
				read_out[x-1+y*n_skinpts][3] = float(fePress[y][x-1])
				
				f = open(folder_path+'/AeroCoeff/p.txt','a')  
				f.write('x = '+ str(read_out[x-1+y*n_skinpts][0])+ '\t')
				f.write('y = '+ str(read_out[x-1+y*n_skinpts][1])+ '\t')
				f.write('z = '+ str(read_out[x-1+y*n_skinpts][2])+ '\t')
				f.write('p = '+ str(read_out[x-1+y*n_skinpts][3])+ '\n')
				f.close()	
	

    read_out= tuple(read_out)

    # aeroModel.MappedField(name='AD_field3', xyzPointData=read_out )
    aeroModel.MappedField(name='AD_field-'+str(stepNo), xyzPointData=read_out )
    
#    factor_pressure = 1.0*(stepNo*t_step)
    factor_pressure = t_step
    factor_pressure = -factor_pressure #Changed, Alejandro Valverde
    if factor_pressure > 1:
		factor_pressure = 1
	
    if stepNo == 1:
		pNorm = [0]*(len(aerofoilCoords[0])-1)
		for i in range(len(aerofoilCoords[0])-1):
			xmina = aerofoilCoords[0][i][0]
			xmaxa = aerofoilCoords[0][i+1][0]
			zmina = aerofoilCoords[0][i][1]
			zmaxa = aerofoilCoords[0][i+1][1]
			ymin = -dtolSurf
			ymax = depthProfile+dtolSurf
			xmin = min(xmina, xmaxa)
			xmax = max(xmina, xmaxa)
			zmin = min(zmina, zmaxa)
			zmax = max(zmina, zmaxa)
			aeroAssembly.Surface( name='pSurf-'+str(i), 
				side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )
			
			# pSset for ODB PDLOAD calculation
			aeroAssembly.Surface( name='pSset-'+str(i), 
				side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )	

			#if ((zmin >= 0.0) and ((xmaxa-xmina) != 0)):
			#	pNorm[i] = [-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			#if ((zmin < 0.0) and ((xmaxa-xmina) != 0)):
			#	pNorm[i] = [pi-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((numpy.interp((xmaxa+xmina)/2, xCL, yCL) <= (zmaxa+zmina)/2) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((numpy.interp((xmaxa+xmina)/2, xCL, yCL) > (zmaxa+zmina)/2) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [pi-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((xmaxa-xmina) == 0):
				pNorm[i] = [-pi/2,sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
					
			if i > 0:
				aeroAssembly.SurfaceByBoolean(name='pSurf-0', surfaces=(aeroAssembly.surfaces['pSurf-0'], 
					aeroAssembly.surfaces['pSurf-'+str(i)], ))
				del aeroAssembly.surfaces['pSurf-'+str(i)]
			# ### OPTIONAL
			# # remove double selection with flange
			# aeroAssembly.Surface(name='pSurf-0', 
				# side2Elements=aeroInstance.elements.sequenceFromLabels(
				# labels=range(833,841), ))   
			# f = open('C:\Temp\script_opt1_ML2\scripts\pODB.txt','a')
			# f.write(str(pNorm[i])+ '\n')
			# f.close()
    else:	
		aeroModel.loads['press'+str(stepNo-1)].deactivate(currentStep)
		
    # aeroModel.Pressure(name='press'+str(stepNo), 
		# createStepName=currentStep, distributionType=FIELD, field='AD_field3', magnitude=factor_pressure,
		# region=aeroAssembly.surfaces['pSurf-0'])
    aeroModel.Pressure(name='press'+str(stepNo), 
		createStepName=currentStep, distributionType=FIELD, field='AD_field-'+str(stepNo), magnitude=factor_pressure,
		region=aeroAssembly.instances['Wing-final-1'].surfaces['Surf-Skin'])#['pSurf-0']) #Alejandro Valverde

    return read_out, pNorm
	
	
#######################################################################
# PRESSURE INTERPOLATION AND APPLICATION ROOT CLAMPING -> PRESSURE ON D-SPAR

def pressureApplyEV(topCp,bottomCp,aerofoilCoords,depthProfile,stepNo,dynPressure,aeroModel,aeroAssembly,aeroInstance,lldis,chord,folder_path,t_step):

    dtolSurf = 0.0001

    currentStep = 'Step-EV'

    # slice through lists to get vectors: xfoil
    xfoilX_top = []
    xfoilCp_top = []
    xfoilX_bott = []
    xfoilCp_bott = []
	
	# append empty lists in first lldis indexes
    for i in range(lldis):
		xfoilX_top.append([])
		xfoilCp_top.append([])
		xfoilX_bott.append([])
		xfoilCp_bott.append([])

	# append x and cp from xfoil data
    for j in range(lldis):
		for item in topCp[j]:
			xfoilX_top[j].append(item[0])
			xfoilCp_top[j].append(item[1])
		for item in bottomCp[j]:
			xfoilX_bott[j].append(item[0])
			xfoilCp_bott[j].append(item[1])
 
    # slice through lists to get vectors: fe mesh
    # aerofoil data
    feX_top = []
    feX_bott = []
    feZ_top = []
    feZ_bott = []
	
    for i in range(lldis):
		feX_top.append([])
		feX_bott.append([])
		feZ_top.append([])
		feZ_bott.append([])
		# s_top = 0
		# s_bott = 0
		for j in range(len(aerofoilCoords[i])):
			if aerofoilCoords[i][j][1] > 0:
				feX_top[i].append(aerofoilCoords[i][j][0])
				feZ_top[i].append(aerofoilCoords[i][j][1])
				# s_top += 1
			else:
				feX_bott[i].append(aerofoilCoords[i][j][0])
				feZ_bott[i].append(aerofoilCoords[i][j][1])
				# s_bott += 1

			
	# interpolate cp
    feCp_top = []
    feCp_bott = []
	
    fePress_top = []
    fePress_bott =[]
	
    for j in range(lldis):
		feCp_top.append([])
		feCp_bott.append([])
		feCp_top[j] = numpy.interp(feX_top[j], xfoilX_top[j], xfoilCp_top[j], left=None, right=None)
		feCp_bott[j] = numpy.interp(feX_bott[j], xfoilX_bott[j], xfoilCp_bott[j], left=None, right=None)

		# change Cp to pressure
		fePress_top.append([])
		fePress_bott.append([])
		fePress_top[j] = [x * dynPressure for x in feCp_top[j]]
		fePress_bott[j] = [x * dynPressure for x in feCp_bott[j]]

    # define one list with top and bottom data
    fePress = []
    aerofoilDataX = []
    aerofoilDataZ = []
	
    for j in range(lldis):
		fePress.append([])
		aerofoilDataX.append([])
		aerofoilDataZ.append([])
		for item in fePress_top[j]:
			fePress[j].append(item)
		for item in fePress_bott[j]:
			fePress[j].append(item)
		for item in feX_top[j]:
			aerofoilDataX[j].append(item)
		for item in feX_bott[j]:
			aerofoilDataX[j].append(item)	
		for item in feZ_top[j]:
			aerofoilDataZ[j].append(item)
		for item in feZ_bott[j]:
			aerofoilDataZ[j].append(item)	
    
	#-----------------------------------------------------
    # apply interpolated pressure to pressure surfaces
    n_skinpts = len(fePress[0])   

    read_out = [[0 for x in range(4)] for x in range(n_skinpts*lldis)]

	# pressure in cae is applied on undeformed geometry, therefore the geometry data only needs to be calculated in step 1
    if stepNo == 1:
		for y in range(lldis):
			for x in range(1,n_skinpts+1):
				read_out[x-1+y*n_skinpts][0] = float(aerofoilDataX[y][x-1])
				read_out[x-1+y*n_skinpts][2] = float(depthProfile/lldis*(y + 0.5))
				read_out[x-1+y*n_skinpts][1] = float(aerofoilDataZ[y][x-1])
				read_out[x-1+y*n_skinpts][3] = float(fePress[y][x-1])
				
				# f = open(folder_path+'/AeroCoeff/p.txt','a')  
				# f.write('x = '+ str(read_out[x-1+y*n_skinpts][0])+ '\t')
				# f.write('y = '+ str(read_out[x-1+y*n_skinpts][1])+ '\t')
				# f.write('z = '+ str(read_out[x-1+y*n_skinpts][2])+ '\t')
				# f.write('p = '+ str(read_out[x-1+y*n_skinpts][3])+ '\n')
				# f.close()	
	

    read_out= tuple(read_out)

    # aeroModel.MappedField(name='AD_field3', xyzPointData=read_out )
    aeroModel.MappedField(name='AD_fieldEV', xyzPointData=read_out )
    
    factor_pressure = 1.0
	
    if stepNo == 1:
		pNorm = [0]*(len(aerofoilCoords[0])-1)
		for i in range(len(aerofoilCoords[0])-1):
			xmina = aerofoilCoords[0][i][0]
			xmaxa = aerofoilCoords[0][i+1][0]
			zmina = aerofoilCoords[0][i][1]
			zmaxa = aerofoilCoords[0][i+1][1]
			ymin = -dtolSurf
			ymax = depthProfile+dtolSurf
			xmin = min(xmina, xmaxa)
			xmax = max(xmina, xmaxa)
			zmin = min(zmina, zmaxa)
			zmax = max(zmina, zmaxa)
			aeroAssembly.Surface( name='pSurf-'+str(i), 
				side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )
			# aeroAssembly.Surface( name='pSset-'+str(i), 
				# side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )	

			if ((zmin >= 0.0) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((zmin < 0.0) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [pi-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((xmaxa-xmina) == 0):
				pNorm[i] = [-pi/2,sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
					
			if i > 0:
				aeroAssembly.SurfaceByBoolean(name='pSurf-0', surfaces=(aeroAssembly.surfaces['pSurf-0'], 
					aeroAssembly.surfaces['pSurf-'+str(i)], ))
				del aeroAssembly.surfaces['pSurf-'+str(i)]
			# ### OPTIONAL
			# # remove double selection with flange
			# aeroAssembly.Surface(name='pSurf-0', 
				# side2Elements=aeroInstance.elements.sequenceFromLabels(
				# labels=range(833,841), ))   
			# f = open('C:\Temp\script_opt1_ML2\scripts\pODB.txt','a')
			# f.write(str(pNorm[i])+ '\n')
			# f.close()
		
    aeroModel.Pressure(name='press-EV', 
		createStepName=currentStep, distributionType=FIELD, field='AD_fieldEV', magnitude=factor_pressure,
		region=aeroAssembly.instances['Wing-final-1'].surfaces['Surf-Skin'])#['pSurf-0']) #Alejandro Valverde
 

#######################################################################
# PRESSURE INTERPOLATION AND APPLICATION ROOT CLAMPING -> PRESSURE ON D-SPAR  numpy.interp(value[0], xCL, yCL) <= value[1]:

def pressureApplyEVcamber(topCp,bottomCp,aerofoilCoords,depthProfile,stepNo,dynPressure,aeroModel,aeroAssembly,aeroInstance,lldis,chord,folder_path,t_step, xCL, yCL):

    dtolSurf = 0.0001

    currentStep = 'Step-EV'

    # slice through lists to get vectors: xfoil
    xfoilX_top = []
    xfoilCp_top = []
    xfoilX_bott = []
    xfoilCp_bott = []
	
	# append empty lists in first lldis indexes
    for i in range(lldis):
		xfoilX_top.append([])
		xfoilCp_top.append([])
		xfoilX_bott.append([])
		xfoilCp_bott.append([])

	# append x and cp from xfoil data
    for j in range(lldis):
		for item in topCp[j]:
			xfoilX_top[j].append(item[0])
			xfoilCp_top[j].append(item[1])
		for item in bottomCp[j]:
			xfoilX_bott[j].append(item[0])
			xfoilCp_bott[j].append(item[1])
 
    # slice through lists to get vectors: fe mesh
    # aerofoil data
    feX_top = []
    feX_bott = []
    feZ_top = []
    feZ_bott = []
	
    for i in range(lldis):
		feX_top.append([])
		feX_bott.append([])
		feZ_top.append([])
		feZ_bott.append([])
		# s_top = 0
		# s_bott = 0
		for j in range(len(aerofoilCoords[i])):
			#if aerofoilCoords[i][j][1] > 0:
			#	feX_top[i].append(aerofoilCoords[i][j][0])
			#	feZ_top[i].append(aerofoilCoords[i][j][1])
			#	# s_top += 1
			#else:
			#	feX_bott[i].append(aerofoilCoords[i][j][0])
			#	feZ_bott[i].append(aerofoilCoords[i][j][1])
			#	# s_bott += 1
			if numpy.interp(aerofoilCoords[i][j][0], xCL, yCL) <= aerofoilCoords[i][j][1]:
				feX_top[i].append(aerofoilCoords[i][j][0])
				feZ_top[i].append(aerofoilCoords[i][j][1])
				# s_top += 1
			else:
				feX_bott[i].append(aerofoilCoords[i][j][0])
				feZ_bott[i].append(aerofoilCoords[i][j][1])
				# s_bott += 1				

			
	# interpolate cp
    feCp_top = []
    feCp_bott = []
	
    fePress_top = []
    fePress_bott =[]
	
    for j in range(lldis):
		feCp_top.append([])
		feCp_bott.append([])
		feCp_top[j] = numpy.interp(feX_top[j], xfoilX_top[j], xfoilCp_top[j], left=None, right=None)
		feCp_bott[j] = numpy.interp(feX_bott[j], xfoilX_bott[j], xfoilCp_bott[j], left=None, right=None)

		# change Cp to pressure
		fePress_top.append([])
		fePress_bott.append([])
		fePress_top[j] = [x * dynPressure for x in feCp_top[j]]
		fePress_bott[j] = [x * dynPressure for x in feCp_bott[j]]

    # define one list with top and bottom data
    fePress = []
    aerofoilDataX = []
    aerofoilDataZ = []
	
    for j in range(lldis):
		fePress.append([])
		aerofoilDataX.append([])
		aerofoilDataZ.append([])
		for item in fePress_top[j]:
			fePress[j].append(item)
		for item in fePress_bott[j]:
			fePress[j].append(item)
		for item in feX_top[j]:
			aerofoilDataX[j].append(item)
		for item in feX_bott[j]:
			aerofoilDataX[j].append(item)	
		for item in feZ_top[j]:
			aerofoilDataZ[j].append(item)
		for item in feZ_bott[j]:
			aerofoilDataZ[j].append(item)	
    
	#-----------------------------------------------------
    # apply interpolated pressure to pressure surfaces
    n_skinpts = len(fePress[0])   

    read_out = [[0 for x in range(4)] for x in range(n_skinpts*lldis)]

	# pressure in cae is applied on undeformed geometry, therefore the geometry data only needs to be calculated in step 1
    if stepNo == 1:
		for y in range(lldis):
			for x in range(1,n_skinpts+1):
				read_out[x-1+y*n_skinpts][0] = float(aerofoilDataX[y][x-1])
				read_out[x-1+y*n_skinpts][2] = float(depthProfile/lldis*(y + 0.5))
				read_out[x-1+y*n_skinpts][1] = float(aerofoilDataZ[y][x-1])
				read_out[x-1+y*n_skinpts][3] = float(fePress[y][x-1])
				
				# f = open(folder_path+'/AeroCoeff/p.txt','a')  
				# f.write('x = '+ str(read_out[x-1+y*n_skinpts][0])+ '\t')
				# f.write('y = '+ str(read_out[x-1+y*n_skinpts][1])+ '\t')
				# f.write('z = '+ str(read_out[x-1+y*n_skinpts][2])+ '\t')
				# f.write('p = '+ str(read_out[x-1+y*n_skinpts][3])+ '\n')
				# f.close()	
	

    read_out= tuple(read_out)

    # aeroModel.MappedField(name='AD_field3', xyzPointData=read_out )
    aeroModel.MappedField(name='AD_fieldEV', xyzPointData=read_out )
    
    factor_pressure = 1.0
	
    if stepNo == 1:
		pNorm = [0]*(len(aerofoilCoords[0])-1)
		for i in range(len(aerofoilCoords[0])-1):
			xmina = aerofoilCoords[0][i][0]
			xmaxa = aerofoilCoords[0][i+1][0]
			zmina = aerofoilCoords[0][i][1]
			zmaxa = aerofoilCoords[0][i+1][1]
			ymin = -dtolSurf
			ymax = depthProfile+dtolSurf
			xmin = min(xmina, xmaxa)
			xmax = max(xmina, xmaxa)
			zmin = min(zmina, zmaxa)
			zmax = max(zmina, zmaxa)
			aeroAssembly.Surface( name='pSurf-'+str(i), 
				side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )
			# aeroAssembly.Surface( name='pSset-'+str(i), 
				# side2Elements=aeroInstance.elements.getByBoundingBox(xmin-dtolSurf,zmin-dtolSurf,ymin,xmax+dtolSurf,zmax+dtolSurf,ymax) )	

			if ((zmin >= 0.0) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((zmin < 0.0) and ((xmaxa-xmina) != 0)):
				pNorm[i] = [pi-atan((zmaxa-zmina)/(xmaxa-xmina)),sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
			if ((xmaxa-xmina) == 0):
				pNorm[i] = [-pi/2,sqrt((zmax-zmin)**2+(xmax-xmin)**2)]
					
			if i > 0:
				aeroAssembly.SurfaceByBoolean(name='pSurf-0', surfaces=(aeroAssembly.surfaces['pSurf-0'], 
					aeroAssembly.surfaces['pSurf-'+str(i)], ))
				del aeroAssembly.surfaces['pSurf-'+str(i)]
			# ### OPTIONAL
			# # remove double selection with flange
			# aeroAssembly.Surface(name='pSurf-0', 
				# side2Elements=aeroInstance.elements.sequenceFromLabels(
				# labels=range(833,841), ))   
			# f = open('C:\Temp\script_opt1_ML2\scripts\pODB.txt','a')
			# f.write(str(pNorm[i])+ '\n')
			# f.close()
		
    aeroModel.Pressure(name='press-EV', 
		createStepName=currentStep, distributionType=FIELD, field='AD_fieldEV', magnitude=factor_pressure,
		region=aeroAssembly.instances['Wing-final-1'].surfaces['Surf-Skin'])#['pSurf-0']) #Alejandro Valverde

#######################################################################
# NEW STEP GENERATION AND RESTART DEFINITION

def aeroLoopStepsNL(jobName,jobNo,stepNo,currentModel,dampingCoeff,initialIncrement,varsFieldOut,varsHistoryOut,stabMag):
    
    currentStep = 'Step-' + str(stepNo)

    if stepNo == 1:
        previousStep = 'Initial'
    else:
        previousStep = 'Step-' + str(stepNo-1)
        previousJobName =  jobName + '_' + str(jobNo) +'-' + str(stepNo-1)
        #-----------------------------------------------------
        # define where to restart analysis for higher steps
        currentModel.setValues( restartJob = previousJobName, restartStep = previousStep )

    #-----------------------------------------------------
    # create next step and restart
    stepLoop = currentModel.StaticStep(
           description='weakly coupled static aeroelastic analysis with xfoil', 
    ##           adaptiveDampingRatio=dampingCoeff,
           # continueDampingFactors=False,
           # stabilizationMagnitude=stabMag, 
           # stabilizationMethod=DAMPING_FACTOR,#DISSIPATED_ENERGY_FRACTION
           # initialInc=0.05,#initialIncrement,
           # maxInc=0.05,
           # minInc=10e-10,
           # maxNumInc=3000, 
           name=currentStep, 
           nlgeom=OFF, 
           previous=previousStep )
            
    #-----------------------------------------------------
    # restart
    stepLoop.Restart(frequency=1, numberIntervals=0, 
        overlay=ON, timeMarks=OFF)

    #-----------------------------------------------------
    # limit output requests
    if stepNo == 1:
        try:
            del currentModel.fieldOutputRequests['F-Output-1']
        except:
            pass
        currentModel.FieldOutputRequest(createStepName=currentStep, 
            name='FieldOutput', variables=varsFieldOut)
        mdb.models['Model-SparAngle-1'].fieldOutputRequests['FieldOutput'].setValues(sectionPoints=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
        # history output
        try:
            del currentModel.historyOutputRequests['H-Output-1']
        except:
            pass
        currentModel.HistoryOutputRequest(createStepName=currentStep, 
            name='HistoryOutput', variables=varsHistoryOut)			
			
			
#######################################################################
# DEFORMED AEROFOIL COORDINATES: DATA OUTPUT FROM *.ODB TO TEXT FILE

def coordsDeformedAerofoil(pathXfoilData,jobName,jobNumber,jobIndex,chordLength,dpos):

    # open output database
    odbName = jobName + '.odb'
    odb = visualization.openOdb(odbName)

    stepName = 'Step-%s' %jobIndex
	
    setName = 'SECTION%s' %dpos

    #-----------------------------------------------------
    # access relevant node sets
    # aerofoilY0 = odb.rootAssembly.instances['AEROFOILINSTANCE-1'].nodeSets['SETAEROFOILOUTEDGEY01']
    aerofoilY0 = odb.rootAssembly.instances['WING-FINAL-1'].nodeSets[setName]
	
    aerofoilCoordsX = []
    aerofoilCoordsZ = []

    # nodesAerofoil=aerofoilY0.nodes
    for k in aerofoilY0.nodes:
        aerofoilCoordsX.append(k.coordinates[0])
        aerofoilCoordsZ.append(k.coordinates[1])
          
    displAerofoilX = []
    displAerofoilZ = []

    dispField = odb.steps[stepName].frames[-1].fieldOutputs['U']
    dispSubField = dispField.getSubset(region=aerofoilY0)
    displValues = dispSubField.values
    for node in displValues:
        displAerofoilX.append(node.data[0])
        displAerofoilZ.append(node.data[1])
        
    odb.close()
    #-----------------------------------------------------
        
    aerofoilData = zip(aerofoilCoordsX,displAerofoilX,aerofoilCoordsZ,displAerofoilZ)

    aerofoilData_top = []
    aerofoilData_bottom = []

    # for value in aerofoilData:
        # if value[2] >= 0.0:
            # aerofoilData_top.append(value)
        # else:
            # aerofoilData_bottom.append(value)

    for value in aerofoilData:
		if (not(1e-07 > value[2] > -1e-07 and value[0] > 0.2)):
			if value[2] >= 0.0:
				aerofoilData_top.append(value)
			else:
				aerofoilData_bottom.append(value)		
			
    aerofoilData_top.sort(key=lambda x: x[0], reverse=True)
    aerofoilData_bottom.sort(key=lambda x: x[0])
    
    #-----------------------------------------------------
    # find trailing edge displacement of bottom node of TE gap
    displTE = aerofoilData_bottom[-1][3]
    displTEx = aerofoilData_bottom[-1][1]

    # find leading edge displacement of bottom node of LE
    displLE = aerofoilData_top[-1][3]
    displLEx = aerofoilData_top[-1][1]
	
    # twist
    twist = numpy.arctan((displTE-displLE)/(chordLength-displLEx+displTEx))*180/pi

    aerofoilDisplacedX = []
    aerofoilDisplacedZ = []

    for line in aerofoilData_top:
        aerofoilDisplacedX.append( line[0]+line[1] )
        aerofoilDisplacedZ.append( line[2]+line[3] )
    for line in aerofoilData_bottom:
        aerofoilDisplacedX.append( line[0]+line[1] )
        aerofoilDisplacedZ.append( line[2]+line[3] )

    #-----------------------------------------------------
    # normalize to unit chord
    aerofoilDisplaced_norm = [map(lambda x: x/chordLength, p) for p in zip(aerofoilDisplacedX,aerofoilDisplacedZ)]   

    #-----------------------------------------------------
    # write *.dat file with coordinates of deformed aerofoil for XFOIL input
    coord = file(os.path.join(pathXfoilData,'coords_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat'),'w')
    coord.write(jobName+'     end of '+stepName+'\n')
    for line in aerofoilDisplaced_norm:
        coord.write('     '.join(str(round(x,7)) for x in line) + '\n')
    coord.close()
    
    return twist, displLE, displTE, zip(aerofoilDisplacedX,aerofoilDisplacedZ)

def coordsDeformedAerofoilCAMBER(pathXfoilData,jobName,jobNumber,jobIndex,chordLength,dpos,xCL,yCL):

    # open output database
    odbName = jobName + '.odb'
    odb = visualization.openOdb(odbName)

    stepName = 'Step-%s' %jobIndex
	
    setName = 'SECTION%s' %dpos

    #-----------------------------------------------------
    # access relevant node sets
    # aerofoilY0 = odb.rootAssembly.instances['AEROFOILINSTANCE-1'].nodeSets['SETAEROFOILOUTEDGEY01']
    aerofoilY0 = odb.rootAssembly.instances['WING-FINAL-1'].nodeSets[setName]
	
    aerofoilCoordsX = []
    aerofoilCoordsZ = []

    # nodesAerofoil=aerofoilY0.nodes
    for k in aerofoilY0.nodes:
        aerofoilCoordsX.append(k.coordinates[0])
        aerofoilCoordsZ.append(k.coordinates[1])
          
    displAerofoilX = []
    displAerofoilZ = []

    dispField = odb.steps[stepName].frames[-1].fieldOutputs['U']
    dispSubField = dispField.getSubset(region=aerofoilY0)
    displValues = dispSubField.values
    for node in displValues:
        displAerofoilX.append(node.data[0])
        displAerofoilZ.append(node.data[1])
        
    odb.close()
    #-----------------------------------------------------
        
    aerofoilData = zip(aerofoilCoordsX,displAerofoilX,aerofoilCoordsZ,displAerofoilZ)

    aerofoilData_top = []
    aerofoilData_bottom = []

    # for value in aerofoilData:
        # if value[2] >= 0.0:
            # aerofoilData_top.append(value)
        # else:
            # aerofoilData_bottom.append(value)

    for value in aerofoilData:
		#if (not(1e-07 > value[2] > -1e-07 and value[0] > 0.2)):
			#if value[2] >= 0.0:
			#	aerofoilData_top.append(value)
			#else:
			#	aerofoilData_bottom.append(value)
        if (not(1e-07 > value[2] > -1e-07 and value[0] > 0.1)):
            if numpy.interp(value[0], xCL, yCL) <= value[2]:
            	aerofoilData_top.append(value)
            else:
            	aerofoilData_bottom.append(value)				
			
    aerofoilData_top.sort(key=lambda x: x[0], reverse=True)
    aerofoilData_bottom.sort(key=lambda x: x[0])
    
    #-----------------------------------------------------
    # find trailing edge displacement of bottom node of TE gap
    displTE = aerofoilData_bottom[-1][3]
    displTEx = aerofoilData_bottom[-1][1]

    # find leading edge displacement of bottom node of LE
    displLE = aerofoilData_top[-1][3]
    displLEx = aerofoilData_top[-1][1]
	
    # twist
    twist = numpy.arctan((displTE-displLE)/(chordLength-displLEx+displTEx))*180/pi

    aerofoilDisplacedX = []
    aerofoilDisplacedZ = []

    for line in aerofoilData_top:
        aerofoilDisplacedX.append( line[0]+line[1] )
        aerofoilDisplacedZ.append( line[2]+line[3] )
    for line in aerofoilData_bottom:
        aerofoilDisplacedX.append( line[0]+line[1] )
        aerofoilDisplacedZ.append( line[2]+line[3] )

    #-----------------------------------------------------
    # normalize to unit chord
    aerofoilDisplaced_norm = [map(lambda x: x/chordLength, p) for p in zip(aerofoilDisplacedX,aerofoilDisplacedZ)]   

    #-----------------------------------------------------
    # write *.dat file with coordinates of deformed aerofoil for XFOIL input
    coord = file(os.path.join(pathXfoilData,'coords_'+str(jobNumber)+'-'+str(jobIndex)+'-'+str(dpos)+'.dat'),'w')
    coord.write(jobName+'     end of '+stepName+'\n')
    for line in aerofoilDisplaced_norm:
        coord.write('     '.join(str(round(x,7)) for x in line) + '\n')
    coord.close()
    
    return twist, displLE, displTE, zip(aerofoilDisplacedX,aerofoilDisplacedZ)

#######################################################################
# MAX STRESS

def maxSTRESS(jobName,jobIndex):
	odbName = jobName + '.odb'
	odb = visualization.openOdb(odbName)

	stepName = 'Step-%s' %jobIndex

	mStr=0
	frameValueOld=-0.05
	for iterations in odb.steps[stepName].frames:
		if iterations.frameValue > frameValueOld+0.049:
			#print(frameValueOld)
			frameValueOld=frameValueOld+0.049

			mStress = iterations.fieldOutputs['MSTRS']
			#mStressTsaiW= iterations.fieldOutputs['TSAIW']
			#mStrain= iterations.fieldOutputs['MSTRN']
			
			#mStr=0
			#mTsaiW=0
			#meps=0
			for it in xrange(1,len(mStress.values)):
				if mStress.values[it].data > mStr:
					mStr=mStress.values[it].data
				#if mStressTsaiW.values[it].data > mTsaiW:
				#    mTsaiW=mStressTsaiW.values[it].data
				#if mStrain.values[it].data > meps:
				#    meps=mStressTsaiW.values[it].data
	odb.close()

	return mStr
	
	
#######################################################################
# INPUT FILE MODIFICATION

def GetBlockPosition(modelName, blockPrefix):
    if blockPrefix == '':
        return len(mdb.models[modelName].keywordBlock.sieBlocks)-1
    pos = 0
    for block in mdb.models[modelName].keywordBlock.sieBlocks:
        if block[0:len(blockPrefix)].lower()==blockPrefix.lower():
            return pos
        pos=pos+1
    return -1
	

