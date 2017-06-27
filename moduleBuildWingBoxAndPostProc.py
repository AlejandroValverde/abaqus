#Functions for abaqus

#Import all modules here so that they are accessible by all the upcoming functions
import math
import numpy as np
from shutil import copyfile
import os, re
import platform
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
import regionToolset
import copy

#For the post-processing
from abaqus import *
from abaqusConstants import *

#User-defined modules
from moduleCommon import *

def loadParameters(paraRead, fileName):

	file = open(fileName, 'r')

	lines = file.readlines()

	for i in range(int((len(lines))/2)):

		nameParater = lines[(i*2)]
		valueParater = lines[(2*i)+1]

		if platform.system() == 'Linux':

			valueParater = valueParater.replace('\r\n','')
			nameParater = nameParater.replace('\r\n','')			

		elif platform.system() == 'Windows':

			valueParater = valueParater.replace('\n','')
			nameParater = nameParater.replace('\n','')

		else:

			valueParater = valueParater.replace('\r\n','')
			nameParater = nameParater.replace('\r\n','')

			print('OS not recognized, assumed unix based')

		setattr(paraRead, nameParater, valueParater)

	file.close()

	return paraRead

def loadMaterials(model, design, load, mat):

	model.Material(name='ABS')
	model.materials['ABS'].Elastic(table=((mat.E_chiral, mat.v_chiral), ))
	model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
	    integrationRule=SIMPSON, material='ABS', name='Section-ABS', numIntPts=5, 
	    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
	    thickness=design.tChiral, thicknessField='', thicknessModulus=None, thicknessType=
	    UNIFORM, useDensity=OFF)

	model.Material(name='Alu')
	model.materials['Alu'].Elastic(table=((mat.E1, mat.v1), ))
	model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
	    integrationRule=SIMPSON, material='Alu', name='Section-Alu', numIntPts=5, 
	    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
	    thickness=2, thicknessField='', thicknessModulus=None, thicknessType=
	    UNIFORM, useDensity=OFF)
	model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
	    integrationRule=SIMPSON, material='ABS', name='Section-Ellipse', numIntPts=5, 
	    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
	    thickness=design.tChiral, thicknessField='', thicknessModulus=None, thicknessType=
	    UNIFORM, useDensity=OFF)

	#For the ribs
	model.Material(name='Alu_rib')
	model.materials['Alu_rib'].Elastic(table=((mat.E_rib, mat.v_rib), ))

	if load.typeAbaqus == 'Explicit': #Necessary to define material density when computing Abaqus Explicit

		model.materials['Alu'].Density(table=((2.7E-6, ), )) #kg /mm^3
		model.materials['ABS'].Density(table=((1.07E-6, ), )) #kg /mm^3
		model.materials['Alu_rib'].Density(table=((2.7E-6, ), )) #kg /mm^3

def buildBasicChiral(model, design):

	## diese parameter funktionieren

	r = design.r  #ringradius
	B = design.B  #ringtiefe
	L = design.L #halblaenge
	e=design.eOverB*B

	M=8 # anzahl zellen nach oben - not useful in this script
	N=8 # anzahl zellen nach rechts - not useful in this script

	#Build sketches
	model.ConstrainedSketch(name='__profile__', sheetSize=20.0)
	model.sketches['__profile__'].CircleByCenterPerimeter(center=(
	    0.0, 0.0), point1=(0.0, r))

	x2center=2*L
	y2center=2*r
	x2Rand=2*L
	y2Rand=r

	model.sketches['__profile__'].CircleByCenterPerimeter(center=(
	    x2center, y2center), point1=(x2Rand, y2Rand))

	math.radians(60)

	x3center=2*L-cos(math.radians(30))*r-cos(math.radians(60))*2*L-cos(math.radians(30))*r
	y3center=2*r-sin(math.radians(30))*r+sin(math.radians(60))*2*L-sin(math.radians(30))*r

	##x3Rand=200-cos(0.523)*r-cos(1.04718)*2*L
	##y3Rand=2*r-sin(0.523)*r+sin(1.04718)*2*L
	##
	##model.sketches['__profile__'].CircleByCenterPerimeter(center=(
	##    x3center, y3center), point1=(x3Rand, y3Rand))
	##
	model.Part(dimensionality=THREE_D, name='Part-1', type=
	    DEFORMABLE_BODY)
	model.parts['Part-1'].BaseShellExtrude(depth=B, sketch=
	    model.sketches['__profile__'])
	del model.sketches['__profile__']

	vertplane=model.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=L, principalPlane=YZPLANE)
	cyl0=model.parts['Part-1'].DatumAxisByCylFace(face=model.parts['Part-1'].faces[0])
	cyl1=model.parts['Part-1'].DatumAxisByCylFace(face=model.parts['Part-1'].faces[1])
	##cyl2=model.parts['Part-1'].DatumAxisByCylFace(face=model.parts['Part-1'].faces[2])

	#####tapespring left
	##
	model.ConstrainedSketch(gridSpacing=15.96, name='__profile__', 
	    sheetSize=638.64, transform=
	    model.parts['Part-1'].MakeSketchTransform(
	    sketchPlane=model.parts['Part-1'].datums[2], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['Part-1'].datums[4], 
	    sketchOrientation=LEFT, origin=(100.0, 0.0, 0.0)))
	model.parts['Part-1'].projectReferencesOntoSketch(filter=
	    COPLANAR_EDGES, sketch=model.sketches['__profile__'])
	if e==0:
	    model.sketches['__profile__'].Line(point1=(r, 0.0), point2=(r, B))
	else:
	    model.sketches['__profile__'].Arc3Points(point1=(r, 0.0), point2=(r, B), point3=(r-e, B/2))

	model.parts['Part-1'].ShellExtrude(flipExtrudeDirection=ON, 
	    sketch=model.sketches['__profile__'], sketchOrientation=
	    LEFT, sketchPlane=model.parts['Part-1'].datums[2], 
	    sketchPlaneSide=SIDE1, sketchUpEdge=
	    model.parts['Part-1'].datums[4], upToFace=
	    model.parts['Part-1'].faces[1])
	del model.sketches['__profile__']

	###tapespring right
	##
	model.ConstrainedSketch(gridSpacing=15.96, name='__profile__', 
	    sheetSize=638.64, transform=
	    model.parts['Part-1'].MakeSketchTransform(
	    sketchPlane=model.parts['Part-1'].datums[2], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['Part-1'].datums[4], 
	    sketchOrientation=LEFT, origin=(100.0, 0.0, 0.0)))
	model.parts['Part-1'].projectReferencesOntoSketch(filter=
	    COPLANAR_EDGES, sketch=model.sketches['__profile__'])

	if e==0:
	    model.sketches['__profile__'].Line(point1=(r, 0.0), point2=(r, B))
	else:
	    model.sketches['__profile__'].Arc3Points(point1=(r, 0.0), point2=(r, B), point3=(r+e, B/2))

	model.parts['Part-1'].ShellExtrude(flipExtrudeDirection=OFF, 
	    sketch=model.sketches['__profile__'], sketchOrientation=
	    LEFT, sketchPlane=model.parts['Part-1'].datums[2], 
	    sketchPlaneSide=SIDE1, sketchUpEdge=
	    model.parts['Part-1'].datums[4], upToFace=
	    model.parts['Part-1'].faces[0])
	del model.sketches['__profile__']
	#
	model.parts['Part-1'].DatumPlaneByPrincipalPlane(offset=r, 
	    principalPlane=XZPLANE)

	#
	#####close tapesprings 
	##
	if e!= 0:
	    model.ConstrainedSketch(gridSpacing=16.74, name='__profile__', 
	        sheetSize=669.72, transform=
	        model.parts['Part-1'].MakeSketchTransform(
	        sketchPlane=model.parts['Part-1'].datums[7], 
	        sketchPlaneSide=SIDE1, 
	        sketchUpEdge=model.parts['Part-1'].edges[2], 
	        sketchOrientation=RIGHT, origin=(0.0, 10.0, 0.0)))
	    model.parts['Part-1'].projectReferencesOntoSketch(filter=
	        COPLANAR_EDGES, sketch=model.sketches['__profile__'])
	    model.sketches['__profile__'].Line(point1=(L, 0.0), point2=
	        (L, -B))
	    model.parts['Part-1'].ShellExtrude(flipExtrudeDirection=OFF, 
	        sketch=model.sketches['__profile__'], sketchOrientation=
	        RIGHT, sketchPlane=model.parts['Part-1'].datums[7], 
	        sketchPlaneSide=SIDE1, sketchUpEdge=
	        model.parts['Part-1'].edges[2], upToFace=
	        model.parts['Part-1'].faces[0])
	    del model.sketches['__profile__']
	    ##
	    #####close tapespring unten
	    ##
	    model.ConstrainedSketch(gridSpacing=16.74, name='__profile__', 
	        sheetSize=669.72, transform=
	        model.parts['Part-1'].MakeSketchTransform(
	        sketchPlane=model.parts['Part-1'].datums[7], 
	        sketchPlaneSide=SIDE1, 
	        sketchUpEdge=model.parts['Part-1'].edges[0], 
	        sketchOrientation=RIGHT, origin=(0.0, 10.0, 0.0)))
	    model.parts['Part-1'].projectReferencesOntoSketch(filter=
	        COPLANAR_EDGES, sketch=model.sketches['__profile__'])
	    model.sketches['__profile__'].Line(point1=(L, 0.0), point2=
	        (L, -B))
	    model.parts['Part-1'].ShellExtrude(flipExtrudeDirection=ON, 
	        sketch=model.sketches['__profile__'], sketchOrientation=
	        RIGHT, sketchPlane=model.parts['Part-1'].datums[7], 
	        sketchPlaneSide=SIDE1, sketchUpEdge=
	        model.parts['Part-1'].edges[0], upToFace=
	        model.parts['Part-1'].faces[2])
	    del model.sketches['__profile__']
	    ##

	### untere tapespring fertig

	### Halbligament erzeugen

	model.Part(name='Part-Halbligament', objectToCopy=model.parts['Part-1'])
	plane_halb=model.parts['Part-Halbligament'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XYPLANE)
	model.ConstrainedSketch(gridSpacing=7.18, name='__profile__', 
	    sheetSize=287.44, transform=
	    model.parts['Part-Halbligament'].MakeSketchTransform(
	    sketchPlane=model.parts['Part-Halbligament'].datums[plane_halb.id], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['Part-Halbligament'].edges[12], 
	    sketchOrientation=RIGHT, origin=(50.0, 10.0, 0.0)))
	model.parts['Part-Halbligament'].projectReferencesOntoSketch(filter=
	    COPLANAR_EDGES, sketch=model.sketches['__profile__'])
	model.sketches['__profile__'].Line(point1=(-r, r), point2=(0.0, 0.0))
	model.sketches['__profile__'].Line(point1=(0,0), point2=(r,-r))
	model.sketches['__profile__'].Line(point1=(r,-r), point2=(2*L,-r))
	model.sketches['__profile__'].Line(point1=(2*L,-r), point2=(2*L,3*r))
	model.sketches['__profile__'].Line(point1=(2*L,3*r), point2=(-r,3*r))
	model.sketches['__profile__'].Line(point1=(-r,3*r), point2=(-r,r))

	model.parts['Part-Halbligament'].CutExtrude(flipExtrudeDirection=ON, 
	    sketch=model.sketches['__profile__'], sketchOrientation=
	    RIGHT, sketchPlane=model.parts['Part-Halbligament'].datums[plane_halb.id], 
	    sketchPlaneSide=SIDE1, sketchUpEdge=model.parts['Part-Halbligament'].edges[12])
	del model.sketches['__profile__']

def internalParameters(model, design):

	# Parameters from the basic chiral structure
	r = design.r  #ringradius
	B = design.B  #ringtiefe
	L = design.L #halblaenge
	e=0.01*B

	# Unit Cell

	rotangle = 60
	drehwinkel=math.degrees(atan(r/L))
	AbstandMittelpunkte=math.sqrt( math.pow((2*r),2) + math.pow((2*L),2)) #Distance center points
	HoeheDreieck=math.sqrt( 0.75*math.pow((AbstandMittelpunkte),2) ) #Height triangle

	# Print out final dimensions for the lattice
	design.L1 = AbstandMittelpunkte*design.N #Lattice total length in the spanwise direction
	design.L2 = 2*HoeheDreieck*design.M #Lattice total length in the transversal direction

	#Extract important parameters
	design.heightTriangle = HoeheDreieck
	design.distanceCenterPoints = AbstandMittelpunkte

	return design

def buildLattice(model, design):
	

	# Parameters from the basic chiral structure
	r = design.r  #ringradius
	B = design.B  #ringtiefe
	L = design.L #halblaenge
	e=0.01*B

	# Unit Cell

	rotangle = 60
	drehwinkel=math.degrees(atan(r/L))
	AbstandMittelpunkte=math.sqrt( math.pow((2*r),2) + math.pow((2*L),2)) #Distance center points
	HoeheDreieck=math.sqrt( 0.75*math.pow((AbstandMittelpunkte),2) ) #Height triangle

	print('Lattice total length in the spanwise direction: ' + str(design.L1))
	print('Lattice total length in the transversal direction: ' + str(design.L2))

	#Rotate basic structure and obtain an initial approach to the unit cell
	for i in range(1,7):
	    model.rootAssembly.Instance(dependent=ON, name='Ligament-'+str(i), part=model.parts['Part-1'])
	    model.rootAssembly.rotate(angle=drehwinkel, axisDirection=(0.0, 0.0, -1.0), axisPoint=(0, 0, 0.0), instanceList=('Ligament-'+str(i), ))
	    model.rootAssembly.rotate(angle=i*60.0, axisDirection=(0.0, 0.0, -1), axisPoint=(AbstandMittelpunkte ,0, 0.0), instanceList=('Ligament-'+str(i), ))
	    model.rootAssembly.Instance(dependent=ON, name='Ligament-out-'+str(i), part=model.parts['Part-1'])
	    model.rootAssembly.rotate(angle=drehwinkel, axisDirection=(0.0, 0.0, -1.0), axisPoint=(0, 0, 0.0), instanceList=('Ligament-out-'+str(i), ))

	model.rootAssembly.translate(instanceList=('Ligament-out-1', ), vector=((AbstandMittelpunkte/2, -HoeheDreieck,0.0)))
	model.rootAssembly.translate(instanceList=('Ligament-out-2', ), vector=((AbstandMittelpunkte/2, HoeheDreieck,0.0)))

	model.rootAssembly.rotate(angle=-60, axisDirection=(0.0, 0.0, -1), axisPoint=(0 ,0, 0.0), instanceList=('Ligament-out-3', ))
	model.rootAssembly.rotate(angle=60, axisDirection=(0.0, 0.0, -1), axisPoint=(0 ,0, 0.0), instanceList=('Ligament-out-4', ))

	model.rootAssembly.rotate(angle=-60, axisDirection=(0.0, 0.0, -1), axisPoint=(AbstandMittelpunkte ,0, 0.0), instanceList=('Ligament-out-5', ))
	model.rootAssembly.rotate(angle=60, axisDirection=(0.0, 0.0, -1), axisPoint=(AbstandMittelpunkte ,0, 0.0), instanceList=('Ligament-out-6', ))
	model.rootAssembly.translate(instanceList=('Ligament-out-5', ), vector=((AbstandMittelpunkte, 0,0.0)))
	model.rootAssembly.translate(instanceList=('Ligament-out-6', ), vector=((AbstandMittelpunkte, 0.0,0.0)))


	model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
	    instances=(model.rootAssembly.instances['Ligament-1'], 
	    model.rootAssembly.instances['Ligament-out-1'], 
	    model.rootAssembly.instances['Ligament-2'], 
	    model.rootAssembly.instances['Ligament-out-2'], 
	    model.rootAssembly.instances['Ligament-3'], 
	    model.rootAssembly.instances['Ligament-out-3'], 
	    model.rootAssembly.instances['Ligament-4'], 
	    model.rootAssembly.instances['Ligament-out-4'], 
	    model.rootAssembly.instances['Ligament-5'], 
	    model.rootAssembly.instances['Ligament-out-5'], 
	    model.rootAssembly.instances['Ligament-6'], 
	    model.rootAssembly.instances['Ligament-out-6']), name=
	    'UnitCell', originalInstances=SUPPRESS)


	# Locals COSY

	##LocCYL=model.rootAssembly.DatumCsysByThreePoints(coordSysType=
	##    CYLINDRICAL, name='Datum-CYL', origin=(AbstandMittelpunkte, 0, 0.0), point1=(AbstandMittelpunkte, -1, 0.0), point2=(AbstandMittelpunkte+1, 0, 0.0))

	# hier sollte es geschnitten werden!

	cuttingplane=model.parts['UnitCell'].DatumPlaneByPrincipalPlane(offset=0.0, 
	    principalPlane=XYPLANE)
	cuttingedge=model.parts['UnitCell'].DatumAxisByPrincipalAxis(principalAxis= YAXIS)

	model.ConstrainedSketch(gridSpacing=24.11, name='__profile__', 
	    sheetSize=964.7, transform=
	    model.parts['UnitCell'].MakeSketchTransform(
	    sketchPlane=model.parts['UnitCell'].datums[cuttingplane.id], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['UnitCell'].datums[cuttingedge.id], 
	    sketchOrientation=LEFT, origin=(AbstandMittelpunkte, 0.0, 0.0)))

	# contruction lines

	cline1=model.sketches['__profile__'].ConstructionLine(angle=30.0, point1=(-AbstandMittelpunkte/2, HoeheDreieck))
	cline2=model.sketches['__profile__'].ConstructionLine(angle=-30.0, point1=(AbstandMittelpunkte/2, HoeheDreieck))
	cline3=model.sketches['__profile__'].ConstructionLine(angle=90.0, point1=(AbstandMittelpunkte, 0.0))
	cline4=model.sketches['__profile__'].ConstructionLine(angle=30.0, point1=(AbstandMittelpunkte/2, -HoeheDreieck))
	cline5=model.sketches['__profile__'].ConstructionLine(angle=-30.0, point1=(-AbstandMittelpunkte/2, -HoeheDreieck))

	# Ligamentpositionen in Polarkoordinaten

	RMidLig=AbstandMittelpunkte/(2*tan(math.radians(30)))                   # Mittelpunkt Ligament
	winkelConLig=360/(2*math.pi)*tan(r/AbstandMittelpunkte)                 # Winkel zu Ligamentschnittpunkten
	RConLig=(AbstandMittelpunkte-sin(math.radians(drehwinkel))*r)/cos(math.radians(winkelConLig))   # Verbingungspunkt Ligament/Kreis

	winkelConLigRest=120-(90-drehwinkel)
	breiteKreisDreieck=AbstandMittelpunkte-r*cos(math.radians(winkelConLigRest))
	hoeheKreisDreieck=r*sin(math.radians(winkelConLigRest))

	RConLigOther=math.sqrt( math.pow((breiteKreisDreieck),2) + math.pow((hoeheKreisDreieck),2))
	winkelConLigOther=360/(2*math.pi)*tan(hoeheKreisDreieck/breiteKreisDreieck)                 # Winkel zu Ligamentschnittpunkten

	# Ligamente Innenseite (verbunden zu MidCircle)

	xRadius=AbstandMittelpunkte-(r*cos(math.radians(120-drehwinkel-90)))
	yRadius=r*sin(math.radians(120-drehwinkel-90))
	RadiusLigClose=math.sqrt( math.pow((xRadius),2) + math.pow((yRadius),2))

	winkelLigClose=math.degrees(atan(yRadius/xRadius))

	# Ligamentpositionen in karthesischen Koordinaten

	xMidLig = []
	yMidLig = []

	xMidLig_aft = []
	yMidLig_aft = []

	xRadius=AbstandMittelpunkte-(r*cos(math.radians(120-drehwinkel-90)))
	yRadius=r*sin(math.radians(120-drehwinkel-90))
	RadiusLigClose=math.sqrt( math.pow((xRadius),2) + math.pow((yRadius),2))

	xMidLig_front = []
	yMidLig_front = []

	xConLig = []                                                                # das sind die karthesischen Koordinaten zum SchnittpunktLigament - Kreis
	yConLig = []

	xConLigOther = []                                                           # andere Seite der Kreisverbindung
	yConLigOther = []

	CircleMx = []                                                           # andere Seite der Kreisverbindung
	CircleMy = []

	##CircleMx_assembly = []                                                           # andere Seite der Kreisverbindung
	##CircleMy_assembly = []

	LigHalfx = []                                                           # andere Seite der Kreisverbindung
	LigHalfy = []

	xLigClose = []                                                           # andere Seite der Kreisverbindung
	yLigClose = []

	xLigHalf_front =[]
	yLigHalf_front =[]

	xLigHalf_aft =[]
	yLigHalf_aft =[]

	winkel1=120-drehwinkel-90
	    
	for i in range(1,7):

	    xcur=RMidLig*cos(math.radians(30+(i-1)*60))
	    ycur=RMidLig*sin(math.radians(30+(i-1)*60))
	    xMidLig.append(xcur)
	    yMidLig.append(ycur)
	    
	    xcur_aft=0.9*RMidLig*cos(math.radians(25+(i-1)*60))
	    ycur_aft=0.9*RMidLig*sin(math.radians(25+(i-1)*60))
	    xMidLig_aft.append(xcur_aft)
	    yMidLig_aft.append(ycur_aft)
	    
	    xcur_front=1.1*RMidLig*cos(math.radians(35+(i-1)*60))
	    ycur_front=1.1*RMidLig*sin(math.radians(35+(i-1)*60))
	    xMidLig_front.append(xcur_front)
	    yMidLig_front.append(ycur_front)
	    
	    xConLigCur=RConLig*cos(math.radians(-winkelConLig+(i-1)*60))
	    yConLigCur=RConLig*sin(math.radians(-winkelConLig+(i-1)*60))
	    xConLig.append(xConLigCur)
	    yConLig.append(yConLigCur)
	    
	    xConLigOtherCur=RConLigOther*cos(math.radians(winkelConLigOther+(i-1)*60))
	    yConLigOtherCur=RConLigOther*sin(math.radians(winkelConLigOther+(i-1)*60))
	    xConLigOther.append(xConLigOtherCur)
	    yConLigOther.append(yConLigOtherCur)
	    
	    CircleMxCur=AbstandMittelpunkte*cos(math.radians((i-1)*60))
	    CircleMyCur=AbstandMittelpunkte*sin(math.radians((i-1)*60))
	##    CircleMx_assemblyCur=AbstandMittelpunkte*cos(math.radians((i-1)*60))+AbstandMittelpunkte
	##    CircleMy_assemblyCur=AbstandMittelpunkte*sin(math.radians((i-1)*60))
	    CircleMx.append(CircleMxCur)
	    CircleMy.append(CircleMyCur)
	##    CircleMx_assembly.append(CircleMx_assemblyCur)
	##    CircleMy_assembly.append(CircleMy_assemblyCur)
	    
	    LigHalfxCur=AbstandMittelpunkte/2*cos(math.radians((i-1)*60))
	    LigHalfyCur=AbstandMittelpunkte/2*sin(math.radians((i-1)*60))
	    LigHalfx.append(LigHalfxCur)
	    LigHalfy.append(LigHalfyCur)
	    
	    LigHalfx_aftCur=0.9*AbstandMittelpunkte/2*cos(math.radians(10+(i-1)*60))
	    LigHalfy_aftCur=0.9*AbstandMittelpunkte/2*sin(math.radians(10+(i-1)*60))
	    xLigHalf_aft.append(LigHalfx_aftCur)
	    yLigHalf_aft.append(LigHalfy_aftCur)
	    
	    LigHalfx_frontCur=1.1*AbstandMittelpunkte/2*cos(math.radians(-10+(i-1)*60))
	    LigHalfy_frontCur=1.1*AbstandMittelpunkte/2*sin(math.radians(-10+(i-1)*60))
	    xLigHalf_front.append(LigHalfx_frontCur)
	    yLigHalf_front.append(LigHalfy_frontCur)
	    
	    LigClosexCur=RadiusLigClose*cos(math.radians(winkelLigClose-(i-1)*60))
	    LigCloseyCur=-RadiusLigClose*sin(math.radians(winkelLigClose-(i-1)*60))
	    xLigClose.append(LigClosexCur)
	    yLigClose.append(LigCloseyCur)


	##model.ConstrainedSketch(gridSpacing=24.11, name='__profile__', 
	##    sheetSize=964.7, transform=
	##    model.parts['UnitCell'].MakeSketchTransform(
	##    sketchPlane=model.parts['UnitCell'].datums[cuttingplane.id], 
	##    sketchPlaneSide=SIDE1, 
	##    sketchUpEdge=model.parts['UnitCell'].datums[cuttingedge.id], 
	##    sketchOrientation=LEFT, origin=(AbstandMittelpunkte, 0.0, 0.0)))

	model.sketches['__profile__'].Line(point1=(xMidLig_front[5],yMidLig_front[5]), point2=(xLigClose[4]+AbstandMittelpunkte,yLigClose[4]))
	model.sketches['__profile__'].Line(point1=(xLigClose[4]+AbstandMittelpunkte,yLigClose[4]), point2=(CircleMx[5], CircleMy[5]))
	model.sketches['__profile__'].Line(point1=(CircleMx[5], CircleMy[5]), point2=(xConLig[5], yConLig[5]))
	model.sketches['__profile__'].Line(point1=(xConLig[5], yConLig[5]), point2=(xMidLig_aft[4], yMidLig_aft[4]))
	model.sketches['__profile__'].Line(point1=(xMidLig_aft[4], yMidLig_aft[4]), point2=(xMidLig[4], yMidLig[4]))
	model.sketches['__profile__'].Line(point1=(xMidLig[4], yMidLig[4]), point2=(xMidLig_front[4], yMidLig_front[4]))
	model.sketches['__profile__'].Line(point1=(xMidLig_front[4], yMidLig_front[4]), point2=(xConLig[2], yConLig[2]-2*CircleMy[2]))
	model.sketches['__profile__'].Line(point1=(xConLig[2], yConLig[2]-2*CircleMy[2]), point2=(CircleMx[4],CircleMy[4]))
	model.sketches['__profile__'].Line(point1=(CircleMx[4],CircleMy[4]), point2=(xLigClose[4],yLigClose[4]))
	model.sketches['__profile__'].Line(point1=(xLigClose[4],yLigClose[4]), point2=(xLigHalf_aft[4],yLigHalf_aft[4]))
	model.sketches['__profile__'].Line(point1=(xLigHalf_aft[4],yLigHalf_aft[4]), point2=(LigHalfx[4],LigHalfy[4]))
	model.sketches['__profile__'].Line(point1=(LigHalfx[4],LigHalfy[4]), point2=(xLigHalf_front[4],yLigHalf_front[4]))
	model.sketches['__profile__'].Line(point1=(xLigHalf_front[4],yLigHalf_front[4]), point2=(LigHalfx[3],LigHalfy[3]))
	model.sketches['__profile__'].Line(point1=(LigHalfx[3],LigHalfy[3]), point2=(xMidLig_aft[2],yMidLig_aft[2]))
	model.sketches['__profile__'].Line(point1=(xMidLig_aft[2],yMidLig_aft[2]), point2=(xMidLig[2],yMidLig[2]))
	model.sketches['__profile__'].Line(point1=(xMidLig[2],yMidLig[2]), point2=(xMidLig_front[2],yMidLig_front[2]))
	model.sketches['__profile__'].Line(point1=(xMidLig_front[2],yMidLig_front[2]), point2=(xLigClose[1]-AbstandMittelpunkte,yLigClose[1]))
	model.sketches['__profile__'].Line(point1=(xLigClose[1]-AbstandMittelpunkte,yLigClose[1]), point2=(CircleMx[2],CircleMy[2]))
	model.sketches['__profile__'].Line(point1=(CircleMx[2],CircleMy[2]), point2=(xConLig[2],yConLig[2]))
	model.sketches['__profile__'].Line(point1=(xConLig[2],yConLig[2]), point2=(xMidLig_aft[1],yMidLig_aft[1]))
	model.sketches['__profile__'].Line(point1=(xMidLig_aft[1],yMidLig_aft[1]), point2=(xMidLig[1],yMidLig[1]))
	model.sketches['__profile__'].Line(point1=(xMidLig[1],yMidLig[1]), point2=(xMidLig_front[1],yMidLig_front[1]))
	model.sketches['__profile__'].Line(point1=(xMidLig_front[1],yMidLig_front[1]), point2=(xConLig[5],yConLig[5]+2*CircleMy[1]))
	model.sketches['__profile__'].Line(point1=(xConLig[5],yConLig[5]+2*CircleMy[1]), point2=(CircleMx[1],CircleMy[1]))
	model.sketches['__profile__'].Line(point1=(CircleMx[1],CircleMy[1]), point2=(xLigClose[1],yLigClose[1]))
	model.sketches['__profile__'].Line(point1=(xLigClose[1],yLigClose[1]), point2=(xLigHalf_aft[1],yLigHalf_aft[1]))
	model.sketches['__profile__'].Line(point1=(xLigHalf_aft[1],yLigHalf_aft[1]), point2=(LigHalfx[1],LigHalfy[1]))
	model.sketches['__profile__'].Line(point1=(LigHalfx[1],LigHalfy[1]), point2=(xLigHalf_front[1],yLigHalf_front[1]))
	# aussen
	model.sketches['__profile__'].Line(point1=(xLigHalf_front[1],yLigHalf_front[1]), point2=(xMidLig[0],yMidLig[0]))
	model.sketches['__profile__'].Line(point1=(xMidLig[0],yMidLig[0]), point2=(xMidLig_front[5],yLigHalf_front[1]))
	model.sketches['__profile__'].Line(point1=(xMidLig_front[5],yLigHalf_front[1]), point2=(xMidLig_front[5],CircleMy[1]+2*r))
	model.sketches['__profile__'].Line(point1=(xMidLig_front[5],CircleMy[1]+2*r), point2=(CircleMx[3]-2*r,CircleMy[1]+2*r))
	model.sketches['__profile__'].Line(point1=(CircleMx[3]-2*r,CircleMy[1]+2*r), point2=(CircleMx[3]-2*r,CircleMy[4]-2*r))
	model.sketches['__profile__'].Line(point1=(CircleMx[3]-2*r,CircleMy[4]-2*r), point2=(xMidLig_front[5],CircleMy[4]-2*r))
	model.sketches['__profile__'].Line(point1=(xMidLig_front[5],CircleMy[4]-2*r), point2=(xMidLig_front[5],yMidLig_front[5]))

	model.parts['UnitCell'].CutExtrude(flipExtrudeDirection=ON, 
	    sketch=model.sketches['__profile__'], sketchOrientation=
	    LEFT, sketchPlane=model.parts['UnitCell'].datums[cuttingplane.id], 
	    sketchPlaneSide=SIDE1, sketchUpEdge=
	    model.parts['UnitCell'].datums[cuttingedge.id])
	del model.sketches['__profile__']


	# zweitextrude:

	model.ConstrainedSketch(gridSpacing=24.11, name='__profile__', 
	    sheetSize=964.7, transform=
	    model.parts['UnitCell'].MakeSketchTransform(
	    sketchPlane=model.parts['UnitCell'].datums[cuttingplane.id], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['UnitCell'].datums[cuttingedge.id], 
	    sketchOrientation=LEFT, origin=(AbstandMittelpunkte, 0.0, 0.0)))

	model.sketches['__profile__'].Line(point1=(xLigHalf_front[1],yLigHalf_front[1]), point2=(LigHalfx[0],LigHalfy[0]))
	model.sketches['__profile__'].Line(point1=(LigHalfx[0],LigHalfy[0]), point2=(xMidLig_aft[5],yMidLig_aft[5]))
	model.sketches['__profile__'].Line(point1=(xMidLig_aft[5],yMidLig_aft[5]), point2=(xMidLig[5],yMidLig[5]))
	model.sketches['__profile__'].Line(point1=(xMidLig[5],yMidLig[5]), point2=(xMidLig_front[5],yMidLig_front[5]))
	model.sketches['__profile__'].Line(point1=(xMidLig_front[5],yMidLig_front[5]), point2=(CircleMx[0]+2*r,yMidLig_front[5]))
	model.sketches['__profile__'].Line(point1=(CircleMx[0]+2*r,yMidLig_front[5]), point2=(CircleMx[0]+2*r,yLigHalf_front[1]))
	model.sketches['__profile__'].Line(point1=(CircleMx[0]+2*r,yLigHalf_front[1]), point2=(xMidLig[0],yMidLig[0]))
	model.sketches['__profile__'].Line(point1=(xMidLig[0],yMidLig[0]), point2=(xLigHalf_front[1],yLigHalf_front[1]))


	model.parts['UnitCell'].CutExtrude(flipExtrudeDirection=ON, 
	    sketch=model.sketches['__profile__'], sketchOrientation=
	    LEFT, sketchPlane=model.parts['UnitCell'].datums[cuttingplane.id], 
	    sketchPlaneSide=SIDE1, sketchUpEdge=
	    model.parts['UnitCell'].datums[cuttingedge.id])
	del model.sketches['__profile__']

	############


	##########################################################################################################################################################

	#Create instances from parts

	M=design.M # cellsup
	for i in range(1,M):
	    model.rootAssembly.Instance(dependent=ON, name='UnitCell-'+str(i+1), part=model.parts['UnitCell'])
	    model.rootAssembly.instances['UnitCell-'+str(i+1)].translate(vector=(0.0, i*HoeheDreieck*2, 0.0))

	#Create tuple from instances for all the UnitCells in the Up direction
	cellsup=[]
	for i in xrange(1,M+1): #xrange is the same as range
	        cellsup.append(model.rootAssembly.instances['UnitCell-'+str(i)])
	cellsup=tuple(cellsup)

	#Merge instances in the Up direction
	model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, instances=( cellsup ), name='UPs', originalInstances=SUPPRESS)

	N=int(design.N) # cellsright
	for i in range(1,N):
	    model.rootAssembly.Instance(dependent=ON, name='UPs-'+str(i+1), part=model.parts['UPs'])
	    model.rootAssembly.instances['UPs-'+str(i+1)].translate(vector=(i*AbstandMittelpunkte, 0.0, 0.0))

	cellslow=[]
	for i in xrange(1,N+1):
	        cellslow.append(model.rootAssembly.instances['UPs-'+str(i)])

	cellslow=tuple(cellslow)

	#Create instance from all
	model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, instances=( cellslow ), name='All', originalInstances=SUPPRESS)

	#Create set for the part 'All' - Modified by Alex (getByBoundingBox)
	model.parts['All'].Set(faces=model.parts['All'].faces.getByBoundingBox(0, -100, -5, design.L1 + 200, design.L2 + 200, design.B + 20), name='LatticeAll')

	#Assign material and section to lattice - Modified by Alex (getByBoundingBox)
	model.parts['All'].SectionAssignment(offset=0.0, offsetField=''
	    , offsetType=MIDDLE_SURFACE, region=model.parts['All'].sets['LatticeAll'], sectionName='Section-ABS', 
	    thicknessAssignment=FROM_SECTION)
	########################################################################################################################################################
	
def cutLattice(model, design):

	#Create Datum plane
	datumPlane = model.parts['All'].DatumPlaneByPrincipalPlane(offset=0.0, 
	    principalPlane=XYPLANE)

	#Axis for sketch
	point1 = model.parts['All'].DatumPointByCoordinate(coords=(-1.0, design.cutUp, 
	    1.0))
	point2 = model.parts['All'].DatumPointByCoordinate(coords=(1.0, design.cutUp, 
	    1.0))
	datumAxis = model.parts['All'].DatumAxisByTwoPoint(point1=
	    model.parts['All'].datums[point1.id], point2=
	    model.parts['All'].datums[point2.id])

	#### Cutting
	margin = 200

	#Define cutting sketch - UP
	model.ConstrainedSketch(gridSpacing=75.3, name='__profile__', 
	    sheetSize=3012.29, transform=
	    model.parts['All'].MakeSketchTransform(
	    sketchPlane=model.parts['All'].datums[datumPlane.id], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['All'].datums[datumAxis.id], 
	    sketchOrientation=TOP, origin=(0.0, 0.0, 0.0)))
	model.sketches['__profile__'].rectangle(point1=(design.L1 + margin, design.cutUp), 
	    point2=(-margin, design.cutUp + margin))

	#Cut - UP
	model.parts['All'].CutExtrude(flipExtrudeDirection=ON, sketch=
	    model.sketches['__profile__'], sketchOrientation=TOP, 
	    sketchPlane=model.parts['All'].datums[datumPlane.id], sketchPlaneSide=
	    SIDE1, sketchUpEdge=model.parts['All'].datums[datumAxis.id])
	del model.sketches['__profile__']

	#Define cutting sketch - DOWN
	model.ConstrainedSketch(gridSpacing=75.3, name='__profile__', 
	    sheetSize=3012.29, transform=
	    model.parts['All'].MakeSketchTransform(
	    sketchPlane=model.parts['All'].datums[datumPlane.id], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['All'].datums[datumAxis.id], 
	    sketchOrientation=TOP, origin=(0.0, 0.0, 0.0)))
	model.sketches['__profile__'].rectangle(point1=(design.L1 + margin, design.cutDown), 
	    point2=(-margin, design.cutDown - margin))

	#Cut - DOWN
	model.parts['All'].CutExtrude(flipExtrudeDirection=ON, sketch=
	    model.sketches['__profile__'], sketchOrientation=TOP, 
	    sketchPlane=model.parts['All'].datums[datumPlane.id], sketchPlaneSide=
	    SIDE1, sketchUpEdge=model.parts['All'].datums[datumAxis.id])
	del model.sketches['__profile__']

	#Define cutting sketch - wingRoot
	model.ConstrainedSketch(gridSpacing=75.3, name='__profile__', 
	    sheetSize=3012.29, transform=
	    model.parts['All'].MakeSketchTransform(
	    sketchPlane=model.parts['All'].datums[datumPlane.id], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['All'].datums[datumAxis.id], 
	    sketchOrientation=TOP, origin=(0.0, 0.0, 0.0)))
	model.sketches['__profile__'].rectangle(point1=(design.cutWingRoot, design.cutDown - margin), 
	    point2=(design.cutWingRoot - margin, design.cutUp + margin))

	#Cut - wingRoot
	model.parts['All'].CutExtrude(flipExtrudeDirection=ON, sketch=
	    model.sketches['__profile__'], sketchOrientation=TOP, 
	    sketchPlane=model.parts['All'].datums[datumPlane.id], sketchPlaneSide=
	    SIDE1, sketchUpEdge=model.parts['All'].datums[datumAxis.id])
	del model.sketches['__profile__']

	#Define cutting sketch - wingTip
	model.ConstrainedSketch(gridSpacing=75.3, name='__profile__', 
	    sheetSize=3012.29, transform=
	    model.parts['All'].MakeSketchTransform(
	    sketchPlane=model.parts['All'].datums[datumPlane.id], 
	    sketchPlaneSide=SIDE1, 
	    sketchUpEdge=model.parts['All'].datums[datumAxis.id], 
	    sketchOrientation=TOP, origin=(0.0, 0.0, 0.0)))
	model.sketches['__profile__'].rectangle(point1=(design.cutWingTip, design.cutDown - margin), 
	    point2=(design.cutWingTip + margin, design.cutUp + margin))

	#Cut - wingTip
	model.parts['All'].CutExtrude(flipExtrudeDirection=ON, sketch=
	    model.sketches['__profile__'], sketchOrientation=TOP, 
	    sketchPlane=model.parts['All'].datums[datumPlane.id], sketchPlaneSide=
	    SIDE1, sketchUpEdge=model.parts['All'].datums[datumAxis.id])
	del model.sketches['__profile__']

def buildBox(model, design, mesh):

	#Build wing box around lattice

	#Input parameters
	# design.M ->Number of unit cells in transversal direction
	# design.N ->Number of unit cells in spanwise direction
	# design.r ->Node radius 
	# design.B ->Node depth
	# design.L ->half length

	#Calculate dimensions for the box
	design.C2 = (design.cutUp + abs(design.cutDown)) # + design.Ct #Interference on purpose = t 
	design.C1 = design.cutWingTip - design.cutWingRoot

	#Create sketch palette with dimensions 150% bigger and the maximum length in the lattice
	model.ConstrainedSketch(name='__profile__', sheetSize=1.5 * max(design.L1, design.L2))

	#Create C shape
	model.sketches['__profile__'].Line(point1=(design.C3, design.cutUp), 
	    point2=(0.0, design.cutUp))
	model.sketches['__profile__'].HorizontalConstraint(
	    addUndoState=False, entity=
	    model.sketches['__profile__'].geometry[2])
	model.sketches['__profile__'].Line(point1=(0.0, design.cutUp), 
	    point2=(0.0, design.cutDown))
	model.sketches['__profile__'].VerticalConstraint(addUndoState=
	    False, entity=model.sketches['__profile__'].geometry[3])
	model.sketches['__profile__'].PerpendicularConstraint(
	    addUndoState=False, entity1=
	    model.sketches['__profile__'].geometry[2], entity2=
	    model.sketches['__profile__'].geometry[3])
	model.sketches['__profile__'].Line(point1=(0.0, design.cutDown), 
	    point2=(design.C3, design.cutDown))
	model.sketches['__profile__'].HorizontalConstraint(
	    addUndoState=False, entity=
	    model.sketches['__profile__'].geometry[4])
	model.sketches['__profile__'].PerpendicularConstraint(
	    addUndoState=False, entity1=
	    model.sketches['__profile__'].geometry[3], entity2=
	    model.sketches['__profile__'].geometry[4])
	model.sketches['__profile__'].EqualLengthConstraint(entity1=
	    model.sketches['__profile__'].geometry[2], entity2=
	    model.sketches['__profile__'].geometry[4])

	#Shell extrude
	model.Part(dimensionality=THREE_D, name='C-box', type=
	    DEFORMABLE_BODY)
	model.parts['C-box'].BaseShellExtrude(depth=design.C1, sketch=
	    model.sketches['__profile__'])
	del model.sketches['__profile__']

	#Create variable from part


	#Create set from part
	model.parts['C-box'].Set(faces=
	    model.parts['C-box'].faces.findAt(((design.C3/2,design.cutUp,design.cutWingTip/2),), 
	    ((0.0,design.cutUp/2,design.cutWingTip/2),),
	    ((design.C3/2,design.cutDown,design.cutWingTip/2),),), name='C-box')

	#Create special set for meshing
	point1 = model.parts['C-box'].DatumPointByCoordinate(coords=(design.C3 - mesh.d, 
	    design.cutUp, 0.0))
	point2 = model.parts['C-box'].DatumPointByCoordinate(coords=(design.C3 - mesh.d, 
	    design.cutUp, design.C1))
	datumAxis = model.parts['C-box'].DatumAxisByTwoPoint(point1=
	    model.parts['C-box'].datums[point1.id], point2=
	    model.parts['C-box'].datums[point2.id])
	point3 = model.parts['C-box'].DatumPointByCoordinate(coords=(design.C3 - mesh.d, 
	    design.cutUp + 100, design.C1/2))
	datumPlane = model.parts['C-box'].DatumPlaneByLinePoint(line=
	    model.parts['C-box'].datums[datumAxis.id], point=
	    model.parts['C-box'].datums[point3.id])

	#Partition
	model.parts['C-box'].PartitionFaceByDatumPlane(datumPlane=
	    model.parts['C-box'].datums[datumPlane.id], faces=
	    model.parts['C-box'].faces.findAt(((design.C3/2,design.cutUp,design.cutWingTip/2),),
	    ((design.C3/2,design.cutDown,design.cutWingTip/2),),) )

	#Sets
	model.parts['C-box'].Set(faces=
	    model.parts['C-box'].faces.findAt(((design.C3 - 1,design.cutUp,design.cutWingTip/2),), 
	    ((design.C3 - 1,design.cutDown,design.cutWingTip/2),),), name='setCloseToLattice')

	model.parts['C-box'].Set(faces=
	    model.parts['C-box'].faces.findAt(((design.C3/2,design.cutUp,design.cutWingTip/2),), 
	    ((0.0,design.cutUp/2,design.cutWingTip/2),),
	    ((design.C3/2,design.cutDown,design.cutWingTip/2),),), name='notCloseToLattice')

	#Section
	model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
	    integrationRule=SIMPSON, material='Alu', name='Section-Cbox', numIntPts=5, 
	    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
	    thickness=design.Ct, thicknessField='', thicknessModulus=None, thicknessType=
	    UNIFORM, useDensity=OFF)
	model.parts['C-box'].SectionAssignment(offset=0.0, offsetField=
	    '', offsetType=MIDDLE_SURFACE, region=
	    model.parts['C-box'].sets['C-box'], sectionName=
	    'Section-Cbox', thicknessAssignment=FROM_SECTION)

	#### Instance operations ####
	model.rootAssembly.DatumCsysByDefault(CARTESIAN)
	model.rootAssembly.Instance(dependent=ON, name='C-box-inst', part=
	    model.parts['C-box'])
	model.rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 1.0, 
	    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('C-box-inst', ))
	model.rootAssembly.translate(instanceList=('C-box-inst', ), 
	    vector=(0.0, 0.0, design.C3))

	#Spanwise shifting
	model.rootAssembly.translate(instanceList=('C-box-inst', ), 
	    vector=(design.cutWingRoot, 0.0, 0.0))

	#Transversal shifting - Not necessary transversal shifting in this direction
	# model.rootAssembly.translate(instanceList=('C-box-inst', ), 
	#     vector=(0.0, -HoeheDreieck, 0.0))

def buildRib(model, design, typeOfRib, typeOfRib2):

	#build open Rib

	#Dimensions 
	design.rib1 = design.C3
	design.rib2 = design.C2 #Lattice dimension in the transversal direction
	# design.ribt must be defined as well

	#SKETCH DESIGN

	#Create sketch palette with dimensions 150% bigger and the maximum length in the lattice
	model.ConstrainedSketch(name='__profile__', sheetSize=1.5 * max(design.L1, design.L2))

	if typeOfRib == 'inner_ribs' or (typeOfRib == 'outer_ribs' and typeOfRib2 == 'open'):

		if typeOfRib == 'inner_ribs':

			design.rib1 = design.C3 - design.B - design.innerRibs_gap #Lattice dimension in the transversal direction

		elif typeOfRib == 'outer_ribs':

			design.rib1 = design.C3

		design.b = design.rib1 - (design.a * 3/2) #Internal parameter, the angle at the extreme of the rib is 45deg.

		#Drawing
		#
		#          ____ -----------------------
		#          ^   |          ^            |
		#          |   |          |a           |
		#          |   |          v            /  -> 45deg 
		#          |   |      -----------------
		#          |   |     |
		#          |   |     |
		#          |   |  a  |
		#    rib_1 |   |<--->|
		#          |   |     |       b
		#          |   |     |<-------------->|
		#          |   |     |                |
		#          |   |      -----------------
		#          |   |                       \  -> 45deg
		#          |   |                        |
		#          v   |                        |
		#          ____ ------------------------
		#              |                        |
		#              |         rib_1          |
		#              |<---------------------->|
		#
		#

		#Create shape
		model.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
		    design.rib1, 0.0))
		model.sketches['__profile__'].HorizontalConstraint(
		    addUndoState=False, entity=
		    model.sketches['__profile__'].geometry[2])
		model.sketches['__profile__'].Line(point1=(design.rib1, 0.0), point2=(
		    design.rib1, design.a/2))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
		    False, entity=model.sketches['__profile__'].geometry[3])
		model.sketches['__profile__'].PerpendicularConstraint(
		    addUndoState=False, entity1=
		    model.sketches['__profile__'].geometry[2], entity2=
		    model.sketches['__profile__'].geometry[3])
		model.sketches['__profile__'].Line(point1=(design.rib1, design.a/2), point2=
		    (design.b + design.a, design.a))
		model.sketches['__profile__'].Line(point1=(design.b + design.a, design.a), point2=
		    (design.a, design.a))
		model.sketches['__profile__'].HorizontalConstraint(
		    addUndoState=False, entity=
		    model.sketches['__profile__'].geometry[5])
		model.sketches['__profile__'].Line(point1=(design.a, design.a), point2=(
		    design.a, design.rib2 - design.a))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
		    False, entity=model.sketches['__profile__'].geometry[6])
		model.sketches['__profile__'].PerpendicularConstraint(
		    addUndoState=False, entity1=
		    model.sketches['__profile__'].geometry[5], entity2=
		    model.sketches['__profile__'].geometry[6])
		model.sketches['__profile__'].Line(point1=(design.a, design.rib2 - design.a), 
		    point2=(design.b + design.a, design.rib2 - design.a))
		model.sketches['__profile__'].HorizontalConstraint(
		    addUndoState=False, entity=
		    model.sketches['__profile__'].geometry[7])
		model.sketches['__profile__'].PerpendicularConstraint(
		    addUndoState=False, entity1=
		    model.sketches['__profile__'].geometry[6], entity2=
		    model.sketches['__profile__'].geometry[7])
		model.sketches['__profile__'].Line(point1=(design.b + design.a, design.rib2 - design.a), 
		    point2=(design.rib1, design.rib2 - (design.a/2)))
		model.sketches['__profile__'].Line(point1=(design.rib1, design.rib2 - (design.a/2)), point2=
		    (design.rib1, design.rib2))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
		    False, entity=model.sketches['__profile__'].geometry[9])
		model.sketches['__profile__'].Line(point1=(design.rib1, 
		    design.rib2), point2=(0.0, design.rib2))
		model.sketches['__profile__'].HorizontalConstraint(
		    addUndoState=False, entity=
		    model.sketches['__profile__'].geometry[10])
		model.sketches['__profile__'].PerpendicularConstraint(
		    addUndoState=False, entity1=
		    model.sketches['__profile__'].geometry[9], entity2=
		    model.sketches['__profile__'].geometry[10])
		model.sketches['__profile__'].Line(point1=(0.0, 
		    design.rib2), point2=(0.0, 0.0))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
		    False, entity=model.sketches['__profile__'].geometry[11])
		model.sketches['__profile__'].PerpendicularConstraint(
		    addUndoState=False, entity1=
		    model.sketches['__profile__'].geometry[10], entity2=
		    model.sketches['__profile__'].geometry[11])

	elif typeOfRib == 'outer_ribs' and typeOfRib2 == 'closed':

		#Create shape
		#1 -> 2
		model.sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
				    design.rib1, 0.0))
		model.sketches['__profile__'].HorizontalConstraint(
				    addUndoState=False, entity=
				    model.sketches['__profile__'].geometry[2])
		#2 -> 3
		model.sketches['__profile__'].Line(point1=(design.rib1, 0.0), point2=(
				    design.rib1, design.rib2))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
				    False, entity=model.sketches['__profile__'].geometry[3])
		model.sketches['__profile__'].PerpendicularConstraint(
				    addUndoState=False, entity1=
				    model.sketches['__profile__'].geometry[2], entity2=
				    model.sketches['__profile__'].geometry[3])

		#3 -> 4
		model.sketches['__profile__'].Line(point1=(design.rib1, design.rib2), point2=(
				    0.0, design.rib2))
		model.sketches['__profile__'].HorizontalConstraint(
				    addUndoState=False, entity=
				    model.sketches['__profile__'].geometry[4])

		#4 -> 1
		model.sketches['__profile__'].Line(point1=(0.0, design.rib2), point2=(
				    0.0, 0.0))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
				    False, entity=model.sketches['__profile__'].geometry[5])
		model.sketches['__profile__'].PerpendicularConstraint(
				    addUndoState=False, entity1=
				    model.sketches['__profile__'].geometry[4], entity2=
				    model.sketches['__profile__'].geometry[5])

		#5 -> 6
		model.sketches['__profile__'].Line(point1=(design.a, design.a), point2=(
				    design.rib1 - design.a, design.a))
		model.sketches['__profile__'].HorizontalConstraint(
				    addUndoState=False, entity=
				    model.sketches['__profile__'].geometry[6])

		#6 -> 7
		model.sketches['__profile__'].Line(point1=(design.rib1 - design.a, design.a), point2=(
				    design.rib1 - design.a, design.rib2 - design.a))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
				    False, entity=model.sketches['__profile__'].geometry[7])
		model.sketches['__profile__'].PerpendicularConstraint(
				    addUndoState=False, entity1=
				    model.sketches['__profile__'].geometry[6], entity2=
				    model.sketches['__profile__'].geometry[7])

		#7 -> 8
		model.sketches['__profile__'].Line(point1=(design.rib1 - design.a, design.rib2 - design.a), point2=(
				    design.a, design.rib2 - design.a))
		model.sketches['__profile__'].HorizontalConstraint(
				    addUndoState=False, entity=
				    model.sketches['__profile__'].geometry[8])

		#8 -> 5
		model.sketches['__profile__'].Line(point1=(design.a, design.rib2 - design.a), point2=(
				    design.a, design.a))
		model.sketches['__profile__'].VerticalConstraint(addUndoState=
				    False, entity=model.sketches['__profile__'].geometry[9])
		model.sketches['__profile__'].PerpendicularConstraint(
				    addUndoState=False, entity1=
				    model.sketches['__profile__'].geometry[8], entity2=
				    model.sketches['__profile__'].geometry[9])

	# PART AND INSTANCE OPERATIONS

	if typeOfRib == 'outer_ribs':

		#Shell creation
		model.Part(dimensionality=THREE_D, name='Rib-outer', type=
		    DEFORMABLE_BODY)
		model.parts['Rib-outer'].BaseShell(sketch=
		    model.sketches['__profile__'])
		del model.sketches['__profile__']

		#General set "all" for the part
		model.parts['Rib-outer'].Set(faces=
		    model.parts['Rib-outer'].faces.findAt(((design.a/2,design.a/2,0.0),),)
		    , name='rib-all') 

		#Section
		model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
		    integrationRule=SIMPSON, material='Alu_rib', name='Section-Rib', numIntPts=5, 
		    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
		    thickness=design.ribt, thicknessField='', thicknessModulus=None, thicknessType=
		    UNIFORM, useDensity=OFF)

		#Section assignment, assign material to part. This part is not a shell, therefore, no thickness is defined when specifying the section.
		model.parts['Rib-outer'].SectionAssignment(offset=0.0, offsetField=''
		    , offsetType=MIDDLE_SURFACE, region=
		    model.parts['Rib-outer'].sets['rib-all'], sectionName=
		    'Section-Rib', thicknessAssignment=FROM_SECTION)

		#### Instance operations ####
		model.rootAssembly.DatumCsysByDefault(CARTESIAN)

		#Wing root rib
		model.rootAssembly.Instance(dependent=ON, name='Rib-root-inst', part=
		    model.parts['Rib-outer']) #Create instance
		model.rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 1.0, 
		    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Rib-root-inst', ))
		model.rootAssembly.translate(instanceList=('Rib-root-inst', ), 
		    vector=(0.0, 0.0, design.C3))
		model.rootAssembly.translate(instanceList=('Rib-root-inst', ), 
		    vector=(0.0, design.cutDown, 0.0))
		model.rootAssembly.translate(instanceList=('Rib-root-inst', ), 
		    vector=(design.cutWingRoot, 0.0, 0.0))

		#Wing tip rib
		model.rootAssembly.Instance(dependent=ON, name='Rib-tip-inst', part=
		    model.parts['Rib-outer']) #Create instance
		model.rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 1.0, 
		    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Rib-tip-inst', ))
		model.rootAssembly.translate(instanceList=('Rib-tip-inst', ), 
		    vector=(0.0, 0.0, design.C3))
		model.rootAssembly.translate(instanceList=('Rib-tip-inst', ), 
		    vector=(0.0, design.cutDown, 0.0))
		model.rootAssembly.translate(instanceList=('Rib-tip-inst', ), 
		    vector=(design.cutWingTip, 0.0, 0.0))

	elif typeOfRib == 'inner_ribs':

		#Shell creation
		model.Part(dimensionality=THREE_D, name='Rib-inner', type=
		    DEFORMABLE_BODY)
		model.parts['Rib-inner'].BaseShell(sketch=
		    model.sketches['__profile__'])
		del model.sketches['__profile__']

		#General set "all" for the part
		model.parts['Rib-inner'].Set(faces=
		    model.parts['Rib-inner'].faces.findAt(((design.a/2,design.a/2,0.0),),)
		    , name='rib-all') 

		#Section
		model.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
		    integrationRule=SIMPSON, material='Alu_rib', name='Section-Rib-inner', numIntPts=5, 
		    poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
		    thickness=design.ribt_inner, thicknessField='', thicknessModulus=None, thicknessType=
		    UNIFORM, useDensity=OFF)

		#Section assignment, assign material to part. This part is not a shell, therefore, no thickness is defined when specifying the section.
		model.parts['Rib-inner'].SectionAssignment(offset=0.0, offsetField=''
		    , offsetType=MIDDLE_SURFACE, region=
		    model.parts['Rib-inner'].sets['rib-all'], sectionName=
		    'Section-Rib-inner', thicknessAssignment=FROM_SECTION)

		#Wing root rib
		totalLength = design.cutWingTip - design.cutWingRoot
		stepInX = totalLength / (design.innerRibs_n + 1)

		design.innerRibsXpos = [design.cutWingRoot + (i*stepInX) for i in range(1, int(design.innerRibs_n) + 1)]

		instances_ribs = ()

		for r in range(int(design.innerRibs_n)):

			instanceRib = model.rootAssembly.Instance(dependent=ON, name='Rib-inner-inst-'+str(r+1), part=
			    model.parts['Rib-inner']) #Create instance
			model.rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 1.0, 
			    0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('Rib-inner-inst-'+str(r+1), ))
			model.rootAssembly.translate(instanceList=('Rib-inner-inst-'+str(r+1), ), 
			    vector=(0.0, 0.0, design.C3))
			model.rootAssembly.translate(instanceList=('Rib-inner-inst-'+str(r+1), ), 
			    vector=(0.0, design.cutDown, 0.0))
			model.rootAssembly.translate(instanceList=('Rib-inner-inst-'+str(r+1), ), 
			    vector=(design.innerRibsXpos[r], 0.0, 0.0))

			instances_ribs += (instanceRib, )

		return instances_ribs



def meshing(design, mesh, partToMesh):


	#Global Mesh - fine mesh
	partToMesh.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=mesh.fineSize) #fine mesh


	#Define x,y,z positions

	#x
	xPos1 = design.cutWingTip

	xPos2 = design.cutWingRoot

	xPos3 = [
			design.cutWingTip/2, 
			design.cutWingTip/2]

	#y

	yPos1to2 = [design.cutUp,
			design.cutUp - design.a,
			(design.cutUp + abs(design.cutDown))/2,
			(design.cutUp + abs(design.cutDown))/2,
			design.cutDown + design.a,
			design.cutDown]

	yPos3 = [design.cutUp,
			design.cutDown]

	#z

	zPos1to2 = [design.C3/2,
			design.C3/2,
			design.C3 - design.a,
			design.C3,
			design.C3/2,
			design.C3/2]

	zPos3 = design.C3

	zPos4 = mesh.d #Position to locate edge in the border with the lattice

	#Search edges for course mesh 
	edges_tuple = ()

	for i in range(len(yPos1to2)): #Rib at wing tip
		
		foundEdge = partToMesh.edges.findAt(((xPos1,yPos1to2[i],zPos1to2[i]),),)

		edges_tuple += (foundEdge,)

	for i in range(len(yPos1to2)): #Rib at wing root

		foundEdge = partToMesh.edges.findAt(((xPos2,yPos1to2[i],zPos1to2[i]),),)

		edges_tuple += (foundEdge,)

	if design.innerRibs_n == 0: #If there aren't ribs in the middle

		for i in range(len(yPos3)):

			foundEdge = partToMesh.edges.findAt(((xPos3[i],yPos3[i],zPos3),),) #Edge far from the lattice

			foundEdge2 = partToMesh.edges.findAt(((xPos3[i],yPos3[i],zPos4),),) #Edge close to the lattice

			edges_tuple += (foundEdge, foundEdge2, )

	else: #If there are ribs in the middle

		index = 0

		for xPosInner in design.innerRibsXpos:

			if xPosInner == design.innerRibsXpos[0]: #For the rib close to the root

				positionInX = xPosInner - ((xPosInner - design.cutWingRoot)/2)

			else:

				positionInX = xPosInner - ((xPosInner - design.innerRibsXpos[index-1])/2)

			for i in range(len(yPos3)):

				foundEdge = partToMesh.edges.findAt(((positionInX,yPos3[i],zPos3),),) #Edge far from the lattice

				foundEdge2 = partToMesh.edges.findAt(((positionInX,yPos3[i],zPos4),),) #Edge close to the lattice
		
				edges_tuple += (foundEdge, foundEdge2, )

			index += 1

		#For the last rib
		positionInX = design.cutWingTip - ((design.cutWingTip - design.innerRibsXpos[index-1])/2)

		for i in range(len(yPos3)):

			foundEdge = partToMesh.edges.findAt(((positionInX,yPos3[i],zPos3),),) #Edge far from the lattice

			foundEdge2 = partToMesh.edges.findAt(((positionInX,yPos3[i],zPos4),),) #Edge close to the lattice
			
			edges_tuple += (foundEdge, foundEdge2, )

	#Apply course mesh on inner ribs
	if design.innerRibs_n != 0:

		edges_tuple_InnerRibs = ()

		for xPosInner in design.innerRibsXpos:

			for i in range(len(zPos1to2)):

				foundEdge = partToMesh.edges.findAt(((xPosInner,yPos1to2[i],zPos1to2[i]),),)

				edges_tuple_InnerRibs += (foundEdge,)

	#Create set from 
	partToMesh.Set(edges=edges_tuple, name='CourseMeshEdges')
	if design.innerRibs_n != 0:
		partToMesh.Set(edges=edges_tuple_InnerRibs, name='CourseMeshEdgesForInnerRibs')

	#Apply course Mesh
	for j in range(len(edges_tuple)):
		partToMesh.seedEdgeBySize(constraint = FINER, deviationFactor=0.1, edges=edges_tuple[j], size=mesh.courseSize)

	if design.innerRibs_n != 0:
		for j in range(len(edges_tuple_InnerRibs)):
			partToMesh.seedEdgeBySize(constraint = FINER, deviationFactor=0.1, edges=edges_tuple_InnerRibs[j], size=mesh.courseSize)

	if mesh.ElemType == 'quad':

		partToMesh.setElementType(elemTypes=(
		    ElemType(elemCode=S8R, elemLibrary=STANDARD), ElemType(elemCode=STRI65, 
		    elemLibrary=STANDARD)), regions=(partToMesh.faces.getByBoundingBox(0, -10, -10, design.cutWingTip+10, design.cutUp+10, design.C3 +10), ))

	# #Generate mesh
	partToMesh.generateMesh()

def mergeInstances(model, instancesToMerge, newName):

	model.rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
	    instances=instancesToMerge, name=newName, originalInstances=SUPPRESS)

def PostProc_linear(iterStr, design, load):

	#j index unused

	nameToSave = iterStr + '-'

	#Postproccesing

	#Get current folder
	cwd = os.getcwd()

	o3 = session.openOdb(name='Job_current.odb')
	session.viewports['Viewport: 1'].setValues(displayedObject=o3)
	a = mdb.models['Model-1'].rootAssembly #in case we want to see the assembly
	odb = session.odbs[cwd + '\\Job_current.odb']

	#Camera control
	session.viewports['Viewport: 1'].view.setValues(nearPlane=3648.64, 
	    farPlane=6264.32, width=1361.87, height=762.166, cameraPosition=(5052.35, 
	    1765.35, -2135.95), cameraUpVector=(-0.488045, 0.779554, 0.392565), 
	    cameraTarget=(1162.05, 116.585, 555.428), viewOffsetX=-67.4425, 
	    viewOffsetY=-38.5525)

	#################################

	#Show von misses contours
	session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
	    DEFORMED, ))
	session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
	    CONTOURS_ON_DEF, ))
	session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
	    deformationScaling=UNIFORM, uniformScaleFactor=10) #Uniform deformation

	#Print deformations plus von Misses
	session.printToFile(fileName='.\\postProc\\'+nameToSave+'def_misses.tiff', format=TIFF, canvasObjects=(
	    session.viewports['Viewport: 1'], ))

	#################################

	# Plot UR3
	#Camera control
	session.viewports['Viewport: 1'].view.setValues(nearPlane=3956.75, 
	    farPlane=6168.99, width=1518.21, height=730.322, viewOffsetX=282.193, 
	    viewOffsetY=-14.7636)
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='UR', outputPosition=NODAL, refinement=(COMPONENT, 'UR3'), )

	#Print deformations plus von Misses
	session.printToFile(fileName='.\\postProc\\'+nameToSave+'ur3.tiff', format=TIFF, canvasObjects=(
	    session.viewports['Viewport: 1'], ))

	##################################
	correction = 0.0

	#Path to analyze rotation of points
	session.Path(name='pth_ur1', type=POINT_LIST, expression=((design.cutWingRoot, design.cutUp - correction, 
	    design.C3/2), (design.cutWingTip, design.cutUp - correction, design.C3/2)))
	pth_ur1 = session.paths['pth_ur1']

	#Path to analyze displacement of points
	session.Path(name='pth_u2_z', type=POINT_LIST, expression=((design.cutWingTip, design.cutUp - correction, 
				0.0), (design.cutWingTip, design.cutUp - correction, design.C3)))
	pth_u2_z = session.paths['pth_u2_z']

	session.Path(name='pth_u2_x', type=POINT_LIST, expression=((design.cutWingRoot, design.cutUp - correction, 
				load.zPos*design.C3), (design.cutWingTip, design.cutUp - correction,load.zPos*design.C3)))
	pth_u2_x = session.paths['pth_u2_x']

	#Obtain data from path - UR1
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
	    variableLabel='UR', outputPosition=NODAL, refinement=(COMPONENT, 'UR1'))

	session.XYDataFromPath(name='XYData_ur1', path=pth_ur1, includeIntersections=False, 
	    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=20, 
	    projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

	#Obtain data from path - U2
	session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
	    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))

	session.XYDataFromPath(name='XYData_u2_z', path=pth_u2_z, includeIntersections=False, 
	    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=20, 
	    projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

	session.XYDataFromPath(name='XYData_u2_x', path=pth_u2_x, includeIntersections=False, 
	    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=20, 
	    projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

	#Reaction force at reference point
	session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('REFERENCEPOINT', ))

	session.xyDataObjects.changeKey(fromName='RF:RF2 PI: ASSEMBLY N: 1', 
	     toName='XYData_rf2')

	#Move to simulation results folder
	globalChangeDir(cwd, '-postProc-'+iterStr)

	#Write ur1 data to XY file
	x0_ur1 = session.xyDataObjects['XYData_ur1']
	session.writeXYReport(fileName=nameToSave+'ur1_xOverL.rpt', xyData=(x0_ur1, ), appendMode=OFF)

	#Write u2_z data to XY file
	x0_u2_z = session.xyDataObjects['XYData_u2_z']
	session.writeXYReport(fileName=nameToSave+'u2_zOverC3.rpt', xyData=(x0_u2_z, ), appendMode=OFF)

	#Write u2_z data to XY file
	x0_u2_x = session.xyDataObjects['XYData_u2_x']
	session.writeXYReport(fileName=nameToSave+'u2_xOverL.rpt', xyData=(x0_u2_x, ), appendMode=OFF)

	#Write rf2 data to XY file
	x0_rf2 = session.xyDataObjects['XYData_rf2']
	session.writeXYReport(fileName=nameToSave+'rf2.rpt', xyData=(x0_rf2, ), appendMode=OFF)

	#Return to original working folder
	globalChangeDir(cwd, '.')

def PostProc_nonlinear(iterStr, design, load):

	#Postproccesing

	#Get current folder
	cwd = os.getcwd()

	o3 = session.openOdb(name='Job_current.odb')
	session.viewports['Viewport: 1'].setValues(displayedObject=o3)
	a = mdb.models['Model-1'].rootAssembly #in case we want to see the assembly
	odb = session.odbs[cwd + '\\Job_current.odb']

	#Camera control
	session.viewports['Viewport: 1'].view.setValues(nearPlane=3648.64, 
	    farPlane=6264.32, width=1361.87, height=762.166, cameraPosition=(5052.35, 
	    1765.35, -2135.95), cameraUpVector=(-0.488045, 0.779554, 0.392565), 
	    cameraTarget=(1162.05, 116.585, 555.428), viewOffsetX=-67.4425, 
	    viewOffsetY=-38.5525)

	#Show von misses contours
	session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
	    DEFORMED, ))
	session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
	    CONTOURS_ON_DEF, ))

	#Path to analyze rotation of points
	posMeasureU2 = np.linspace( design.cutWingRoot, design.cutWingTip, 10, endpoint = True)
	correction = 0.0

	nameStep=o3.steps['load']

	#Move to simulation results folder
	globalChangeDir(cwd, '-postProc-'+iterStr)

	for frameID in range(len(nameStep.frames)):

		session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=frameID)

		#U2
		session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
		    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))

		for x in posMeasureU2:

			session.Path(name='pth_u2_frame'+str(frameID)+'_x'+str(round(x,2)), type=POINT_LIST, expression=((x, design.cutUp - correction, 
						design.C3), (x, design.cutUp - correction, 0.0)))
			pth_u2_x = session.paths['pth_u2_frame'+str(frameID)+'_x'+str(round(x,2))]

			session.XYDataFromPath(name='u2_fr'+str(frameID)+'_x'+str(round(x,0)), path=pth_u2_x, includeIntersections=False, 
			    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=20, 
			    projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

			#Write u2 data to XY file
			x0_u2 = session.xyDataObjects['u2_fr'+str(frameID)+'_x'+str(round(x,0))]
			session.writeXYReport(fileName='u2_frame'+str(frameID)+'_x'+str(round(x,2))+'.rpt', xyData=(x0_u2, ), appendMode=OFF)

		##### UR1
		session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
		    variableLabel='UR', outputPosition=NODAL, refinement=(COMPONENT, 'UR1'))

		#Upper
		session.Path(name='pth_ur1_frame'+str(frameID)+'_up', type=POINT_LIST, expression=((design.cutWingRoot, design.cutUp - correction, 
		    design.C3/2), (design.cutWingTip, design.cutUp - correction, design.C3/2)))
		pth_ur1_up = session.paths['pth_ur1_frame'+str(frameID)+'_up']

		#Down
		session.Path(name='pth_ur1_frame'+str(frameID)+'_dn', type=POINT_LIST, expression=((design.cutWingRoot, design.cutDown + correction, 
		    design.C3/2), (design.cutWingTip, design.cutDown + correction, design.C3/2)))
		pth_ur1_dn = session.paths['pth_ur1_frame'+str(frameID)+'_dn']

		#Upper
		session.XYDataFromPath(name='ur1_fr'+str(frameID)+'_up', path=pth_ur1_up, includeIntersections=False, 
		    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=20, 
		    projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

		#Down
		session.XYDataFromPath(name='ur1_fr'+str(frameID)+'_dn', path=pth_ur1_dn, includeIntersections=False, 
		    projectOntoMesh=False, pathStyle=UNIFORM_SPACING, numIntervals=20, 
		    projectionTolerance=0, shape=UNDEFORMED, labelType=NORM_DISTANCE)

		#Write ur1 data to XY file
		#Up
		x0_ur1_up = session.xyDataObjects['ur1_fr'+str(frameID)+'_up']
		session.writeXYReport(fileName='ur1_up_frame'+str(frameID)+'.rpt', xyData=(x0_ur1_up, ), appendMode=OFF)

		#Down
		x0_ur1_dn = session.xyDataObjects['ur1_fr'+str(frameID)+'_dn']
		session.writeXYReport(fileName='ur1_dn_frame'+str(frameID)+'.rpt', xyData=(x0_ur1_dn, ), appendMode=OFF)

	#Return to original working folder
	globalChangeDir(cwd, '.')