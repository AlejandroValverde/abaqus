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

import numpy as np
import math

def defineBCs(model, design, instanceToApplyLoadAndBC, typeBC):

	# BC specification

	def createCouplingLatticeNodeByRFpoint(xNode, yNode, angles, design, instanceToApplyLoadAndBC, model, facesNotFoundCounter):

		z = design.B / 2

		nameSet = 'set_x'+str(int(xNode))+'_y'+str(int(yNode))

		nameConstraint = 'constraint_x'+str(int(xNode))+'_y'+str(int(yNode))

		rf = model.rootAssembly.ReferencePoint(point=(xNode, yNode, z))
		rfRegion = regionToolset.Region(referencePoints = (model.rootAssembly.referencePoints[rf.id], ))

		#Find faces

		faces_list=[]
		for angle in angles: #xrange is the same as range
			xFace = xNode + (design.r * math.cos(angle * (math.pi/180)) )
			yFace = yNode + (design.r * math.sin(angle * (math.pi/180)) )

			face_found = instanceToApplyLoadAndBC.faces.findAt(((xFace, yFace, z),))

			if not face_found: #if empty, face couldn't be found
				facesNotFoundCounter += 1
			
			if face_found not in faces_list: #if the face hasn't already been found
				faces_list.append(face_found)
			else:
				print('Face found is being neglected, it was already found')
		
		faces_tuple = tuple(faces_list)

		model.rootAssembly.Set(faces = faces_tuple, name=nameSet)
		# faceRegion = regionToolset.Region(faces = faces_tuple)

		#Enable coupling condition
		model.Coupling(controlPoint= rfRegion, couplingType=
		    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
		    nameConstraint, surface= model.rootAssembly.sets[nameSet], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

		# model.rootAssembly.sets[nameSet]

		return facesNotFoundCounter

	if typeBC == 'clamped':

		#Create set for fixed condition
		model.rootAssembly.Set(faces = instanceToApplyLoadAndBC.faces.findAt(((design.cutWingRoot, design.cutDown + 1, design.C3/2),),), name='fixed')

		#Assign fixed condition on border line
		model.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
		    distributionType=UNIFORM, fieldName='', localCsys=None, name='fixed', 
		    region=model.rootAssembly.sets['fixed'], u1=SET, u2=SET, 
		    u3=SET, ur1=SET, ur2=SET, ur3=SET)

	elif typeBC == 'coupling' or typeBC == 'encastre': #One of the ribs faces constrained to one point in 6 dof

		rf = model.rootAssembly.ReferencePoint(point=(design.cutWingRoot, (design.cutUp+abs(design.cutDown))/2, design.C3/2))
		model.rootAssembly.Set(name='referencePoint', referencePoints=(model.rootAssembly.referencePoints[rf.id], ))

		#Face
		if design.typeOfModel == 'simpleModel' or design.typeOfModel == 'onlyLattice':
			model.rootAssembly.Set(edges = instanceToApplyLoadAndBC.edges.findAt(((design.cutWingRoot, design.cutDown, design.C3/2),),
																				((design.cutWingRoot, design.cutUp, design.C3/2),),
																				((design.cutWingRoot, design.cutUp/2, design.C3),),
																				((design.cutWingRoot, design.cutUp/2, 0.0),),
																				), name='fixed')
		else:
			faces_list=[]

			#First face
			face_found = instanceToApplyLoadAndBC.faces.findAt(((design.cutWingRoot, design.cutDown + 1, design.C3/2),),)
			faces_list.append(face_found)

			#Search more faces
			facesNotFoundCounter = 0

			yIncrement = (abs(design.cutDown) + design.cutUp)/50

			y = design.cutDown

			for i in range(50):

				y = y + yIncrement
				face_found = instanceToApplyLoadAndBC.faces.findAt(((design.cutWingRoot, y, design.B/2),),)
				if not face_found: #if empty, face couldn't be found
					facesNotFoundCounter += 1
				
				if face_found not in faces_list: #if the face hasn't already been found
					faces_list.append(face_found)
				else:
					print('Face found is being neglected, it was already found')

			faces_tuple = tuple(faces_list)
			model.rootAssembly.Set(faces = faces_tuple, name='fixed')

		if typeBC == 'coupling':
			#Enable coupling condition
			model.Coupling(controlPoint= model.rootAssembly.sets['referencePoint'], couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    'Clamped through RF point and coupling at root', surface= model.rootAssembly.sets['fixed'], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

			#Fix reference point
			model.DisplacementBC(name='fixed', createStepName='Initial', 
			    region=model.rootAssembly.sets['referencePoint'], u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
			    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
			    localCsys=None)

		elif typeBC == 'encastre':

			model.EncastreBC(createStepName='Initial', localCsys=None, 
			    name='Encastre', region=mdb.models['Model-1'].rootAssembly.sets['fixed'])

	elif typeBC == 'couplingAtLatticeNodes':

		totalLength = design.cutWingTip - design.cutWingRoot

		totalHeight = design.cutUp + abs(design.cutDown)

		#The vector Q iterates through the first set of columns of chiral nodes
		Q_i = np.arange(design.distanceCenterPoints, totalLength + design.distanceCenterPoints, design.distanceCenterPoints)
		Q_j = np.arange(0.0, totalHeight, 2 * design.heightTriangle)

		#The vector P iterates through the second set of columns of chiral nodes
		P_i = np.arange(design.distanceCenterPoints * 3/2, totalLength + design.distanceCenterPoints, design.distanceCenterPoints)
		P_j = np.arange(design.heightTriangle, totalHeight, 2 * design.heightTriangle)

		facesNotFoundCounter = 0

		#Iterate through Q

		for q_i in Q_i:

			for q_j in Q_j:

				if (q_j == Q_j[0] and q_i == Q_i[0]) or (q_j == Q_j[-1] and q_i == Q_i[0]):
					
					print('Avoiding constraints interference')

				else:

					if q_j == Q_j[0]: #Lower node, half cut
						angles = [45, 90, 135]

					elif q_j == Q_j[-1]: #Upper node, half cut
						angles = [225, 270, 315]

					elif q_i == Q_i[0] and design.rootRibShape == 'closed' and not (q_j == Q_j[0] or q_j == Q_j[-1]): #Avoids interface with nodes that are coupled to the reference point used to clamp the root
						angles = [0, 45, 315]

					else: #Nodes in the middle
						angles = [0, 45, 90, 135, 180, 225, 270, 315]

					facesNotFoundCounter = createCouplingLatticeNodeByRFpoint(q_i, q_j, angles, design, instanceToApplyLoadAndBC, model, facesNotFoundCounter)

		#Iterate through P

		for p_i in P_i:

			for p_j in P_j:

				facesNotFoundCounter = createCouplingLatticeNodeByRFpoint(p_i, p_j, [0, 45, 90, 135, 180, 225, 270, 315], design, instanceToApplyLoadAndBC, model, facesNotFoundCounter)

		if facesNotFoundCounter > 10:

			raise ValueError('ERROR: More than 10 faces could not be found when applying coupling conditions at the chiral nodes')

	else:

		raise ValueError('Not correct option chosen for boundary condition definition') 


def loads(model, design, mesh, load, instanceToApplyLoadAndBC, typeLoad, typeAnalysis, typeAbaqus):
	"""loads(model, typeLoad)

	Input:

		-

	"""

	def giveXPosOfRefPoint(design, pos):

		totalLength = design.cutWingTip - design.cutWingRoot

		return design.cutWingRoot + (totalLength*pos)

	def setCameraBeforeDefiningSet(model):

		model.rootAssembly.regenerate()
		session.viewports['Viewport: 1'].setValues(displayedObject=model.rootAssembly)
		session.viewports['Viewport: 1'].assemblyDisplay.setValues(
		    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
		session.viewports['Viewport: 1'].view.setValues(nearPlane=3926.85, 
		    farPlane=6034.91, width=925.052, height=444.988, viewOffsetX=-310.927, 
		    viewOffsetY=337.564)

	def searchNodesForSequenceOfx(xList, yPos, zPos, tupleOfNodes, instanceToApplyLoadAndBC):

		#Search nodes on upper part
		for i in range(len(xList)):

			radius = 1

			FlagSearch = True

			while FlagSearch:

				node = instanceToApplyLoadAndBC.nodes.getByBoundingSphere((xList[i], yPos, zPos), radius)
				
				if node: #Continue if node found

					#Apply force just to the first node found
					nodeTuple = node.getMask()
					maskStr0 = nodeTuple[0].split('#')
					maskStr = '[#'+str(maskStr0[1])+'#'+str(maskStr0[2]) + ']'
					tupleOfNodes += (instanceToApplyLoadAndBC.nodes.getSequenceFromMask(mask=(maskStr,),), )

					FlagSearch = False

				else:
					print('Node to apply force not found, radius increased')
					radius += 1

		return tupleOfNodes

	# Load specification

	#Create step for load
	if typeAnalysis == 'linear' and typeAbaqus.lower() == 'standard':
		model.StaticStep(description=
		    'Step for load, standard, linear', name='load', previous=
		    'Initial', nlgeom = OFF)

	elif typeAnalysis == 'nonlinear' and typeAbaqus.lower() == 'standard':
		model.StaticStep(description=
		    'Step for load, standard, nonlinear', name='load', previous=
		    'Initial', nlgeom = ON, initialInc=0.001, maxInc=load.maxTimeIncrement, minInc=load.minTimeIncrement, maxNumInc=load.maxNumInc)

	elif typeAbaqus.lower() == 'explicit':

		model.ExplicitDynamicsStep(description='Step for load, explicit, nonlinear', name=
		    'load', previous='Initial')

	if load.dampFlag:

		model.steps['load'].setValues(continueDampingFactors=
		    False, stabilizationMagnitude=load.damp, stabilizationMethod=DAMPING_FACTOR)

	if typeLoad == 'displacement':

		#Load type 1

		#Create a set where a the displacement is imposed 
		model.rootAssembly.Set(name='pointLoad', vertices=
		    instanceToApplyLoadAndBC.vertices.findAt(((design.cutWingTip,design.cutUp,0.0),),) )


		#Define displacement condition
		model.DisplacementBC(amplitude=UNSET, createStepName='load', 
		    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
		    'displacement', region=
		    model.rootAssembly.sets['pointLoad'], u1=UNSET, u2=load.displ, 
		    u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

	elif typeLoad == 'moment':

		model.rootAssembly.Set(edges = instanceToApplyLoadAndBC.edges.findAt(((design.cutWingTip, design.cutUp - design.a, design.C3/2),), 
		    ((design.cutWingTip, (design.cutUp+abs(design.cutDown))/2, design.C3 - design.a),), 
		    ((design.cutWingTip, design.cutDown + design.a, design.C3/2),),
		    ), name='momentEdgesOnRib')

		rf2 = model.rootAssembly.ReferencePoint(point=(design.cutWingTip, (design.cutUp+abs(design.cutDown))/2, design.C3/2))
		model.rootAssembly.Set(name='referencePointMoment', referencePoints=(model.rootAssembly.referencePoints[rf2.id], ))

		model.Coupling(controlPoint=
		    model.rootAssembly.sets['referencePointMoment'], 
		    couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, 
		    name='Constraint-moment', surface=
		    model.rootAssembly.sets['momentEdgesOnRib'], u1=OFF, u2=OFF, u3=
		    OFF, ur1=ON, ur2=OFF, ur3=OFF)

		model.Moment(cm1=load.momentMagnitude, createStepName='load', 
		    distributionType=UNIFORM, field='', localCsys=None, name='moment', region=
		    model.rootAssembly.sets['referencePointMoment'])

	elif typeLoad == 'force1' or typeLoad == 'force2': #Load on upper wing tip edge

		#Load type 1

		#Create a set where a the displacement is imposed
		if typeLoad == 'force1':
			model.rootAssembly.Set(name='pointLoad', vertices=
				instanceToApplyLoadAndBC.vertices.findAt(((design.cutWingTip,design.cutUp,0.0),),) )

		elif typeLoad == 'force2':
			model.rootAssembly.Set(name='pointLoad', vertices=
				instanceToApplyLoadAndBC.vertices.findAt(((design.cutWingTip,design.cutDown,design.C3),),) )

		#Define displacement condition
		model.ConcentratedForce(cf2=load.ForceMagnitude, createStepName='load', 
		    distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
		    model.rootAssembly.sets['pointLoad'])

	elif typeLoad == 'singleForceOnLastRib_upper' or typeLoad == 'singleForceOnLastRib_lower':

		setCameraBeforeDefiningSet(model)

		#Create set
		tupleOfNodes = ()

		if typeLoad == 'singleForceOnLastRib_upper':
			tupleOfNodes = searchNodesForSequenceOfx([design.cutWingTip], design.cutUp - (design.a/2), (load.zPos*design.C3), tupleOfNodes, instanceToApplyLoadAndBC)
		elif typeLoad == 'singleForceOnLastRib_lower':
			tupleOfNodes = searchNodesForSequenceOfx([design.cutWingTip], design.cutDown + (design.a/2), (load.zPos*design.C3), tupleOfNodes, instanceToApplyLoadAndBC)

		model.rootAssembly.Set(name='setForSingleLoadOnOuterRib', nodes=tupleOfNodes)

		load.force = load.ForceMagnitude #Single force

		model.ConcentratedForce(cf2=load.force, createStepName='load', 
		    distributionType=UNIFORM, field='', localCsys=None, name='singleLoadOnOuterRib', 
		    region=model.rootAssembly.sets['setForSingleLoadOnOuterRib'])


	elif (typeLoad == 'linForceInnerRibs_upper' or typeLoad == 'linForceInnerRibs_middle' or typeLoad == 'linForceInnerRibs_upper_down') and (design.innerRibs_n != 0):

		setCameraBeforeDefiningSet(model)

		#Create set
		tupleOfNodes = ()

		#Add force on the wing tip rib
		forceXPos = design.innerRibsXpos + [design.cutWingTip] #[ design.innerRibsXpos[-1] + design.innerRibsXpos[1] - design.innerRibsXpos[0]]
		
		#Search nodes on upper part
		if typeLoad == 'linForceInnerRibs_middle':
			tupleOfNodes = searchNodesForSequenceOfx(forceXPos, (design.cutUp - design.cutDown) / 2.0, (load.zPos*design.C3), tupleOfNodes, instanceToApplyLoadAndBC)
		else: #If 'upper' or 'upper_down'
			tupleOfNodes = searchNodesForSequenceOfx(forceXPos, design.cutUp - (design.a/2), (load.zPos*design.C3), tupleOfNodes, instanceToApplyLoadAndBC)

		#Search nodes on lower part, if specified, more nodes will be added
		if typeLoad == 'linForceInnerRibs_upper_down':
			tupleOfNodes = searchNodesForSequenceOfx(forceXPos, design.cutDown + (design.a/2), (load.zPos*design.C3), tupleOfNodes, instanceToApplyLoadAndBC)

		model.rootAssembly.Set(name='setForDistributedLoad', nodes=tupleOfNodes)

		#Modify load force magnitude
		if typeLoad == 'linForceInnerRibs_upper_down':
			load.force = load.ForceMagnitude / (2*len(forceXPos))
		else:
			load.force = load.ForceMagnitude / len(forceXPos)

		model.ConcentratedForce(cf2=load.force, createStepName='load', 
		    distributionType=UNIFORM, field='', localCsys=None, name='distributedLoad', 
		    region=model.rootAssembly.sets['setForDistributedLoad'])

	elif typeLoad == 'linForce':

		setCameraBeforeDefiningSet(model)

		#Create set
		tupleOfNodes = ()
		print('Searching nodes for load')
		for i in range(len(load.points)):

			FlagSearch = True
			radius = 1

			while FlagSearch:

				node = instanceToApplyLoadAndBC.nodes.getByBoundingSphere((giveXPosOfRefPoint(design, load.points[i]), design.cutUp, (load.zPos*design.C3)), radius)

				if node:
					tupleOfNodes += (node, )

					FlagSearch = False

				else:
					print('Node to apply force not found, radius increased')
					radius += 1

			print('Final radius utilized for finding the node: '+str(radius))

		model.rootAssembly.Set(name='setForDistributedLoad', nodes=tupleOfNodes)

		model.ConcentratedForce(cf2=load.force, createStepName='load', 
		    distributionType=UNIFORM, field='', localCsys=None, name='distributedLoad', 
		    region=model.rootAssembly.sets['setForDistributedLoad'])

	else:

		raise ValueError('Not correct option chosen for load definition') 