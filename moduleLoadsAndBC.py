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

#User-defined modules
from moduleCommon import *

def defineBCs(model, design, instanceToApplyLoadAndBC, load, typeBC):

	# BC specification

	def createCouplingLatticeNode(xNode, yNode, angles, design, instanceToApplyLoadAndBC, model, facesNotFoundCounter):

		z = design.B / 2

		nameSet = 'set_x'+str(int(xNode))+'_y'+str(int(yNode))

		nameConstraint = 'Local_sys_constraint_x'+str(int(xNode))+'_y'+str(int(yNode))

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

		#Enable coupling condition
		if load.conditionNodesInnerLattice == 'couplingThroughRF':
			model.Coupling(controlPoint= rfRegion, couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    nameConstraint, surface= model.rootAssembly.sets[nameSet], u1=OFF, u2=OFF, u3=ON, ur1=ON, ur2=ON, ur3=ON)#u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

		elif load.conditionNodesInnerLattice == 'couplingThroughCilSYS':

			#Local SYS creation
			localsys = model.rootAssembly.DatumCsysByThreePoints(coordSysType=
			    CYLINDRICAL, name='Local_sys'+'_x'+str(int(xNode))+'_y'+str(int(yNode)), origin=model.rootAssembly.referencePoints[rf.id], point1=
			    (xNode + (design.r * math.cos(0 * (math.pi/180)) ), yNode + (design.r * math.sin(0 * (math.pi/180)) ), 
			    z), point2=(xNode + (design.r * math.cos(90 * (math.pi/180)) ), yNode + (design.r * math.sin(90 * (math.pi/180)) ), 
			    z))

			model.Coupling(controlPoint=rfRegion, couplingType=KINEMATIC, influenceRadius=
			    WHOLE_SURFACE, localCsys=model.rootAssembly.datums[localsys.id], 
			    name=nameConstraint, surface= model.rootAssembly.sets[nameSet], u1=ON, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)



		return facesNotFoundCounter

	def searchNodesForSequenceOfZ(xRange, yRange, zList, tupleOfNodes, instanceToApplyLoadAndBC):

		#Search nodes on upper part
		for y in yRange:
			for x in xRange:
				for i in range(len(zList)):

					radius = 1

					FlagSearch = True

					while FlagSearch:

						node = instanceToApplyLoadAndBC.nodes.getByBoundingSphere((x, y, zList[i]), radius)
						
						if node: #Continue if node found

							#Get a single node from the MeshNodeArray object found (node)
							tupleOfNodes += (instanceToApplyLoadAndBC.nodes.sequenceFromLabels(labels=(node[0].label,),), )

							FlagSearch = False

						elif radius < 50:
							print('Node to apply force not found, radius increased')
							radius += 1

						else:
							raise ValueError('ERROR: Radius bigger than 50')

					print('Final radius utilized for finding the node: '+str(radius))


		return tupleOfNodes, radius

	def createCouplingLatticeWithSkin(xNode, xRange, yNode, ySkin, design, instanceToApplyLoadAndBC, model, load):

		z = design.B / 2

		nameSet = 'set_skin_x'+str(int(xNode))+'_y'+str(int(yNode))

		nameConstraint = 'constraint_extra_x'+str(int(xNode))+'_y'+str(int(yNode))

		rf = model.rootAssembly.ReferencePoint(point=(xNode, yNode, z))
		rfRegion = regionToolset.Region(referencePoints = (model.rootAssembly.referencePoints[rf.id], ))

		#Search nodes on skin
		#Create set
		tupleOfNodes = ()
		zList = np.linspace(0.0, design.B, 10)

		tupleOfNodes, finalRadius = searchNodesForSequenceOfZ(xRange, ySkin, zList, tupleOfNodes, instanceToApplyLoadAndBC)

		model.rootAssembly.Set(name=nameSet, nodes=tupleOfNodes)

		#Enable coupling condition
		if load.additionalBC == 'couplingNodesUp':
			model.Coupling(controlPoint= rfRegion, couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    nameConstraint, surface= model.rootAssembly.sets[nameSet], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF) #influenceRadius?? = WHOLE_SURFACE or float?

		elif load.additionalBC == 'couplingNodesUp_x1_free':
			model.Coupling(controlPoint= rfRegion, couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    nameConstraint, surface= model.rootAssembly.sets[nameSet], u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF) #influenceRadius?? = WHOLE_SURFACE or float?

	def createDoubleCouplingLatticeWithSkin(xNode, xRange, yNode, ySkin, design, angles, instanceToApplyLoadAndBC, model, load):

		z = design.B / 2

		nameSet = 'set_double_x'+str(int(xNode))+'_y'+str(int(yNode))

		nameSetSkin = 'set_double_skin_x'+str(int(xNode))+'_y'+str(int(yNode))

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
		    nameConstraint, surface= model.rootAssembly.sets[nameSet], u1=OFF, u2=OFF, u3=ON, ur1=ON, ur2=ON, ur3=ON)

		#Search nodes on skin
		#Create set
		tupleOfNodes = ()
		zList = np.linspace(0.0, design.B, 10)

		tupleOfNodes, finalRadius = searchNodesForSequenceOfZ(xRange, ySkin, zList, tupleOfNodes, instanceToApplyLoadAndBC)

		model.rootAssembly.Set(name=nameSetSkin, nodes=tupleOfNodes)

		#Enable coupling condition
		if load.additionalBC == 'couplingNodesUp':
			model.Coupling(controlPoint= rfRegion, couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    'skin_'+nameConstraint, surface= model.rootAssembly.sets[nameSetSkin], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF) #influenceRadius?? = WHOLE_SURFACE or float?

		elif load.additionalBC == 'couplingNodesUp_x1_free':
			model.Coupling(controlPoint= rfRegion, couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    'skin_'+nameConstraint, surface= model.rootAssembly.sets[nameSetSkin], u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF) #influenceRadius?? = WHOLE_SURFACE or float?


	def createCouplingLatticeWithSkin2(xNode, xRange, yNode, ySkin, design, instanceToApplyLoadAndBC, model, load):

		z = design.B / 2

		nameSet = 'set_skin_x'+str(int(xNode))+'_y'+str(int(yNode))
		nameConstraint = 'constraint_extra_x'+str(int(xNode))+'_y'+str(int(yNode))
		nameSetLattice = 'set_lattice_nodes_x'+str(int(xNode))+'_y'+str(int(yNode))

		#Reference -> nodes on lattice
		#Find nodes in faces

		#Search nodes
		zList = np.linspace(0.0, design.B, 10)

		tupleOfNodes = ()
		angles = [0, 45, 90, 135, 180, 225, 270, 315]
		for angle in angles: #xrange is the same as range
			xFace = xNode + (design.r * math.cos(angle * (math.pi/180)) )
			yFace = yNode + (design.r * math.sin(angle * (math.pi/180)) )


			tupleOfNodes, finalRadius = searchNodesForSequenceOfZ([xFace], yFace, zList, tupleOfNodes, instanceToApplyLoadAndBC)

		model.rootAssembly.Set(nodes = tupleOfNodes, name=nameSetLattice)

		######################

		#Search nodes on skin
		#Create set
		tupleOfNodes = ()
		zList = np.linspace(0.0, design.B, 10)

		tupleOfNodes, finalRadius = searchNodesForSequenceOfZ(xRange, ySkin, zList, tupleOfNodes, instanceToApplyLoadAndBC)

		model.rootAssembly.Set(name=nameSet, nodes=tupleOfNodes)

		#Enable coupling condition
		if load.additionalBC == 'couplingNodesUp':
			model.Coupling(controlPoint= model.rootAssembly.sets[nameSet], couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    nameConstraint, surface= model.rootAssembly.sets[nameSetLattice], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF) #influenceRadius?? = WHOLE_SURFACE or float?

		elif load.additionalBC == 'couplingNodesUp_x1_free':
			model.Coupling(controlPoint= model.rootAssembly.sets[nameSet], couplingType=
			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
			    nameConstraint, surface= model.rootAssembly.sets[nameSetLattice], u1=ON, u2=OFF, u3=ON, ur1=ON, ur2=ON, ur3=OFF) #influenceRadius?? = WHOLE_SURFACE or float?
	
	def createCouplingWithTyre(xNode, yNode, xSkin, ySkin, Range, design, instanceToApplyLoadAndBC, model, load):

		z = design.B / 2

		if xNode == xSkin:
			typeConnection = 'upper_and_lower'
		elif yNode == ySkin:
			typeConnection = 'middle'

		nameSetSkin = 'set_aboveTyre_skin_x'+str(int(xNode))+'_y'+str(int(yNode))
		nameNodeTyre = 'node_tyre_x'+str(int(xNode))+'_y'+str(int(yNode))

		#Set skin
		tupleOfNodes = ()
		zList = np.linspace(0.0, design.B, 10)
		if typeConnection == 'upper_and_lower' and load.rangeOfNodesOnSkinForCouplingFlag:
			tupleOfNodes, finalRadius = searchNodesForSequenceOfZ(Range, [ySkin], zList, tupleOfNodes, instanceToApplyLoadAndBC)
		elif typeConnection == 'middle' and load.rangeOfNodesOnSkinForCouplingFlag:
			tupleOfNodes, finalRadius = searchNodesForSequenceOfZ([xSkin], Range, zList, tupleOfNodes, instanceToApplyLoadAndBC)
		else:
			tupleOfNodes, finalRadius = searchNodesForSequenceOfZ([xSkin], [ySkin], [z], tupleOfNodes, instanceToApplyLoadAndBC)
		model.rootAssembly.Set(name=nameSetSkin, nodes=tupleOfNodes)

		#Node on tyre
		tupleOfNodes = ()
		tupleOfNodes, finalRadius = searchNodesForSequenceOfZ([xNode], [yNode], [z], tupleOfNodes, instanceToApplyLoadAndBC)
		model.rootAssembly.Set(name=nameNodeTyre, nodes=tupleOfNodes)


		#Equation condition
		# if load.additionalBC == 'couplingNodesUp':
		# 	dof_vect = [1, 2, 3, 4, 5]
		# elif load.additionalBC == 'couplingNodesUp_x1_free':
		# 	dof_vect = [2, 3, 4, 5]

		# for dof in dof_vect:
		# 	model.Equation(name='eq_dof_'+str(dof)+'_x'+str(int(xNode))+'_y'+str(int(yNode)), terms=((1.0, 
		# 	    nameSetSkin, dof), (-1.0, nameNodeTyre, dof)))

		nameCoupling = 'constraint_tyre_x'+str(int(xNode))+'_y'+str(int(yNode))
		if load.additionalBC == 'couplingNodesUp':
			if typeConnection == 'upper_and_lower':
				model.Coupling(controlPoint= model.rootAssembly.sets[nameSetSkin], couplingType=
							    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
							    nameCoupling, surface= model.rootAssembly.sets[nameNodeTyre], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF)
			elif typeConnection == 'middle':
				model.Coupling(controlPoint= model.rootAssembly.sets[nameSetSkin], couplingType=
							    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
							    nameCoupling, surface= model.rootAssembly.sets[nameNodeTyre], u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF)

		elif load.additionalBC == 'couplingNodesUp_x1_free':
			if typeConnection == 'upper_and_lower':
				model.Coupling(controlPoint= model.rootAssembly.sets[nameSetSkin], couplingType=
							    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
							    nameCoupling, surface= model.rootAssembly.sets[nameNodeTyre], u1=OFF, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=OFF)
			elif typeConnection == 'middle':
				model.Coupling(controlPoint= model.rootAssembly.sets[nameSetSkin], couplingType=
							    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
							    nameCoupling, surface= model.rootAssembly.sets[nameNodeTyre], u1=ON, u2=OFF, u3=ON, ur1=ON, ur2=ON, ur3=OFF)
		# nameCoupling = '2_constraint_extra_x'+str(int(xNode))+'_y'+str(int(yNode))
		# model.Coupling(controlPoint= model.rootAssembly.sets[nameNodeTyre], couplingType=
		# 			    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
		# 			    nameCoupling, surface= model.rootAssembly.sets[nameSetSkin], u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)

	def createCouplingWithLocalSYS(xNode, yNode, xSkin, ySkin, Range, design, instanceToApplyLoadAndBC, model, load):

		z = design.B / 2

		nameSet = 'faces_x'+str(int(xNode))+'_y'+str(int(yNode))
		nameSetSkin = 'set_skin_x'+str(int(xNode))+'_y'+str(int(yNode))
		nameConstraint = 'Local_sys_constraint_x'+str(int(xNode))+'_y'+str(int(yNode))
		nameCoupling = 'couplingWithSkin_x'+str(int(xNode))+'_y'+str(int(yNode))
		nameRFset = 'referencePoint_x'+str(int(xNode))+'_y'+str(int(yNode))

		rf = model.rootAssembly.ReferencePoint(point=(xNode, yNode, z))
		rfRegion = regionToolset.Region(referencePoints = (model.rootAssembly.referencePoints[rf.id], ))
		model.rootAssembly.Set(name=nameRFset, referencePoints=(model.rootAssembly.referencePoints[rf.id], ))

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

		#Internal coupling through local SYS 
		localsys = model.rootAssembly.DatumCsysByThreePoints(coordSysType=
		    CYLINDRICAL, name='Local_sys'+'_x'+str(int(xNode))+'_y'+str(int(yNode)), origin=model.rootAssembly.referencePoints[rf.id], point1=
		    (xNode + (design.r * math.cos(0 * (math.pi/180)) ), yNode + (design.r * math.sin(0 * (math.pi/180)) ), 
		    z), point2=(xNode + (design.r * math.cos(90 * (math.pi/180)) ), yNode + (design.r * math.sin(90 * (math.pi/180)) ), 
		    z))

		model.Coupling(controlPoint=rfRegion, couplingType=KINEMATIC, influenceRadius=
		    WHOLE_SURFACE, localCsys=model.rootAssembly.datums[localsys.id], 
		    name=nameConstraint, surface= model.rootAssembly.sets[nameSet], u1=ON, u2=OFF, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)

		#Skin
		tupleOfNodes = ()
		tupleOfNodes, finalRadius = searchNodesForSequenceOfZ([xSkin], [ySkin], [z], tupleOfNodes, instanceToApplyLoadAndBC)
		model.rootAssembly.Set(name=nameSetSkin, nodes=tupleOfNodes)

		#Coupling
		model.Coupling(controlPoint= model.rootAssembly.sets[nameSetSkin], couplingType=
		    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
		    nameCoupling, surface= model.rootAssembly.sets[nameRFset], u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)

	########################################
	#Main code for BC's
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

			yIncrement = (abs(design.cutDown_effective) + design.cutUp_effective)/50

			y = design.cutDown

			for i in range(48):

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

	elif typeBC == 'couplingAtLatticeNodes' or typeBC != 'none':

		Q_i, Q_j, P_i, P_j = getQandPvectors(design)

		if typeBC == 'couplingAtLatticeNodes':

			facesNotFoundCounter = 0

			#Iterate through Q

			for q_i in Q_i:

				for q_j in Q_j:

					if (q_j in [Q_j[0], Q_j[-1]] and q_i in [Q_i[0], Q_i[-1]]) and design.cutGap == 0.0:
						
						print('Avoiding constraints interference')

					else:

						if q_j == Q_j[0]: #Lower node, half cut

							if load.additionalBC != 'none':
								angles = []
							elif design.cutGap == 0.0:
								angles = [45, 90, 135]
							else:
								angles = [0, 45, 90, 135, 180, 225, 270, 315]

						elif q_j == Q_j[-1]: #Upper node, half cut

							if load.additionalBC != 'none':
								angles = []
							elif design.cutGap == 0.0:
								angles = [225, 270, 315]
							else:
								angles = [0, 45, 90, 135, 180, 225, 270, 315]

						elif q_i == Q_i[0] and design.rootRibShape == 'closed' and not (q_j == Q_j[0] or q_j == Q_j[-1]): #Avoids interface with nodes that are coupled to the reference point used to clamp the root
							angles = [0, 45, 315]

						else: #Nodes in the middle
							if (q_i == Q_i[0] or q_i == Q_i[-1]) and design.cutGap == 0.0 and load.additionalBC != 'none': #Node at the extremes with no gap
								angles = []
							else:
								angles = [0, 45, 90, 135, 180, 225, 270, 315]

						if angles: #If vector angles is not empty
							facesNotFoundCounter = createCouplingLatticeNode(q_i, q_j, angles, design, instanceToApplyLoadAndBC, model, facesNotFoundCounter)

			#Iterate through P

			for p_i in P_i:

				for p_j in P_j:

					facesNotFoundCounter = createCouplingLatticeNode(p_i, p_j, [0, 45, 90, 135, 180, 225, 270, 315], design, instanceToApplyLoadAndBC, model, facesNotFoundCounter)

			if facesNotFoundCounter > 10:

				raise ValueError('ERROR: More than 10 faces could not be found when applying coupling conditions at the chiral nodes')

		else: #For coupling at the gap between lattice and skin

			load.rangeOfNodesOnSkinForCouplingFlag = False

			#Iterate through Q

			for q_i in Q_i:

				for q_j in Q_j:

					rangeX = np.linspace(q_i - design.r, q_i + design.r, 10)
					rangeY = np.linspace(q_j - design.r, q_j + design.r, 10)

					if (not (q_j == Q_j[0] or q_j == Q_j[-1])) and (q_i == Q_i[0] or q_i == Q_i[-1]):

						if q_i == Q_i[0]: #Root
							createCouplingWithLocalSYS(q_i, q_j, q_i - (design.r + design.cutGap), q_j, rangeY, design, instanceToApplyLoadAndBC, model, load)

						elif q_i == Q_i[-1]: #Tip
							createCouplingWithLocalSYS(q_i, q_j, q_i + (design.r + design.cutGap), q_j, rangeY, design, instanceToApplyLoadAndBC, model, load)

					else:

						angles = [0, 45, 90, 135, 180, 225, 270, 315]

						if q_j == Q_j[0]: #Lower node, half cut

							createCouplingWithLocalSYS(q_i, q_j, q_i, q_j - (design.r + design.cutGap), rangeX, design, instanceToApplyLoadAndBC, model, load)

						elif q_j == Q_j[-1]: #Upper node, half cut
							
							createCouplingWithLocalSYS(q_i, q_j, q_i, q_j + design.r + design.cutGap, rangeX, design, instanceToApplyLoadAndBC, model, load)

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

					#Get a single node from the MeshNodeArray object found (node)
					tupleOfNodes += (instanceToApplyLoadAndBC.nodes.sequenceFromLabels(labels=(node[0].label,),), )

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

	if typeLoad == 'displacement_tip':

		#Create a set where a the displacement is imposed 
		model.rootAssembly.Set(name='pointLoad', vertices=
		    instanceToApplyLoadAndBC.vertices.findAt(((design.cutWingTip,design.cutUp,0.0),),) )


		#Define displacement condition
		model.DisplacementBC(amplitude=UNSET, createStepName='load', 
		    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
		    'displacement', region=
		    model.rootAssembly.sets['pointLoad'], u1=UNSET, u2=load.displ, 
		    u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

	elif typeLoad == 'displacement_lastRib':

		#Create set
		tupleOfNodes = ()
		tupleOfNodes = searchNodesForSequenceOfx([design.cutWingTip], design.cutUp - (design.a/2), (load.zPos*design.C3), tupleOfNodes, instanceToApplyLoadAndBC)

		#Create a set where a the displacement is imposed 
		model.rootAssembly.Set(name='setForDisplacementOnOuterRib', nodes=tupleOfNodes)

		#Define displacement condition
		model.DisplacementBC(amplitude=UNSET, createStepName='load', 
		    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
		    'displacement', region=
		    model.rootAssembly.sets['setForDisplacementOnOuterRib'], u1=UNSET, u2=load.displ, 
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