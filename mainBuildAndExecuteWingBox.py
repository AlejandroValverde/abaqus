#Code to build wing box with chiral structure in one of the webs

class structtype():
    pass

#############################################
import pdb #pdb.set_trace()
import time

#Import user-defined modules
from moduleBuildWingBoxAndPostProc import *
from moduleLoadsAndBC import *
from moduleCommon import *

######################################################
## Load design parameters for parametric study
paraRead = structtype()

#Read parameter study definition file, to get number of line for each of the parameters
inputFileName = 'inputAbaqus.txt'
paraRead = loadParameters(paraRead, inputFileName)

## Material
mat = structtype()
mat.E1 = 69000.0 #N/mm^2 - For the box, aluminum 
mat.v1 = 0.3269 #For the box, aluminum 
mat.E_chiral = 3100.0 #N/mm^2 - For the chiral structure, ABS
mat.v_chiral = 0.3 #For the chiral structure, ABS
mat.E_rib = 200000 #mat.E1 * float(paraRead.E_ribOverE1) #N/mm^2, for the rib, expressed as a fraction of the main material
mat.v_rib = 0.25
mat.E2 = mat.E1 / float(paraRead.E1OverE2_simpleModel)
mat.v2 = mat.v1

## Design parameters
design = structtype()

## Type of model
design.typeOfModel = paraRead.typeOfModel #'simpleModel' 'simpleModelWithRibs', or 'completeModel'

## Unit chiral structure
design.eOverB = float(paraRead.eOverB) #typical: 0.01
design.tChiral = float(paraRead.tChiral) #typical: 0.5mm

## Chiral lattice
design.M = int(paraRead.M) #Number of unit cells in transversal direction
design.N = int(paraRead.N) #20 - Number of unit cells in spanwise direction
design.r = float(paraRead.r)  #Node radius 
design.B = float(paraRead.B)  #Node depth
design.L = float(paraRead.L) #half length

#Cutting operations - To cut by half of the ligaments
design.cutGap_y = float(paraRead.cutGap_y) #Gap between the lattice and the skin in the y direction
design.cutGap_x = float(paraRead.cutGap_x) #Gap between the lattice and the skin in the x direction

## Box structure
#The rest of the dimensions of the box structure are calculated using the lattice final dimensions as a reference
design.C3 = float(paraRead.C3) #Length in the chordwise direction
design.Ct = float(paraRead.Cbox_t) #20 #C-box wall thickness, in millimeters

## Rib structure
design.a = float(paraRead.rib_a) #60
design.ribt = float(paraRead.rib_t) #10 #Rib thickness, in millimeters
design.ribt_inner = float(paraRead.rib_t_inner) #10 #Inner ribs thickness, in millimeters
design.innerRibs_n = int(paraRead.innerRibs_n) #Choose 0 for not inner ribs
design.innerRibs_gap = 20 #Distance from the lattice, in millimeters
design.rootRibShape = paraRead.rootRibShape #'open' or 'closed'
design.tipRibShape = paraRead.tipRibShape #'open' or 'closed'

## Meshing
mesh = structtype()
mesh.d = 50 #Offset between different meshing regions
mesh.courseSize = float(paraRead.courseSize) #Global seed element size for course mesh
mesh.fineSize = float(paraRead.fineSize) #Local seed element size for fine mesh
mesh.ElemType = '' #'quad'

## Load
load = structtype()
load.typeLoad = paraRead.typeLoad 
#Types of load:
# - 'moment': Moment applied at the rib tip, through coupling and a reference point
#		*Parameters: momentMagnitude
# - 'force1': Single force applied at the tip, where the z coordinate is 0 and the y coordinate is maximum
#		*Parameters: ForceMagnitude
# - 'force2': Single force applied at the tip, where the z coordinate is C3 and the y coordinate is minimum
#		*Parameters: ForceMagnitude
# - 'displacement_tip': Displacement imposed at the tip, at the point where the z coordinate is 0 and the y coordinate is maximum
#		*Parameters: displ
# - 'displacement_lastRib': Displacement imposed at the last rib, upper flange
#		*Parameters: displ, zPos
# - 'linForce': Distributed set of single forces on the wingbox skin
#		*Parameters: ForceMagnitude, zPos, forceXEnd, forceXStart, forceXn
# - 'linForceInnerRibs_upper': Distributed set of single forces applied at each of the inner ribs and the tip outer rib, at the upper flange
#		*Parameters: ForceMagnitude, zPos
# - 'linForceInnerRibs_middle': Distributed set of single forces applied at each of the inner ribs and the tip outer rib, at the middle flange
#		*Parameters: ForceMagnitude, zPos
# - 'linForceInnerRibs_upper_down': Distributed set of single forces applied at each of the inner ribs and the tip outer rib, at the upper and lower flange
#		*Parameters: ForceMagnitude, zPos
# - 'singleForceOnLastRib_upper': Single force applied on the tip outer rib, at the upper flange
#		*Parameters: ForceMagnitude, zPos
# - 'singleForceOnLastRib_bootom': Single force applied on the tip outer rib, at the lower flange
#		*Parameters: ForceMagnitude, zPos
#Step
load.typeAnalysis = paraRead.typeAnalysis #'linear', 'nonlinear' or 'double_linear_nonlinear'
load.typeAbaqus = paraRead.typeAbaqus #'Standard' or 'Explicit'
load.maxTimeIncrement = float(paraRead.maxTimeIncrement) # A Float specifying the maximum time increment allowed. It has to be less than the total time period (1.0)
load.minTimeIncrement = float(paraRead.minTimeIncrement) # A Float specifying the minimum time increment allowed. The default value is the smaller of the suggested initial time increment
load.initialTimeIncrement = float(paraRead.initialTimeIncrement) #A Float specifying the initial time increment. The default value is the total time period for the step.
load.maxNumInc = int(paraRead.maxNumInc) #An Int specifying the maximum number of increments in a step. The default value is 100.

#Displacement
load.displ = float(paraRead.displImposed) #mm

#Moment
load.momentMagnitude = float(paraRead.momentMagnitude) #N mm

#Force
load.ForceMagnitude = float(paraRead.ForceMagnitude)
load.force = load.ForceMagnitude / int(paraRead.forceXn)
#load.force stores the single force F_unit applied at each point i. 
# i = 1:N
# F_unit = ForceMagnitude / N
# If the load is defined as applied on the inner ribs, then the total magnitude for each F_unit is calculated
# as if N = 'design.innerRibs_n'
load.points = np.linspace(float(paraRead.forceXEnd), float(paraRead.forceXStart), int(paraRead.forceXn), endpoint = True)
load.zPos = float(paraRead.forceZPos)
load.dampFlag = (float(paraRead.damp) != 0.0)
load.damp = float(paraRead.damp)

#BCs
load.typeBC = paraRead.typeBC #'clamped', 'coupling', 'encastre'
load.additionalBC = paraRead.additionalBC #'none', 'connection', 'connection_tyre', 'connection_SYS'
load.conditionNodesInnerLattice = paraRead.conditionNodesInnerLattice #'couplingThroughRF', 'tyre', 'couplingThroughCilSYS'
load.dofContraint = paraRead.dofContraint.split(',')

if load.additionalBC != 'none' and design.cutGap_y == 0.0:
	design.cutGap_y = 5.0

## Job
jobDef = structtype()
if 'nonlinear' in paraRead.typeAnalysis:
	jobDef.jobName = paraRead.jobName + '_nonlinear'
else:
	jobDef.jobName = paraRead.jobName + '_linear'
jobDef.saveJob = False
jobDef.numCpus = 4

## Session works
session = structtype()
session.executeJob = (paraRead.executeJob == 'True')
session.executePostProc = (paraRead.executePostProc == 'True')
#############################################
# Part and instances operations

#Model variable
model = mdb.models['Model-1']

#Load parameters from the Chiral design
design = internalParameters(model, design)

#Load materials
loadMaterials(model, design, load, mat)

if design.typeOfModel == 'completeModel': #Standard design

	#Build lattice structure basic elements
	buildBasicChiral(model, design)

	#Build the whole lattice
	latticePartName, latticeInstanceName = buildLattice(model, design)

	#Cut lattice
	design.cutUp = 2 * (design.M - 1) * design.heightTriangle
	design.cutDown = 0 #-60 #Obtain from inspection
	design.cutWingRoot = design.distanceCenterPoints
	design.cutWingTip = design.distanceCenterPoints * design.N

	cutLattice(model, design)

	#Build box
	boxPartName, boxInstanceName = buildBox(model, design, mesh)

	#Marge box and lattice
	mergeInstances(model, (model.rootAssembly.instances[latticeInstanceName], model.rootAssembly.instances[boxInstanceName]), 'BoxPlusLattice')

	#Build Ribs
	# The function 'buildRib' takes two parameters:
	# - typeOfRib: 	Choose 'inner_ribs' to create ribs in between the wing tip and the root, and
	#				Choose 'outer_ribs' to create ribs located at the wing tip and at the root. These can be:
	#				- typeOfRib2: Choose 'open' to make ribs with a open section and closed to make them a close profile

	ribRootInstanceName= buildRib(model, design, 'root_rib', design.rootRibShape)
	ribTipInstanceName= buildRib(model, design, 'tip_rib', design.tipRibShape) 

	if design.innerRibs_n != 0:
		instances_ribs_inner = buildRib(model, design, 'inner_ribs', [])
	else:
		instances_ribs_inner = ()

	#Merge box and lattice to pair of ribs
	mergeInstances(model, (model.rootAssembly.instances['BoxPlusLattice-1'], 
		model.rootAssembly.instances[ribRootInstanceName], model.rootAssembly.instances[ribTipInstanceName], )+instances_ribs_inner, 'RibBoxLattice')

	partToApplyMeshBCsLoads = model.parts['RibBoxLattice']
	instanceToApplyMeshBCsLoads = model.rootAssembly.instances['RibBoxLattice-1']

	#Build tyre for nodes
	if 'tyre' in load.additionalBC or load.conditionNodesInnerLattice == 'tyre':
		instances_tyres = buildTyre(model, design, load, instanceToApplyMeshBCsLoads)
		mergeInstances(model, (instanceToApplyMeshBCsLoads, )+instances_tyres, 'RibBoxLatticeTyres')

		partToApplyMeshBCsLoads = model.parts['RibBoxLatticeTyres']
		instanceToApplyMeshBCsLoads = model.rootAssembly.instances['RibBoxLatticeTyres-1']


elif design.typeOfModel == 'simpleModel':

	boxPartName, boxInstanceName = buildBasicBox(model, design)

	partToApplyMeshBCsLoads = model.parts[boxPartName]
	instanceToApplyMeshBCsLoads = model.rootAssembly.instances[boxInstanceName]

elif design.typeOfModel == 'onlyLattice':

	#Build lattice structure basic elements
	buildBasicChiral(model, design)

	#Build the whole lattice
	latticePartName, latticeInstanceName = buildLattice(model, design)

	#Cut lattice
	cutLattice(model, design)

	#Add supports
	supportsInstanceList = buildAditionalSupportsLattice(model, design, mesh)

	#Merge
	mergeInstances(model, (model.rootAssembly.instances[latticeInstanceName], )+supportsInstanceList, 'LatticeWithSupports')

	partToApplyMeshBCsLoads = model.parts['LatticeWithSupports']
	instanceToApplyMeshBCsLoads = model.rootAssembly.instances['LatticeWithSupports-1']

elif design.typeOfModel == 'simpleModelWithRibs':

	boxPartName, boxInstanceName = buildBasicBox(model, design)

	#Build Ribs
	# The function 'buildRib' takes two parameters:
	# - typeOfRib: 	Choose 'inner_ribs' to create ribs in between the wing tip and the root, and
	#				Choose 'outer_ribs' to create ribs located at the wing tip and at the root. These can be:
	#				- typeOfRib2: Choose 'open' to make ribs with a open section and closed to make them a close profile

	ribRootInstanceName= buildRib(model, design, 'root_rib', design.rootRibShape)
	ribTipInstanceName= buildRib(model, design, 'tip_rib', design.tipRibShape) 

	if design.innerRibs_n != 0:
		instances_ribs_inner = buildRib(model, design, 'inner_ribs', [])
	else:
		instances_ribs_inner = ()

	#Merge box and lattice to pair of ribs
	mergeInstances(model, (model.rootAssembly.instances[boxInstanceName], 
		model.rootAssembly.instances[ribRootInstanceName], model.rootAssembly.instances[ribTipInstanceName], )+instances_ribs_inner, 'RibBox')

	partToApplyMeshBCsLoads = model.parts['RibBox']
	instanceToApplyMeshBCsLoads = model.rootAssembly.instances['RibBox-1']


##########################
#Meshing operations

#Build mesh
meshing(design, mesh, partToApplyMeshBCsLoads)

########################
#Load and BC's

#Boundary conditions and coupling restriction definition
if not design.typeOfModel == 'onlyLattice':
	defineBCs(model, design, instanceToApplyMeshBCsLoads, load, 'coupling') #Type of BC: 'clamped' or 'coupling'

if design.typeOfModel == 'completeModel': #Standard design
	defineBCs(model, design, instanceToApplyMeshBCsLoads, load, load.typeBC)


	defineBCs(model, design, instanceToApplyMeshBCsLoads, load, 'couplingAtLatticeNodes')

	if (design.cutGap_y != 0.0 and load.additionalBC != 'none'):# or load.conditionNodesInnerLattice == 'tyre':

		defineBCs(model, design, instanceToApplyMeshBCsLoads, load, load.additionalBC)


#Load definition
loads(model, design, mesh, load, instanceToApplyMeshBCsLoads, load.typeLoad, load.typeAnalysis, load.typeAbaqus) #Type of load: 'displ', 'force' or 'distributedForce', type of analysis: 'linear' or 'nonlinear'

################################
#Job operations

#Variables
jobCurrentName = jobDef.jobName
jobExecutionFlag = True
modelName = 'Model-1'

while jobExecutionFlag:

	#Create job
	mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
	    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
	    memory=90, memoryUnits=PERCENTAGE, model=modelName, modelPrint=OFF, 
	    multiprocessingMode=DEFAULT, name=jobCurrentName, nodalOutputPrecision=SINGLE, 
	    numCpus=jobDef.numCpus, numDomains=jobDef.numCpus, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
	    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

	#Create job variable
	jobCurrent = mdb.jobs[jobCurrentName]

	#Write job to file
	if jobDef.saveJob:

		print('Job saved...')
		jobCurrent.writeInput(consistencyChecking=OFF)

	#Submit job
	if session.executeJob:

		model.rootAssembly.regenerate()

		try:
			#The job will be submitted here. It the simulation aborts, then another
			executionFlag = True
			print('Job submitted and calculating...')

			jobCurrent.submit(consistencyChecking=OFF)
			jobCurrent.waitForCompletion()

		except Exception as e:
			executionFlag = False

		if executionFlag:
			print('Job successfully completed')
		else:
			print('Job aborted')

	#Post-processing
	if session.executePostProc:
		print('Running post-processing...')

		#Get current folder
		cwd = os.getcwd()

		#Check if postProc folder already exists
		globalCreateDir(cwd, '-postProc')
		
		#Create folder for simulation results
		globalCreateDir(cwd, '-postProc-'+paraRead.Iter)

		if load.typeAnalysis == 'linear':
			PostProc_linear(paraRead.Iter, design, load, jobCurrentName)

		elif 'nonlinear' in load.typeAnalysis:
			PostProc_nonlinear(paraRead.Iter, design, load, jobCurrentName)

		#Copy input file to postProc folder
		globalCopyFile(cwd, cwd+'-postProc', inputFileName, paraRead.Iter + '-' + inputFileName.replace('.txt', '_'+'nonlinear'+'.txt'))

		#Return to original working folder
		globalChangeDir(cwd, '.')

	if 'double' in load.typeAnalysis:
		#Now the job will be submitted for linear analysis
		load.typeAnalysis = 'linear' #The program won't enter again in this if statement

		#New model
		modelName = 'Model-linear'
		mdb.Model(name=modelName, objectToCopy=mdb.models['Model-1'])
		model = mdb.models[modelName]
		
		#Change step
		model.steps['load'].setValues(description='Step for load, standard, linear', nlgeom=OFF,initialInc=1.0, maxInc=1.0, minInc=1e-05)
		if load.dampFlag:
			#Turn off damping for the linear analysis
			model.steps['load'].setValues(adaptiveDampingRatio=None, continueDampingFactors=False, stabilizationMethod=NONE)

		jobCurrentName = jobCurrentName.replace('nonlinear', 'linear')

	else:

		jobExecutionFlag = False #Terminate execution

#Save everything
mdb.saveAs(pathName='model')