#Code to build wing box with chiral structure in one of the webs

class structtype():
    pass

#############################################
import pdb #pdb.set_trace()

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

## Design parameters
design = structtype()

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
design.cutUp = (178 * design.M) - (145+35) #Obtained from inspection
design.cutDown = 0 #-60 #Obtain from inspection
design.cutWingRoot = 65 + 36
design.cutWingTip = 105 * design.N - 50

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

## Meshing
mesh = structtype()
mesh.d = 50 #Offset between different meshing regions
mesh.courseSize = float(paraRead.courseSize) #Global seed element size for course mesh
mesh.fineSize = float(paraRead.fineSize) #Local seed element size for fine mesh
mesh.ElemType = '' #'quad'

## Load
load = structtype()
load.typeLoad = paraRead.typeLoad #'moment', 'force', 'displacement', 'linForce', 'linForceInnerRibs_upper', 'linForceInnerRibs_upper_down', 'linForceInnerRibs_middle', 
#Step
load.typeAnalysis = paraRead.typeAnalysis #'linear' or 'nonlinear'
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

## Job
jobDef = structtype()
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

#Load materials
loadMaterials(model, design, load)

#Build lattice structure basic elements
buildBasicChiral(model, design)

#Build the whole lattice
design = internalParameters(model, design)
buildLattice(model, design)

#Cut lattice
cutLattice(model, design)

#Build box
buildBox(model, design, mesh)

#Marge box and lattice
mergeInstances(model, (model.rootAssembly.instances['All-1'], model.rootAssembly.instances['C-box-inst']), 'BoxPlusLattice')

#Build Ribs
# The function 'buildRib' takes two parameters:
# - typeOfRib: 	Choose 'inner_ribs' to create ribs in between the wing tip and the root, and
#				Choose 'outer_ribs' to create ribs located at the wing tip and at the root. These can be:
#				- typeOfRib2: Choose 'open' to make ribs with a open section and closed to make them a close profile

buildRib(model, design, 'outer_ribs', 'open') 

if design.innerRibs_n != 0:
	instances_ribs_inner = buildRib(model, design, 'inner_ribs', [])
else:
	instances_ribs_inner = ()

#Merge box and lattice to pair of ribs
mergeInstances(model, (model.rootAssembly.instances['BoxPlusLattice-1'], 
	model.rootAssembly.instances['Rib-root-inst'], model.rootAssembly.instances['Rib-tip-inst'], )+instances_ribs_inner, 'RibBoxLattice')

##########################
#Meshing operations

#Build mesh
meshing(design, mesh, model.parts['RibBoxLattice'])

########################
#Load and BC's

#Boundary conditions and coupling restriction definition
defineBCs(model, design, model.rootAssembly.instances['RibBoxLattice-1'], 'withReferencePoint') #Type of BC: 'clamped' or 'withReferencePoint'
defineBCs(model, design, model.rootAssembly.instances['RibBoxLattice-1'], 'couplingAtLatticeNodes')

#Load definition
loads(model, design, mesh, load, model.rootAssembly.instances['RibBoxLattice-1'], load.typeLoad, load.typeAnalysis, load.typeAbaqus) #Type of load: 'displ', 'force' or 'distributedForce', type of analysis: 'linear' or 'nonlinear'

################################
#Job operations

#Create job
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job_current', nodalOutputPrecision=SINGLE, 
    numCpus=jobDef.numCpus, numDomains=jobDef.numCpus, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

#Write job to file
if jobDef.saveJob:

	print('Job saved...')
	mdb.jobs['Job_current'].writeInput(consistencyChecking=OFF)

#Submit job
if session.executeJob:

	model.rootAssembly.regenerate()
	mdb.jobs['Job_current'].submit(consistencyChecking=OFF)

	print('Job submitted and calculating...')
	mdb.jobs['Job_current'].waitForCompletion()

	if str(mdb.jobs['Job_current'].status) == 'COMPLETED':

		print('Job successfully completed')

	else:

		print('Job not successfully completed')

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
		PostProc_linear(paraRead.Iter, design, load)

	elif load.typeAnalysis == 'nonlinear':
		PostProc_nonlinear(paraRead.Iter, design, load)

	#Copy input file to postProc folder
	globalCopyFile(cwd, cwd+'-postProc', inputFileName, paraRead.Iter + '-' + inputFileName)

	#Return to original working folder
	globalChangeDir(cwd, '.')
