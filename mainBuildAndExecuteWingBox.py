#Code to build wing with embedded chiral structure wing-box

#Variables to be already defined when executing this script:
# - rib_box_low
# - rib_box_low
# - mod8_DN
# - mod8_UP
# - mod7_UP
# - mod7_DN

#Import user-defined modules
from moduleBuildWingBoxAndPostProc import *
from moduleLoadsAndBC import *
from moduleCommon import *

class structtype():
    pass

######################################################
## Load design parameters for parametric study
paraRead = structtype()

#Read parameter study definition file, to get number of line for each of the parameters
inputFileName = 'inputAbaqus.txt'
paraRead = loadParameters(paraRead, inputFileName)

## Material, UNITS CHANGES TO m
mat = structtype() 
mat.E1 = 69000000000.0 #N/m^2 - For the box, aluminum 
mat.v1 = 0.3269 #For the box, aluminum 
mat.E_chiral = 3100000000.0 #N/m^2 - For the chiral structure, ABS
mat.v_chiral = 0.3 #For the chiral structure, ABS
mat.E_rib = 200000000000.0 #mat.E1 * float(paraRead.E_ribOverE1) #N/m^2, for the rib, expressed as a fraction of the main material
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

#Cutting operations - To cut by half of the ligaments
design.cutGap_y = float(paraRead.cutGap_y) #Gap between the lattice and the skin in the y direction
design.cutGap_x = float(paraRead.cutGap_x) #Gap between the lattice and the skin in the x direction

###################

## Chiral lattice

design.M = int(paraRead.M) #Number of unit cells in transversal direction
design.N = int(paraRead.N) #20 - Number of unit cells in spanwise direction
design.B = float(paraRead.B)  #Node depth

#From wing
design.width = width
design.rib_box_low = rib_box_low
design.rib_box_up = rib_box_up
design.C2r = mod8_UP[1] - mod8_DN[1]
design.C2f = mod7_UP[1] - mod7_DN[1]
design.C2diff = (design.C2f - design.C2r) / 2
design.s = width
design.C3 = (PositionRearSpar*c_0) - (PositionFrontSpar*c_0)

#Initial values for L and r
design.L0 = float(paraRead.L) #half length
design.r0 = float(paraRead.r)  #Node radius 

dataDict = searchLandrForWing(design)
design.L = dataDict['L']
design.r = dataDict['r']

######################

## Box structure
#The rest of the dimensions of the box structure are calculated using the lattice final dimensions as a reference
design.Ct = float(paraRead.Cbox_t) #20 #C-box wall thickness, in millimeters

## Rib structure
design.a = float(paraRead.rib_a) #60
design.ribt = float(paraRead.rib_t) #Rib thickness, in millimeters
design.ribt_inner = float(paraRead.rib_t) #float(paraRead.rib_t_inner) #10 #Inner ribs thickness, in millimeters
design.innerRibs_n = int(paraRead.innerRibs_n) #Choose 0 for not inner ribs
design.innerRibs_gap = 1.5 * float(paraRead.B) #Distance from the lattice, in millimeters
design.rootRibShape = paraRead.rootRibShape #'open' or 'closed'
design.tipRibShape = paraRead.tipRibShape #'open' or 'closed'

## Meshing
mesh = structtype()
mesh.d = 2*float(paraRead.B) #Offset between different meshing regions
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

if load.additionalBC == 'none' and design.cutGap_y != 0.0:
	design.cutGap_y = 0.0

## Job
jobDef = structtype()
if 'nonlinear' in paraRead.typeAnalysis:
	jobDef.jobName = paraRead.jobName + '_nonlinear'
else:
	jobDef.jobName = paraRead.jobName + '_linear'
jobDef.saveJob = False
jobDef.numCpus = 4

## Post-processing
postProcStruct = structtype()

## Session works
session = structtype()
session.executeJob = (paraRead.executeJob == 'True')
session.executePostProc = (paraRead.executePostProc == 'True')

# model = mdb.models['Model-SparAngle-1'] Model defined in other script

#Load parameters from the Chiral design
design = internalParameters(model, design)

#Load materials
loadMaterials(model, design, load, mat)

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
# boxPartName, boxInstanceName = buildBox(model, design, mesh)
boxPartName, boxInstanceName = buildBox_forWing(model, design, mesh)

#Marge box and lattice
# mergeInstances(model, (model.rootAssembly.instances[latticeInstanceName], model.rootAssembly.instances[boxInstanceName]), 'BoxPlusLattice')

#Build Ribs
# The function 'buildRib' takes two parameters:
# - typeOfRib: 	Choose 'inner_ribs' to create ribs in between the wing tip and the root, and
#				Choose 'outer_ribs' to create ribs located at the wing tip and at the root. These can be:
#				- typeOfRib2: Choose 'open' to make ribs with a open section and closed to make them a close profile

ribRootInstanceName= buildRib_forWing(model, design, pointsSkin, 'root_rib', design.rootRibShape)
ribTipInstanceName= buildRib_forWing(model, design, pointsSkin, 'tip_rib', design.tipRibShape) 

if design.innerRibs_n != 0:
	instances_ribs_inner = buildRib_forWing(model, design, pointsSkin, 'inner_ribs', [])
else:
	instances_ribs_inner = ()

#Merge box and lattice to pair of ribs
mergeInstances(model, (model.rootAssembly.instances[boxInstanceName], 
	model.rootAssembly.instances[ribRootInstanceName], model.rootAssembly.instances[ribTipInstanceName], )+instances_ribs_inner, 'RibBox')

partToApplyMeshBCsLoads = model.parts['RibBox']
instanceToApplyMeshBCsLoads = model.rootAssembly.instances['RibBox-1']
postProcStruct.finalInstanceName = 'RibBox-1'

#Build tyre for nodes
if 'tyre' in load.additionalBC or load.conditionNodesInnerLattice == 'tyre':
	instances_tyres = buildTyre(model, design, load, instanceToApplyMeshBCsLoads)

	#Remove first lattice nodes intersecting with tyres
	model.rootAssembly.InstanceFromBooleanCut(
	    cuttingInstances=instances_tyres, instanceToBeCut=model.rootAssembly.instances[latticeInstanceName], 
	    name='LatticeWithoutTyres', originalInstances=SUPPRESS)

	#Resume tyre instances
	for instance in instances_tyres:
		model.rootAssembly.resumeFeatures((instance.name))

	mergeInstances(model, (model.rootAssembly.instances['LatticeWithoutTyres-1'], )+instances_tyres, 'Lattice')

	partToApplyMeshBCsLoads = model.parts['Lattice']
	instanceToApplyMeshBCsLoads = model.rootAssembly.instances['Lattice-1']
	postProcStruct.finalInstanceName = 'Lattice-1'