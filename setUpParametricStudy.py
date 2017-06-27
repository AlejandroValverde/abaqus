# Define parameters range

rangesDict={'E1overE_rib' : [], 
			'N' : [],#[5, 10, 20, 30, 40, 50],
			'M' : [],
			'r' : [],
			'B' : [],
			'L' : [],
			'Cbox_t' : [8],#[5, 10, 20],#5, 10, 20, 30, 40, 60, 80, 100],
			'rib_t' : [],#[3, 5, 8, 10],#2, 6, 10, 20, 40],
			'rib_t_inner' : [],
			'rib_a' : [],#[50, 100, 150, 200],
			'C3' : [],
			'wingBoxLength' : [], 
			'eOverB' : [],#[0.005, 0.01, 0.03, 0.05, 0.1, 0.15], 
			'tChiral' : [],#[0.05, 0.1, 0.5, 2.0, 2.5],
			'typeLoad' : [],
			'displImposed' : [],
			'ForceMagnitude' : [],
			'momentMagnitude' : [],
			'forceXStart' : [],
			'forceXEnd' : [],
			'forceXn' : [],
			'forceZPos' : [],#[0.50, 0.70, 0.80, 0.90, 0.95],
			'innerRibs_n' : [],#3, 5, 7, 9, 11]}
			'courseSize' : [],
			'fineSize' : [],
			'maxTimeIncrement' : [],
			'initialTimeIncrement' : [],
			'minTimeIncrement' : [],
			'maxNumInc' : [],
			'executeJob' : [],
			'executePostProc' : [],
			'damp' : [],
			'typeAnalysis' : [],
			'typeAbaqus' : []}

nominalDict={'E1overE_rib' : 1, #N/mm^2, for the rib, expressed as a fraction of the main material
			'N' : 10, #Number of unit cells in spanwise direction
			'M' : 3, #Number of unit cells in transversal direction
			'r' : 10.0,#Node radius 
			'B' : 30.0, #Node depth
			'L' : 50.0, #half length
			'Cbox_t' : 8, #C-box wall thickness, $t_{C}$ (mm)
			'rib_t' : 3, #Rib thickness, $t_{rib}$ (mm)
			'rib_t_inner' : 3,
			'rib_a' : 30, #Rib dimension frame width, $a$ (mm)
			'C3' : 400, #C-box length in the chordwise direction (mm)
			'wingBoxLength' : None, #Calculated using "N" as a parameter 
			'eOverB' : 0.01, #Chiral ligament eccentricity, $e/B$ (%)
			'tChiral' : 1.0, #Chiral lattice section thickness, $t_{chiral}$ (mm)
			'typeLoad' : 'linForceInnerRibs_upper', #'moment', 'force', 'displacement', 'linForce', 'linForceInnerRibs_(upper, middle, upper_down)', 'singleForceOnLastRib_(upper, down)'
			'displImposed' : 50,
			'ForceMagnitude' : -8000, #Applied force magnitude  (N)
			'momentMagnitude' : -200000,
			'forceXStart' : 0.1, #Initial x-coordinate of distributed force (mm)
			'forceXEnd' : 1, #Final x-coordinate of distributed force (mm)
			'forceXn' : 2, #Number of points where the force is applied
			'forceZPos' : 0.5, #Z-coordinate of distributed force, nondimensionalized with "C3"
			'innerRibs_n' : 2,
			'courseSize' : 30,
			'fineSize' : 4.0,
			'maxTimeIncrement' : 0.1, #A Float specifying the maximum time increment allowed. It has to be less than the total time period (1.0)
			'initialTimeIncrement' : 0.001, #A Float specifying the initial time increment. The default value is the total time period for the step.
			'minTimeIncrement' : 0.000000000000001,
			'maxNumInc' : 1000,#An Int specifying the maximum number of increments in a step. The default value is 100.
			'executeJob' : True,
			'executePostProc' : True,
			'damp' : 0.7,
			'typeAnalysis' : 'nonlinear',
			'typeAbaqus' : 'Standard'}

parameters=('E1overE_rib', 'N', 'M', 'r', 'B', 'L', 'Cbox_t', 'rib_t', 'rib_t_inner', 'rib_a', 'C3', 'wingBoxLength', 'eOverB', 'tChiral', 
			'typeLoad', 'displImposed', 'ForceMagnitude', 'momentMagnitude', 'forceXStart', 'forceXEnd', 'forceXn', 'forceZPos',
			'innerRibs_n', 'courseSize',	'fineSize', 'maxTimeIncrement', 'initialTimeIncrement', 'minTimeIncrement',
			'maxNumInc', 'executeJob', 'executePostProc', 'damp', 'typeAnalysis', 'typeAbaqus')

iterationIDlimit = -1 #Limit number of iterations, choose -1 for no limit