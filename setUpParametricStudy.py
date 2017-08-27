# Define parameters range
parameters=('typeOfModel', 'jobName', 'E_ribOverE1', 'E1OverE2_simpleModel', 'N', 'M', 'r', 'B', 'L', 'cutGap_x', 'cutGap_y', 
			'Cbox_t', 'rib_t', 'rib_t_inner', 'innerRibs_n', 'rootRibShape', 'tipRibShape', 'rib_a', 'C3', 'wingBoxLength', 
			'eOverB', 'tChiral', 'typeLoad', 'typeBC', 'additionalBC', 'conditionNodesInnerLattice', 'dofContraint', 'displImposed', 'ForceMagnitude', 
			'momentMagnitude', 'forceXStart', 'forceXEnd', 'forceXn', 'forceZPos',
			'courseSize', 'fineSize', 'maxTimeIncrement', 'initialTimeIncrement', 'minTimeIncrement',
			'maxNumInc', 'executeJob', 'executePostProc', 'damp', 'typeAnalysis', 'typeAbaqus')

#Define nominal values of the parameters
nominalDict={'typeOfModel' : 'completeModel',
			'jobName' : 'Job_current_parametric',
			'E_ribOverE1' : 10, #N/mm^2, for the rib, expressed as a fraction of the main material
			'E1OverE2_simpleModel' : 1,
			'N' : 8, #Number of unit cells in spanwise direction
			'M' : 3, #Number of unit cells in transversal direction
			'r' : 10.0,#Node radius 
			'B' : 20.0, #Node depth
			'L' : 50.0, #half length
			'cutGap_y' : 5.0, #Gap between the lattice and the skin (mm) in the y direction. Insert here a value if you want 'cutGap_y'!=5.0 when running script with 'additionalBC'!='none' 
			'cutGap_x' : 0.0, #Gap between the lattice and the skin (mm) in the x direction.
			'Cbox_t' : 0.8, #C-box wall thickness, $t_{C}$ (mm)
			'rib_t' : 2, #Rib thickness, $t_{rib}$ (mm)
			'rib_t_inner' : 2,
			'innerRibs_n' : 2,
			'rootRibShape' : 'closed', #'Shape of the outer rib, 'closed' or 'open'
			'tipRibShape' : 'closed',
			'rib_a' : 40, #Rib dimension frame width, $a$ (mm)
			'C3' : 300, #C-box length in the chordwise direction (mm)
			'wingBoxLength' : None, #Calculated using "N" as a parameter 
			'eOverB' : 0.01, #Chiral ligament eccentricity, $e/B$ (%)
			'tChiral' : 0.5, #Chiral lattice section thickness, $t_{chiral}$ (mm)
			'typeLoad' : 'singleForceOnLastRib_upper', #'moment', 'force1', 'force2', 'displacement', 'linForce', 'linForceInnerRibs_(upper, middle, upper_down)', 'singleForceOnLastRib_(upper, down)'
			'typeBC' : 'coupling', #'coupling', 'encastre'
			'additionalBC' : 'connection_tyre', #'none', 'connection', 'connection_tyre', 'connection_SYS'
			'conditionNodesInnerLattice' : 'tyre', #'tyre', 'couplingThroughRF', 'couplingThroughCilSYS'
			'dofContraint' : '1,2,3',
			'displImposed' : -50,
			'ForceMagnitude' : -1200, #Applied force magnitude  (N)
			'momentMagnitude' : -200000,
			'forceXStart' : 0.1, #Initial x-coordinate of distributed force (mm)
			'forceXEnd' : 1, #Final x-coordinate of distributed force (mm)
			'forceXn' : 2, #Number of points where the force is applied
			'forceZPos' : 0.5, #Z-coordinate of distributed force, nondimensionalized with "C3"
			'courseSize' : 25,
			'fineSize' : 2.5,
			'maxTimeIncrement' : 0.1, #A Float specifying the maximum time increment allowed. It has to be less than the total time period (1.0)
			'initialTimeIncrement' : 0.001, #A Float specifying the initial time increment. The default value is the total time period for the step.
			'minTimeIncrement' : 0.000000000000001,
			'maxNumInc' : 1000,#An Int specifying the maximum number of increments in a step. The default value is 100.
			'executeJob' : True,
			'executePostProc' : True,
			'damp' : 0.0,
			'typeAnalysis' : 'double_linear_nonlinear',
			'typeAbaqus' : 'Standard'}

#Define 
rangesDict={'typeOfModel' : [],
			'jobName' : [], #Name to be assigned to the job
			'E_ribOverE1' : [], 
			'E1OverE2_simpleModel' : [],
			'N' : [], #[6, 8, 10],#[5, 10, 20, 30, 40, 50],
			'M' : [], #[3, 4, 5],
			'r' : [], #[5, 10, 15],
			'B' : [], #[10.0, 20.0, 30.0],
			'L' : [], #[30.0, 50.0, 70.0],
			'cutGap_y' : [],
			'cutGap_x' : [],
			'Cbox_t' : [0.4, 0.6, 0.8],#[5, 10, 20],#5, 10, 20, 30, 40, 60, 80, 100],
			'rib_t' : [],#[3, 5, 8, 10],#2, 6, 10, 20, 40],
			'rib_t_inner' : [],
			'innerRibs_n' : [],
			'rootRibShape' : [], #'Shape of the root rib, 'closed' or 'open'
			'tipRibShape' : [], #'Shape of the tip rib, 'closed' or 'open'
			'rib_a' : [],#[50, 100, 150, 200],
			'C3' : [], #[200, 300, 400],
			'wingBoxLength' : [], 
			'eOverB' : [0.005, 0.01, 0.03, 0.05],#[0.005, 0.01, 0.03, 0.05, 0.1, 0.15], 
			'tChiral' : [0.3, 0.5, 0.7],#[0.05, 0.1, 0.5, 2.0, 2.5],
			'typeLoad' : [],
			'typeBC' : [],
			'additionalBC' : [],
			'conditionNodesInnerLattice' : [],
			'dofContraint' : [],
			'displImposed' : [],
			'ForceMagnitude' : [],
			'momentMagnitude' : [],
			'forceXStart' : [],
			'forceXEnd' : [],
			'forceXn' : [],
			'forceZPos' : [], #[0.50, 0.70, 0.80, 0.90, 0.95],
			'courseSize' : [], #[25, 50, 75],
			'fineSize' : [], #[2.0, 4.0, 6.0, 8.0],
			'maxTimeIncrement' : [],
			'initialTimeIncrement' : [],
			'minTimeIncrement' : [],
			'maxNumInc' : [],
			'executeJob' : [],
			'executePostProc' : [],
			'damp' : [],
			'typeAnalysis' : [],
			'typeAbaqus' : []}

iterationIDlimit = -1 #Limit number of iterations, choose -1 for no limit