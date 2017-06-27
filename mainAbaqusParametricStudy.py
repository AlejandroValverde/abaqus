import os
import pdb #pdb.set_trace()

from moduleCommon import *

#Functions
def writeInputParaToFile(fileName, iter, parameters, valueCurrent, keyCurrent, nominalDict):

	file = open(fileName, 'w')

	file.write('Iter' + '\n')

	file.write(str(iter) + '\n')

	for parameter, value in zip(parameters, [nominalDict[para] for para in parameters]):

		if parameter == keyCurrent:

			file.write(keyCurrent + '\n')

			file.write(str(valueCurrent) + '\n')

		else:

			file.write(parameter + '\n')

			if parameter == 'wingBoxLength':

				if keyCurrent == 'N': #If a range of N is being considered

					file.write(str((105*valueCurrent)-151) + '\n')

				else:

					file.write(str((105*nominalDict['N'])-151) + '\n')

			else:

				file.write(str(value) + '\n')

	file.close()

def writeParametricStudyDeffile(fileName, rangesDict, parameters):

	file = open(fileName, 'w')

	count, lineNumberForRange = 0, 4

	for parameter, rangeCurrent in zip(parameters, [rangesDict[para] for para in parameters]):

		file.write(parameter + '\n')

		parameterStart = count + 1

		parameterEnd = len(rangeCurrent) + count

		if len(rangeCurrent) == 0.0:

			file.write(str(lineNumberForRange)+','+str(0)+','+str(0)+'\n')

		else:

			file.write(str(lineNumberForRange)+','+str(parameterStart)+','+str(parameterEnd)+'\n')

		count += len(rangeCurrent)

		lineNumberForRange += 2

	file.close()

# Define parameters range

rangesDict={'N' : [],#[5, 10, 20, 30, 40, 50],
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

nominalDict={'N' : 10, #Number of unit cells in spanwise direction
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
			'typeLoad' : 'linForceInnerRibs_upper_down', #'moment', 'force', 'displacement', 'linForce', 'linForceInnerRibs_upper_down', 'linForceInnerRibs_upper', 'linForceInnerRibs_middle'
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

parameters=('N', 'M', 'r', 'B', 'L', 'Cbox_t', 'rib_t', 'rib_t_inner', 'rib_a', 'C3', 'wingBoxLength', 'eOverB', 'tChiral', 
			'typeLoad', 'displImposed', 'ForceMagnitude', 'momentMagnitude', 'forceXStart', 'forceXEnd', 'forceXn', 'forceZPos',
			'innerRibs_n', 'courseSize',	'fineSize', 'maxTimeIncrement', 'initialTimeIncrement', 'minTimeIncrement',
			'maxNumInc', 'executeJob', 'executePostProc', 'damp', 'typeAnalysis', 'typeAbaqus')

iterationIDlimit = -1 #Limit number of iterations, choose -1 for no limit

########################################

#Write parameter study definition file
writeParametricStudyDeffile('parametricStudyDef.txt', rangesDict, parameters)

cwd = os.getcwd() #Get working directory

#Study loop
iterationID = 1

for (keyCurrent, rangeCurrent) in rangesDict.items(): #For all the parameters defined

	if rangesDict[keyCurrent]: #Continue if range is not empty

		for valueCurrent in rangeCurrent:

			#Clear files from last iteration if it was unsuccessful 
			for f in os.listdir(cwd):
			    if f.startswith('Job_current'):
			        os.remove(f)

			print('Abaqus parametric study initialized, iteration: ' + str(iterationID))

			#Update Abaqus input file
			print('Updating Abaqus input parameters')
			writeInputParaToFile('inputAbaqus.txt', iterationID, parameters, valueCurrent, keyCurrent, nominalDict)

			#Build geometry, create and execute job, and run post-processing
			print('Building and executing model...')
			os.system('abaqus cae noGUI=mainBuildAndExecuteWingBox.py')

			#Check if postProc folder already exists
			globalCreateDir(cwd, '-postProc')
			
			#Create folder for simulation results
			globalCreateDir(cwd, '-postProc-'+str(iterationID))

			#Move input file to PostProc folder
			# os.chdir(cwd + '\\postProc') #Move to post-processing directory
			# for f in os.listdir():
			#     if re.search(str(iterationID) + '-abaqusInput*', f):
			#         os.remove(f) #Delete "parametricStudyDef" if it already exits from previous studies
			# os.chdir(cwd) #Return to working folder
			# os.rename('inputAbaqus.txt', '.\\postProc\\' + str(iterationID) + '-abaqusInput.txt')

			#Clear files from present iteration
			for f in os.listdir(cwd):
			    if f.startswith('Job_current'):
			        os.remove(f)

			iterationID += 1

			if iterationID == (iterationIDlimit+1): #Stop loop after at specific number of iterations, if required
				break

	if iterationID == (iterationIDlimit+1): #Stop loop after at specific number of iterations, if required
		break

#Move parametric study definition to file to PostProc folder
globalChangeDir(cwd, '-postProc')
for f in os.listdir(os.getcwd()):
    if f.startswith('parametricStudyDef'):
        os.remove(f) #Delete "parametricStudyDef" if it already exits from previous studies
globalChangeDir(cwd, '.') #Return to working folder
globalCopyFile(cwd, cwd+'-postProc', 'parametricStudyDef.txt', 'parametricStudyDef.txt')

print('-> Abaqus parametric study finished')