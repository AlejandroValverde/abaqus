import os
import sys
import platform
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

		else: #Write the parameter that it's being considered for the parametric study iteration

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

if sys.version_info.major == 2:
	execfile('setUpParametricStudy.py') #Load parametric study values
elif sys.version_info.major == 3:
	exec(open("./setUpParametricStudy.py").read())

########################################

#Write parameter study definition file
writeParametricStudyDeffile('parametricStudyDef.txt', rangesDict, parameters)

cwd = os.getcwd() #Get working directory

#Study loop
iterationID = 1
# for (keyCurrent, rangeCurrent) in rangesDict.items(): #For all the parameters defined
for keyCurrent, rangeCurrent in zip(parameters, [rangesDict[para] for para in parameters]):

	if rangesDict[keyCurrent]: #Continue if range is not empty

		for valueCurrent in rangeCurrent:

			print('\n'+'\n'+'### Abaqus parametric study initialized, iteration: ' + str(iterationID))

			#Update name
			jobNameComplete = nominalDict['jobName'] + '_' + nominalDict['typeAnalysis']

			#Clear files from last iteration
			for f in os.listdir(cwd):
			    if f.startswith(jobNameComplete) or f.startswith('abaqus.rpy'):
			        os.remove(f)

			#Update Abaqus input file
			print('-> Updating Abaqus input parameters')
			writeInputParaToFile('inputAbaqus.txt', iterationID, parameters, valueCurrent, keyCurrent, nominalDict)

			#Build geometry, create and execute job, and run post-processing
			print('-> Building and executing model (nonlinear simulation)...')
			os.system('abaqus cae noGUI=mainBuildAndExecuteWingBox.py')

			#Copy job file to specific postproc folder if the program is being run in Linux
			if platform.system() == 'Windows':
				globalCopyFile(cwd, cwd+'-postProc-'+str(iterationID), jobNameComplete+'.odb', jobNameComplete+'.odb')

			#Clear files from last computation
			for f in os.listdir(cwd):
			    if f.startswith(jobNameComplete) or f.startswith('abaqus.rpy'):
			        os.remove(f)

			######################################
			##### Run linear simulation of present iteration
			nominalDict['typeAnalysis'] = 'linear'

			#Update name
			jobNameComplete = nominalDict['jobName'] + '_' + nominalDict['typeAnalysis']

			#Clear files from last iteration
			for f in os.listdir(cwd):
			    if f.startswith(jobNameComplete) or f.startswith('abaqus.rpy'):
			        os.remove(f)

			#Update Abaqus input file
			print('-> Updating Abaqus input parameters for linear simulation')
			writeInputParaToFile('inputAbaqus.txt', iterationID, parameters, valueCurrent, keyCurrent, nominalDict)
			
			#Build geometry, create and execute job, and run post-processing
			print('-> Building and executing model (linear simulation)...')
			os.system('abaqus cae noGUI=mainBuildAndExecuteWingBox.py')

			#Copy job file to specific postproc folder if the program is being run in Linux
			if platform.system() == 'Windows':
				globalCopyFile(cwd, cwd+'-postProc-'+str(iterationID), jobNameComplete+'.odb', jobNameComplete+'.odb')

			#####Return to original configuration
			nominalDict['typeAnalysis'] = 'nonlinear'

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

#Copy and input file to corresponding Post-proc folder
globalCopyFile(cwd, cwd+'-postProc', 'parametricStudyDef.txt', 'parametricStudyDef.txt')

print('---> Abaqus parametric study finished')