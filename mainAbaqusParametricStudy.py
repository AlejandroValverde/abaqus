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

execfile('setUpParametricStudy.py') #Load parametric study values
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
			    if f.startswith('Job_current') or f.startswith('abaqus.rpy'):
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
			    if f.startswith('Job_current') or f.startswith('abaqus.rpy'):
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

#Copy ODB file and input file to corresponding Post-proc folder
globalCopyFile(cwd, cwd+'-postProc', 'parametricStudyDef.txt', 'parametricStudyDef.txt')
globalCopyFile(cwd, cwd+'-postProc', 'Job_current.odb', 'Job_current.odb')

print('-> Abaqus parametric study finished')