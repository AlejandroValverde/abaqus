import os
import sys
import platform
import math
import pdb #pdb.set_trace()
import getopt

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

				#Calculate total length
				AbstandMittelpunkte=math.sqrt( math.pow((2*nominalDict['r']),2) + math.pow((2*nominalDict['L']),2))

				if keyCurrent == 'N': #If a range of N is being considered

					file.write(str((AbstandMittelpunkte*valueCurrent)-AbstandMittelpunkte) + '\n')

				else:

					file.write(str((AbstandMittelpunkte*nominalDict['N'])-AbstandMittelpunkte) + '\n')

			else:

				file.write(str(value) + '\n')

	file.close()

def writeParametricStudyDeffile(fileName, rangesDict, parameters):

	file = open(fileName, 'w')

	count = 0

	for parameter, rangeCurrent in zip(parameters, [rangesDict[para] for para in parameters]):

		file.write(parameter + '\n')

		parameterStart = count + 1

		parameterEnd = len(rangeCurrent) + count

		if len(rangeCurrent) == 0.0:	

			file.write(str(0)+','+str(0)+'\n')

		else:

			file.write(str(parameterStart)+','+str(parameterEnd)+'\n')

		count += len(rangeCurrent)

	file.close()

def readCMDoptions(argv, CMDoptionsDict):

	short_opts = "i:"
	long_opts = ["ifile="]
	try:
		opts, args = getopt.getopt(argv,short_opts,long_opts)
	except getopt.GetoptError:
		raise ValueError('ERROR: Not correct input to script')

	# check input
	if len(opts) != len(long_opts):
		raise ValueError('ERROR: Invalid number of inputs')	

	for opt, arg in opts:

		if opt in ("-i", "--ifile"):
			# postProcFolderName = arg
			CMDoptionsDict['setUpParametricStudyFile'] = arg

	return CMDoptionsDict

########################################

#Read postProc folder name from CMD
CMDoptionsDict = {}
CMDoptionsDict = readCMDoptions(sys.argv[1:], CMDoptionsDict)

if sys.version_info.major == 2:
	execfile(CMDoptionsDict['setUpParametricStudyFile']) #Load parametric study values
elif sys.version_info.major == 3:
	exec(open("./"+CMDoptionsDict['setUpParametricStudyFile']).read())

#Write parameter study definition file
writeParametricStudyDeffile('parametricStudyDef.txt', rangesDict, parameters)

cwd = os.getcwd() #Get working directory

#Study loop
iterationID = 1
currentExecutionFlag = True
for keyCurrent, rangeCurrent in zip(parameters, [rangesDict[para] for para in parameters]):

	if rangesDict[keyCurrent]: #Continue if range is not empty

		for valueCurrent in rangeCurrent:

			print('\n'+'\n'+'### Abaqus parametric study initialized, iteration: ' + str(iterationID))

			#Update job name
			if 'nonlinear' in nominalDict['typeAnalysis']:
				jobNameComplete = nominalDict['jobName'] + '_nonlinear'
			else:
				jobNameComplete = nominalDict['jobName'] + '_linear'

			#Clear files from last iteration
			for f in os.listdir(cwd):
			    if f.startswith(jobNameComplete) or f.startswith('abaqus.rpy'):
			        os.remove(f)

			#Update Abaqus input file
			print('-> Updating Abaqus input parameters')
			writeInputParaToFile('inputAbaqus.txt', iterationID, parameters, valueCurrent, keyCurrent, nominalDict)

			#Build geometry, create and execute job, and run post-processing
			if 'double' in nominalDict['typeAnalysis']: 
				print('-> Building and executing model (nonlinear + linear simulation)...')
			elif 'nonlinear' in nominalDict['typeAnalysis']:
				print('-> Building and executing model (nonlinear simulation)...')
			else:
				print('-> Building and executing model (linear simulation)...')
			os.system('abaqus cae noGUI=mainBuildAndExecuteWingBox.py')

			#Copy job file to specific postproc folder if the program is being run in Linux
			flagCopyJob = True
			additionalLinearSimExecutedFlag = True
			while flagCopyJob:
				if platform.system() == 'Windows':
					globalCopyFile(cwd, cwd+'-postProc-'+str(iterationID), jobNameComplete+'.odb', jobNameComplete+'.odb')

				#Clear files from last computation
				for f in os.listdir(cwd):
				    if f.startswith(jobNameComplete) or f.startswith('abaqus.rpy'):
				        os.remove(f)

				#Another job will be created if a linear simulation was also run
				if 'double' in nominalDict['typeAnalysis'] and additionalLinearSimExecutedFlag:
					jobNameComplete = nominalDict['jobName'] + '_linear'
					additionalLinearSimExecutedFlag = False

				else:
					flagCopyJob = False

			#Check how far the nonlinear simulation went
			globalChangeDir(cwd, '-postProc-'+str(iterationID))
			frameIDs, frameFractions = readFrameInfoFile('frameInfo.txt')
			print('-> Simulation was executed up to a load fraction of: '+str(round(frameFractions[-1], 2)))
			globalChangeDir(cwd, '.') #Return to working folder


			iterationID += 1

			if iterationID == (iterationIDlimit+1): #Stop loop after at specific number of iterations, if required
				break

	if iterationID == (iterationIDlimit+1): #Stop loop after at specific number of iterations, if required
		print('---> Abaqus parametric study manually aborted in iteration: '+str(iterationIDlimit))
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