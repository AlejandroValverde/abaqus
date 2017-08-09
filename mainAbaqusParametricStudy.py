import os
import sys
import platform
import math
import pdb #pdb.set_trace()
import getopt

from moduleCommon import *

#Functions
######################################
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

	short_opts = "i:o:"
	long_opts = ["ifile=","convControl="]
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

		elif opt in ("-o", "--convControl"):
			CMDoptionsDict['convergenceControl'] = arg

	return CMDoptionsDict

def readFrameInfoFile(fileName):

	file = open(fileName, 'r')

	lines = file.readlines()

	frameIDs, frameFractions = [], []

	for i in range(int(len(lines)/2)): #int(len(lines)/2: Number of frames, each frame comprises two lines for name and for range

		frameID = lines[(i*2)]
		frameFraction = lines[(2*i)+1]

		frameID = frameID.replace('\n','')
		frameFraction = frameFraction.replace('\n','')

		frameID = frameID.replace('\r\n','')
		frameFraction = frameFraction.replace('\r\n','')

		frameIDs += [float(frameID)]
		frameFractions += [float(frameFraction)]

	file.close()

	return frameIDs, frameFractions

def checkConvergencyAndReturnFlag(iterationID, current_nominalDict):

	#Check convergence of current simulation
	globalChangeDir(cwd, '-postProc-'+str(iterationID))
	frameIDs, frameFractions = readFrameInfoFile('frameInfo.txt')
	lastTau = frameFractions[-1]
	print('-> Simulation was executed up to a load fraction of: '+str(round(lastTau, 2)))
	globalChangeDir(cwd, '.') #Return to working folder

	criteria = CMDoptionsDict['convergenceControl']
	#Example criteria: 30mesh_60damp95
	criteria
	if 'mesh' in criteria:
		meshTauBorder = float(criteria.split('_')[0][:2])

		if lastTau < meshTauBorder/100:

			current_nominalDict['fineSize'] = current_nominalDict['fineSize'] + 1
			current_nominalDict['courseSize'] = current_nominalDict['courseSize'] + 5

			print('-> Course mesh size increased to '+str(current_nominalDict['courseSize']))
			print('-> Fine mesh size increased to '+str(current_nominalDict['fineSize']))
			flagAnotherJob = True

		elif ((float(criteria.split('_')[1][:2])/100) < lastTau < (float(criteria.split('_')[1][6:])/100)) and 'damp' in criteria:

			if current_nominalDict['damp'] == 0.0:
				current_nominalDict['damp'] = 0.00000002 #2E-8
			else:
				current_nominalDict['damp'] = current_nominalDict['damp'] * 10

			print('-> Damping factor increased to '+str(current_nominalDict['damp']))

			flagAnotherJob = True

		else:

			flagAnotherJob = False
	return flagAnotherJob, current_nominalDict

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

for keyCurrent, rangeCurrent in zip(parameters, [rangesDict[para] for para in parameters]):

	if rangesDict[keyCurrent]: #Continue if range is not empty

		for valueCurrent in rangeCurrent:

			#Initialize iteration parameters and flags for new iteration
			current_nominalDict = nominalDict.copy() #Get the original nominal dict, IMPORTANT: copy() to not change nominal dict
			internalIteration, lastTau = 1, 0
			flagAnotherJob = True

			while flagAnotherJob:

				print('\n'+'\n'+'### Abaqus parametric study initialized, iteration: ' + str(iterationID))

				#Update job name
				if 'nonlinear' in current_nominalDict['typeAnalysis']:
					jobNameComplete = current_nominalDict['jobName'] + '_nonlinear'
				else:
					jobNameComplete = current_nominalDict['jobName'] + '_linear'

				#Clear files from last iteration
				for f in os.listdir(cwd):
				    if f.startswith(jobNameComplete) or f.startswith('abaqus.rpy'):
				        os.remove(f)

				#Update Abaqus input file
				print('-> Updating Abaqus input parameters')
				writeInputParaToFile('inputAbaqus.txt', iterationID, parameters, valueCurrent, keyCurrent, current_nominalDict)

				#Build geometry, create and execute job, and run post-processing
				if 'double' in current_nominalDict['typeAnalysis']: 
					print('-> Building and executing model (nonlinear + linear simulation)...')
				elif 'nonlinear' in current_nominalDict['typeAnalysis']:
					print('-> Building and executing model (nonlinear simulation)...')
				else:
					print('-> Building and executing model (linear simulation)...')
				os.system('abaqus cae noGUI=mainBuildAndExecuteWingBox.py')

				#Copy job file to specific postproc folder if the program is being run in Linux
				#Copy nonlinear file
				if platform.system() == 'Windows':
					globalCopyFile(cwd, cwd+'-postProc-'+str(iterationID), jobNameComplete+'.odb', jobNameComplete+'_damp'+str(current_nominalDict['damp'])+'.odb')
					globalCopyFile(cwd, cwd+'-postProc-'+str(iterationID), 'model.cae', 'model.cae')
					globalCopyFile(cwd, cwd+'-postProc-'+str(iterationID), 'model.jnl', 'model.jnl')

				#Clear files from last computation
				for f in os.listdir(cwd):
				    if f.startswith(jobNameComplete) or f.startswith('abaqus.rpy'):
				        os.remove(f)

				if 'double' in current_nominalDict['typeAnalysis']: #If linear sim was executed
					if platform.system() == 'Windows':
						globalCopyFile(cwd, cwd+'-postProc-'+str(iterationID), jobNameComplete.replace('nonlinear','linear')+'.odb', jobNameComplete.replace('nonlinear','linear')+'.odb')

					#Clear files from last computation
					for f in os.listdir(cwd):
					    if f.startswith(jobNameComplete.replace('nonlinear','linear')) or f.startswith('abaqus.rpy'):
					        os.remove(f)

				#Check how far the nonlinear simulation went
				#Last iteration progress
				last_lastTau = lastTau
				#get info from last iteration 
				
				frameIDs, frameFractions = readFrameInfoFile('frameInfo.txt')
				lastTau = frameFractions[-1]
				rowFromCurrentIter = ['actual', internalIteration, keyCurrent, valueCurrent, lastTau, current_nominalDict['fineSize'], current_nominalDict['courseSize'], current_nominalDict['damp']]
				flagAnotherJob, current_nominalDict = checkConvergencyAndReturnFlag(iterationID, current_nominalDict)

				internalIteration += 1

				if internalIteration >= 3 and current_nominalDict['damp'] == 0.0 and lastTau < (float(CMDoptionsDict['convergenceControl'].split('_')[0][:2])/100):
					#If damp was not already tried and the simulation is still inside the region where control is done through mesh size
					#This will be added to another increment in mesh size
					current_nominalDict['damp'] = 0.00000002
					print('-> Convergence was achieved up to '+str(lastTau)+', try using damping now')
				elif internalIteration >= 4 and flagAnotherJob and (lastTau < last_lastTau):
					print('-> Convergence was achieved up to '+str(lastTau)+', stopping iterations')
					flagAnotherJob = False

				#Table output
				table_status = tableOutput('Simulation summary', ['iteration', 'ID', 'parameter', 'value', 'max Q_fr/Q_to', 'f_mesh', 'c_mesh', 'damp'])
				table_status.printRow(rowFromCurrentIter)
				if flagAnotherJob: #If there is going to be another job
					table_status.printRow(['next', internalIteration, keyCurrent, valueCurrent, '-', current_nominalDict['fineSize'], current_nominalDict['courseSize'], current_nominalDict['damp']])

			#Iteration finished
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

print('\n'+'\n'+'---> Abaqus parametric study finished')