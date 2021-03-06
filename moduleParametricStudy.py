try:
	import matplotlib.pyplot as plt
except Exception as e:
	print('->>>> Failed to import matplotlib')
else:
	print('Matplotlib imported successfully')
import numpy as np
import pdb #pdb.set_trace()
import math
import getopt
import os
import sys

from moduleCommon import *

class tableOutput(object):
	"""docstring for ClassName"""
	def __init__(self, titleStr, columns):
		# super(ClassName, self).__init__()
		#Initialize table to be print in CMD

		#Print title
		print('\n\n'+titleStr+'\n')

		columnWidths = []		
		for columnName in columns:

			columnWidth = len(columnName)+2
			print(str(columnName).rjust(columnWidth), end='')

			columnWidths += [columnWidth]

		self.columnWidths = columnWidths

		print('\n')

	def printRow(self, rowValues):

		i = 0
		for item in rowValues:
			if isinstance(item, str): #Is string
				print(item.rjust(max(self.columnWidths[i], len(item))), end='')

			elif isinstance(item, int):
				formatSpec = '%'+str(self.columnWidths[i])+'.0f'
				print(str(formatSpec % item).rjust(self.columnWidths[i]), end='')

			else: #Is number, float
				formatSpec = '%'+str(self.columnWidths[i])+'.3f'
				print(str(formatSpec % item).rjust(self.columnWidths[i]), end='')

			i += 1

		print('\n')
		
def readCMDoptions(argv, CMDoptionsDict):

    short_opts = "i:p:s:m:t:w:"
    long_opts = ["ifile=", "plotOptions=", "saveFigure=", "plotMean=", "showTitle=", "collectTwist="]
    try:
    	opts, args = getopt.getopt(argv,short_opts,long_opts)
    except getopt.GetoptError:
        raise ValueError('ERROR: Not correct input to script')

    for opt in opts:
    	if '-h' in opt or '-help' in opt:
    		print('-> Input options:'+'\n')
    		for long_opt in long_opts:
    			print(long_opt),

    		sys.exit()

    # check input
    if len(opts) != len(long_opts):
    	raise ValueError('ERROR: Invalid number of inputs')
    # else:
    # 	for opt in opts:
    # 		opt[0][1] in short_opts
    # 		raise ValueError('ERROR: Not found option in short options')	

    for opt, arg in opts:

        if opt in ("-i", "--ifile"):
            # postProcFolderName = arg
            CMDoptionsDict['postProcFolderName'] = arg
        elif opt in ("-p", "--plotOptions"):
            # plotOptString = arg
            CMDoptionsDict['plotOptString'] = arg
        elif opt in ("-s", "--saveFigure"):
        	#Options: 'true' or 'false'
        	if arg.lower() in ('true', 't'):
        		CMDoptionsDict['flagSaveFigure'] = True
        		CMDoptionsDict['showFigures'] = True
        	elif arg.lower() in ('false', 'f'):
        		CMDoptionsDict['flagSaveFigure'] = False
        		CMDoptionsDict['showFigures'] = True
        	elif arg.lower() in ('notshow'):
        		CMDoptionsDict['flagSaveFigure'] = False
        		CMDoptionsDict['showFigures'] = False
        	else:
        		raise ValueError('ERROR: Incorrect option chosen to save figure')
        elif opt in ("-m", "--plotMean"):
        	#Options: 'true' or 'false'
        	if arg.lower() in ('true', 't'):
        		CMDoptionsDict['plotMean'] = True
        	elif arg.lower() in ('false', 'f'):
        		CMDoptionsDict['plotMean'] = False
        	else:
        		raise ValueError('ERROR: Incorrect option chosen to plot mean')
        elif opt in ("-t", "--showTitle"):
        	#Options: 'true' or 'false'
        	if arg.lower() in ('true', 't'):
        		CMDoptionsDict['showTitle'] = True
        	elif arg.lower() in ('false', 'f'):
        		CMDoptionsDict['showTitle'] = False
        	else:
        		raise ValueError('ERROR: Incorrect option chosen to show title or not')
        elif opt in ("-w", "--collectTwist"):
        	#Options: 'true' or 'false'
        	if arg.lower() in ('true', 't'):
        		CMDoptionsDict['collectTwist'] = True
        	elif arg.lower() in ('false', 'f'):
        		CMDoptionsDict['collectTwist'] = False
        	else:
        		raise ValueError('ERROR: Incorrect option chosen to show title or not')

    return CMDoptionsDict

def importParametricStudyDeffile(fileName):

	file = open(fileName, 'r')

	lines = file.readlines()

	dictOut = {'init' : [0, 0, 0]}

	for i in range(int(len(lines)/2)): #int(len(lines)/2: Number of parameters, each parameter comprises two lines for name and for range

		nameParater = lines[(i*2)]
		lineRange = lines[(2*i)+1]

		lineRange = lineRange.replace('\n','')
		nameParater = nameParater.replace('\n','')

		lineRange = lineRange.replace('\r\n','')
		nameParater = nameParater.replace('\r\n','')

		dictOut[nameParater] = [float(lineRange.split(',')[0]), float(lineRange.split(',')[1])]

	file.close()

	return dictOut
		
class caseStudy(object):

	def __init__(self, ID):

		self.id = ID

		
	def importInputData(self, fileName):

		file = open(fileName, 'r')

		lines = file.readlines()

		for i in range(int((len(lines))/2)):

			nameParater = lines[(i*2)]
			valueParater = lines[(2*i)+1]

			valueParater = valueParater.replace('\n','')
			nameParater = nameParater.replace('\n','')

			valueParater = valueParater.replace('\r\n','')
			nameParater = nameParater.replace('\r\n','')

			setattr(self, nameParater, valueParater)

		file.close()


	def import_data_from_path(self, fileName, kindY, kindX):

		file = open(fileName, 'r')

		lines = file.readlines()

		out = np.zeros(21)
		outPos = np.zeros(21)

		lineNumber, i = 1, 0

		for line in lines:

			if 4 <= lineNumber <= 24:

				out[i] = ConvertNumber(line[30:])

				outPos[i] = ConvertNumber(line[:30])

				i += 1

			lineNumber += 1

		file.close()

		setattr(self, kindY+'_'+kindX, out)

		setattr(self, kindX, outPos)

	def obtainFrameLoadFractionInfo(self, fileName):

		file = open(fileName, 'r')

		lines = file.readlines()

		frameInfoList = []

		for i in range(int((len(lines))/2)):

			valueParater = lines[(2*i)+1]

			valueParater = valueParater.replace('\n','')

			valueParater = valueParater.replace('\r\n','')

			frameInfoList += [valueParater]

		file.close()

		return frameInfoList

	def obtainTwistData(self, fileName, simString):

		file = open(fileName, 'r')

		lines = file.readlines()

		frac = []
		twist = []

		lineNumber = 1

		for line in lines:

			if lineNumber >= 7 and line in ['\r\n', '\n', '\r']:
				break

			elif lineNumber >= 6:

				frac += [ConvertNumber(line[:30])]
				twist += [[ConvertNumber(line[30:49]),
							ConvertNumber(line[49:68]),
							ConvertNumber(line[68:87]),
							ConvertNumber(line[87:106]),
							ConvertNumber(line[106:125]),
							ConvertNumber(line[125:144]),
							ConvertNumber(line[144:163]),
							ConvertNumber(line[163:182])]]

			lineNumber += 1

		file.close()

		setattr(self, simString+'twist_frac', frac)
		setattr(self, simString+'twist_tipRib', twist)

	def obtainEnergyData(self, fileName):

		file = open(fileName, 'r')

		lines = file.readlines()

		frac = []
		extWork = []
		estab = []

		lineNumber = 1

		for line in lines:

			if lineNumber >= 5 and line in ['\r\n', '\n', '\r']:
				break

			elif lineNumber >= 4:

				frac += [ConvertNumber(line[:30])]
				extWork += [ConvertNumber(line[30:49])]
				estab += [ConvertNumber(line[49:])]

			lineNumber += 1

		file.close()

		setattr(self, 'energy_frac', frac)
		setattr(self, 'energy_extWork', extWork)
		setattr(self, 'energy_estab', estab)



class dataPerFrame(caseStudy):

	def __init__(self, arg):
		super(dataPerFrame, self).__init__(1)
		self.frameID = arg

class dataPerFramePerX(caseStudy):

	def __init__(self, arg):
		super(dataPerFramePerX, self).__init__(1)
		self.xPos = arg		
		

def applyPlottingSettingsToAxesTicks(ax, plotSettings):

	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(plotSettings['axesTicks']['size']) 

	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(plotSettings['axesTicks']['size']) 

def caseDistintion(data, studyDefDict, plotSettings, CMDoptionsDict, table):

	def figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings):

		if flagDict[keyCurrent]: #Figure initialization
			
			figure, ax = plt.subplots(1, 1)
			ax.grid(which='both', **plotSettings['grid'])
			figure.set_size_inches(10, 6, forward=True)
			axDict[keyCurrent] = ax
			figDict[keyCurrent] = figure
			flagDict[keyCurrent] = False
			axDict[keyCurrent].tick_params(axis='both', **plotSettings['axesTicks'])

		return flagDict, axDict, figDict

	keysAll = [(key) for key in studyDefDict.keys()]
	flagDict = studyDefDict.fromkeys(keysAll, True)
	axDict = studyDefDict.fromkeys(keysAll, 0)
	figDict = studyDefDict.fromkeys(keysAll, 0)
	keysUsed = []
	
	if plotSettings['typeOfPlot'] in ('UR1_tau', 'UR1_frame'):
		
		counterNperKey = studyDefDict.fromkeys(keysAll, 0)
		scatterHandles = studyDefDict.fromkeys(keysAll, [])
		keysWith_UR1_tau_plot = []
		flagDict_stiff = studyDefDict.fromkeys(keysAll, True)
		axDict_stiff = studyDefDict.fromkeys(keysAll, 0)
		figDict_stiff = studyDefDict.fromkeys(keysAll, 0)

	rangeID = 1
	for case in data:

		for (keyCurrent, rangeCurrent) in studyDefDict.items():

			if studyDefDict[keyCurrent][0] <= case.id <= studyDefDict[keyCurrent][1]:

				#Store key
				if not keyCurrent in keysUsed: keysUsed.append(keyCurrent)

				#Plotting operations

				if plotSettings['typeOfPlot'] == 'energy' and case.damp != '0.0':
					flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
					flagDict[keyCurrent] = True
					if 'Force' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + '='+str(getattr(case, keyCurrent))+', $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
					elif 'displacement' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + '='+str(getattr(case, keyCurrent))+', $displ_y$=' + str(case.displ)+'mm', **plotSettings['title'])
					plotEnergy(case, plotSettings, keyCurrent, axDict[keyCurrent])

				elif plotSettings['typeOfPlot'] == 'UR1_frame':

					flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
					if 'Force' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
					elif 'displacement' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $displ_y$=' + str(case.displ)+'mm', **plotSettings['title'])
					scatterHandles[keyCurrent] = plotUR1_frame(case, plotSettings, keyCurrent, axDict[keyCurrent], counterNperKey, scatterHandles)
					counterNperKey[keyCurrent] += 1

					if not keyCurrent in keysWith_UR1_tau_plot: keysWith_UR1_tau_plot.append(keyCurrent)

				elif plotSettings['typeOfPlot'] == 'UR1_tau':

					flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
					flagDict_stiff, axDict_stiff, figDict_stiff = figureInitialization(flagDict_stiff, axDict_stiff, figDict_stiff, keyCurrent, plotSettings)
					if 'Force' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
						axDict_stiff[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
					elif 'displacement' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $displ_y$=' + str(case.displ)+'mm', **plotSettings['title'])
						axDict_stiff[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $displ_y$=' + str(case.displ)+'mm', **plotSettings['title'])

					scatterHandles[keyCurrent] = plotUR1_tau(case, table, plotSettings, keyCurrent, axDict[keyCurrent], axDict_stiff[keyCurrent], counterNperKey, scatterHandles, CMDoptionsDict)
					counterNperKey[keyCurrent] += 1

					if not keyCurrent in keysWith_UR1_tau_plot: keysWith_UR1_tau_plot.append(keyCurrent)

				elif plotSettings['typeOfPlot'] == 'plotU2_z_LastTau':

					flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
					flagDict[keyCurrent] = True
					if 'Force' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + '='+str(getattr(case, keyCurrent))+', $Q_y$=' + str(case.ForceMagnitude)+'N'+'/last frame', **plotSettings['title'])
					elif 'displacement' in case.typeLoad and CMDoptionsDict['showTitle']:
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + '='+str(getattr(case, keyCurrent))+', $displ_y$=' + str(case.displ)+'mm'+'/last frame', **plotSettings['title'])

					plotU2_z_LastTau(case, table, plotSettings, keyCurrent, axDict[keyCurrent])

				#Saving for plots that have more than one figure per parameter used
				if CMDoptionsDict['flagSaveFigure'] and plotSettings['typeOfPlot'] in ['plotU2_z_LastTau']:
					globalCreateDir(os.getcwd(), '-figures') #Create directory if it does not already exists
					axDict[keyCurrent].set_xlim([0.0, 1.0])
					figDict[keyCurrent].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyCurrent+'_'+str(getattr(case, keyCurrent))+'.pdf'))
					figDict[keyCurrent].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyCurrent+'_'+str(getattr(case, keyCurrent))+'.png'))

				if CMDoptionsDict['flagSaveFigure'] and plotSettings['typeOfPlot'] in ['energy']:
					globalCreateDir(os.getcwd(), '-figures') #Create directory if it does not already exists
					axDict[keyCurrent].set_xlim([0.0, 1.0])
					axDict[keyCurrent].set_ylim([0.0, None])
					figDict[keyCurrent].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyCurrent+'_'+str(getattr(case, keyCurrent))+'.pdf'))
					figDict[keyCurrent].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyCurrent+'_'+str(getattr(case, keyCurrent))+'.png'))

	if plotSettings['typeOfPlot'] == 'UR1_tau':

		for key in tuple(keysWith_UR1_tau_plot):

			axDict[key].legend(handles = scatterHandles[key], **plotSettings['legend'])
			axDict_stiff[key].legend(**plotSettings['legend'])

	#Save figures
	if keysUsed and CMDoptionsDict['flagSaveFigure'] and plotSettings['typeOfPlot'] in ['UR1_tau']: #If at least one plot was crated: if keysUsed
		globalCreateDir(os.getcwd(), '-figures') #Create directory if it does not already exists
		for keyUsed in keysUsed:
			axDict[keyUsed].set_xlim([0.0, 1.0])
			axDict[keyUsed].set_ylim([None, 0.0])

			figDict[keyUsed].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyUsed+'.pdf'))
			figDict[keyUsed].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyUsed+'.png'))

def plotU2_z_LastTau(classOfData, table, plotSettings, attr, ax):

	ax.set_xlabel('$z / C_3$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	indexForMaxFrame = classOfData.framesCount.index(max(classOfData.framesCount))
	j=0
	maxXforModel = 105*float(classOfData.N) #Approximate maximum x position for the model
	maxU2, minU2, posZ_maxU2, posZ_minU2, posX_maxU2, posX_minU2 = 0,0,0,0,0,0

	for x in classOfData.dataFrames[indexForMaxFrame].xPosForU2:

		zOverC3_vector = classOfData.dataFrames[indexForMaxFrame].dataU2OverX[j].zOverC3
		u2_zOverC3_vector = classOfData.dataFrames[indexForMaxFrame].dataU2OverX[j].u2_zOverC3
		ax.plot(zOverC3_vector, u2_zOverC3_vector, label='x/L='+str(round(x/maxXforModel,2)), **plotSettings['line'])
		
		#Store info of max and min U2
		indexMaxU2 = np.argmax(u2_zOverC3_vector)
		indexMinU2 = np.argmin(u2_zOverC3_vector)
		# pdb.set_trace()
		if j == 0:
			maxU2, minU2, posZ_maxU2, posZ_minU2, posX_maxU2, posX_minU2 = max(u2_zOverC3_vector),min(u2_zOverC3_vector),zOverC3_vector[indexMaxU2],zOverC3_vector[indexMinU2],x/maxXforModel,x/maxXforModel
		else:
			if max(u2_zOverC3_vector)>maxU2:
				maxU2, posZ_maxU2, posX_maxU2 = max(u2_zOverC3_vector),zOverC3_vector[indexMaxU2],x/maxXforModel
			if min(u2_zOverC3_vector)<minU2:
				minU2, posZ_minU2, posX_minU2 = min(u2_zOverC3_vector),zOverC3_vector[indexMinU2],x/maxXforModel
			# pdb.set_trace()

		j += 1

	ax.legend(**plotSettings['legend'])

	#Output on CMD
	table.printRow([attr, getattr(classOfData, attr), maxU2, posZ_maxU2, posX_maxU2, minU2, posZ_minU2, posX_minU2])


def plotUR1_frame(classOfData, plotSettings, attr, ax, counterNperKey, scatterHandles): #NOT IN USE

	ax.set_xlabel('$frame$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	i = 0
	for frame in classOfData.framesCount:

		dataThisFrame = classOfData.dataFrames[i]

		#To obtain twist from U2 difference
		indexForMaxX = dataThisFrame.xPosForU2.index(max(dataThisFrame.xPosForU2))

		if plotSettings['meanOption']:

			meanTwist = np.mean([dataThisFrame.ur1_xOverL_up[-1], dataThisFrame.ur1_xOverL_dn[-1], dataThisFrame.twistFromU2[indexForMaxX]])
			ax.plot(frame , meanTwist * (180/math.pi), marker = 'o', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line'])

		else:
		
			#UR1 - Up
			ax.plot(frame , dataThisFrame.ur1_xOverL_up[-1] * (180/math.pi), marker = 'o', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line']) #/ max(classOfData.framesCount)
			
			#UR1 - Down
			ax.plot(frame , dataThisFrame.ur1_xOverL_dn[-1] * (180/math.pi), marker = 's', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line']) #/ max(classOfData.framesCount)
			indexForMaxX = dataThisFrame.xPosForU2.index(max(dataThisFrame.xPosForU2))
			ax.plot(frame , dataThisFrame.twistFromU2[indexForMaxX] * (180/math.pi), marker = '+', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line']) #/ max(classOfData.framesCount)
		
		i += 1

	#Linear,  frame= last frame

	#Check error if twist calculated from different parts
	a = (classOfData.linear_u2_zOverC3[-1] - classOfData.linear_u2_zOverC3[0]) / float(classOfData.C3)
	b = classOfData.linear_ur1_xOverL[-1]

	if (abs(abs(a - b) / b)*100) > 5: #If error > 5%
		raise ValueError('ERROR: More than 10 faces could not be found when applying coupling conditions at the chiral nodes')	
	meanTwist_linear = np.mean([a, b])
	ax.plot([0.0, frame] , [0.0, meanTwist_linear * (180/math.pi)], linestyle = '-.', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line'])

	if plotSettings['meanOption']:

		handle1 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='o', linestyle='', label='UR1 mean, '+attr+'='+str(getattr(classOfData, attr)))
		handle2 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='', linestyle='-.', label='UR1 mean linear, '+attr+'='+str(getattr(classOfData, attr)))
				
		scatterHandles[attr] = scatterHandles[attr] + [handle1]
		scatterHandles[attr] = scatterHandles[attr] + [handle2]

	else:

		handle1 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='o', linestyle='', label='UR1 up, '+attr+'='+str(getattr(classOfData, attr)))
		handle2 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='s', linestyle='', label='UR1 down, '+attr+'='+str(getattr(classOfData, attr)))
		handle3 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='+', linestyle='', label='Diff U2 up, '+attr+'='+str(getattr(classOfData, attr)))

		scatterHandles[attr] = scatterHandles[attr] + [handle1]
		scatterHandles[attr] = scatterHandles[attr] + [handle2]
		scatterHandles[attr] = scatterHandles[attr] + [handle3]

	return scatterHandles[attr]


def plotUR1_tau(classOfData, table, plotSettings, attr, ax, ax_stiff, counterNperKey, scatterHandles, CMDoptionsDict):

	ax.set_xlabel('$Q_{frame} / Q_{total}$', **plotSettings['axes_x'])
	ax_stiff.set_xlabel('$Q_{frame} / Q_{total}$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	storeMeans = []
	storeMeans_lin = []
	storeFractions = []
	errorStore = []
	i = 0
	for fraction in classOfData.framesFraction:

		storeFractions += [fraction]

		#Find index for current frame, classOfData.dataFrames contains unsorted data
		indexFrame = classOfData.framesCount.index(i)

		dataThisFrame = classOfData.dataFrames[indexFrame]

		#To obtain twist from U2 difference
		indexForMaxX = dataThisFrame.xPosForU2.index(max(dataThisFrame.xPosForU2))

		#Mean operations
		if plotSettings['twistData']:
			meanTwist = np.mean(classOfData.twist_tipRib[i]+[dataThisFrame.twistFromU2[indexForMaxX]])
			storeMeans += [meanTwist]
			
			twistDegrees = [l * (180/math.pi) for l in classOfData.twist_tipRib[i]+[dataThisFrame.twistFromU2[indexForMaxX]]]
			maxErrorFromMean = maxErrorFromMeanFunction(twistDegrees)
			errorStore += [maxErrorFromMean]

		else:
			meanTwist = np.mean([dataThisFrame.ur1_xOverL_up[-1], dataThisFrame.ur1_xOverL_dn[-1], dataThisFrame.twistFromU2[indexForMaxX]])
			storeMeans += [meanTwist]
			
			maxErrorFromMean = maxErrorFromMeanFunction([dataThisFrame.ur1_xOverL_up[-1]* (180/math.pi), dataThisFrame.ur1_xOverL_dn[-1]* (180/math.pi), dataThisFrame.twistFromU2[indexForMaxX]* (180/math.pi)])
			errorStore += [maxErrorFromMean]

		if plotSettings['meanOption']:

			ax.plot(fraction , meanTwist * (180/math.pi), marker = 'o', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line'])

		else:
		
			#UR1 - Up
			ax.plot(fraction , dataThisFrame.ur1_xOverL_up[-1] * (180/math.pi), marker = 'o', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line']) #/ max(classOfData.framesCount)
			
			#UR1 - Down
			ax.plot(fraction , dataThisFrame.ur1_xOverL_dn[-1] * (180/math.pi), marker = 's', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line']) #/ max(classOfData.framesCount)
			indexForMaxX = dataThisFrame.xPosForU2.index(max(dataThisFrame.xPosForU2))
			ax.plot(fraction , dataThisFrame.twistFromU2[indexForMaxX] * (180/math.pi), marker = '+', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line']) #/ max(classOfData.framesCount)
		
		i += 1

	#Linear,  frame= last frame
	if plotSettings['plotLinear']:

		if plotSettings['twistData']:
			
			meanTwist_linear = np.mean(classOfData.linear_twist_tipRib[1])
			storeMeans_lin += [meanTwist_linear]
			
			twistDegrees_lin = [l * (180/math.pi) for l in classOfData.linear_twist_tipRib[1]]
			errorLinear = maxErrorFromMeanFunction(twistDegrees_lin)

		else:
			#Check error if twist calculated from different parts
			a = (classOfData.linear_u2_zOverC3[-1] - classOfData.linear_u2_zOverC3[0]) / float(classOfData.C3)
			b = classOfData.linear_ur1_xOverL[-1]
			errorLinear = (abs(a - b) / b)*100

			if False: #If error > 5%
				raise ValueError('ERROR: The error in the calculation of the twist from different parts is more than 5%')	
			meanTwist_linear = np.mean([a, b])

		ax.plot([0.0, 1.0] , [0.0, meanTwist_linear * (180/math.pi)], linestyle = '-.', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line'])

	else:

		meanTwist_linear = 0.0
		errorLinear = 0.0

	if plotSettings['meanOption']:

		handle1 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='o', linestyle='', label=str(getattr(classOfData, attr)))#'nonlinear')
		if plotSettings['plotLinear']:
			handle2 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='', linestyle='-.', label=str(getattr(classOfData, attr)))#'linear')
				
		scatterHandles[attr] = scatterHandles[attr] + [handle1]
		if plotSettings['plotLinear']:
			scatterHandles[attr] = scatterHandles[attr] + [handle2]

	else:

		handle1 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='o', linestyle='', label='UR1 up, '+attr+'='+str(getattr(classOfData, attr)))
		handle2 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='s', linestyle='', label='UR1 down, '+attr+'='+str(getattr(classOfData, attr)))
		handle3 = plt.Line2D([],[], color=plotSettings['colors'][counterNperKey[attr]], marker='+', linestyle='', label='Diff U2 up, '+attr+'='+str(getattr(classOfData, attr)))

		scatterHandles[attr] = scatterHandles[attr] + [handle1]
		scatterHandles[attr] = scatterHandles[attr] + [handle2]
		scatterHandles[attr] = scatterHandles[attr] + [handle3]
	
	###########################
	# Stiffness
	if float(max(classOfData.framesFraction)) > 0.0:
		ax_stiff.set_ylabel('$k$ [kN/rad]', **plotSettings['axes_y'])
		# Calculate set of Stiffness
		stiffList = calculateStiffness(storeMeans, storeFractions, classOfData)

		ax_stiff.plot(storeFractions , stiffList, linestyle = plotSettings['linestyles'][counterNperKey[attr]-(4*int(counterNperKey[attr]/4))], c = plotSettings['colors'][int((counterNperKey[attr]-1)/3)], label=str(getattr(classOfData, attr)), **plotSettings['line'])
	
	###############################
	#Output on CMD
	table.printRow([attr, getattr(classOfData, attr), min(storeMeans)*(180/math.pi), max(errorStore), meanTwist_linear * (180/math.pi), errorLinear])
	
	# pdb.set_trace()
	return scatterHandles[attr]

def maxErrorFromMeanFunction(values):

	meanValues = np.mean(values)
	errorList = []	
	for v in values:

		if abs(meanValues) < 0.0000001: #Treat as zero, 1E-7
			error = 0.0
		else:
			error = abs(abs(v - meanValues) / meanValues)*100

		if False:#error > 5: #If error > 5%
			raise ValueError('ERROR: The error in the calculation of the twist from different parts is more than 5%. Error is '+str(error)+'%')

		errorList += [error]


	return max(errorList)

def plotEnergy(classOfData, plotSettings, attr, ax):

	ax.set_xlabel('$Q_{frame} / Q_{total}$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	ax.plot(classOfData.energy_frac , classOfData.energy_extWork, linestyle = '-', c = 'k', label='External work', **plotSettings['line'])
	ax.plot(classOfData.energy_frac , classOfData.energy_estab, linestyle = '-.', c = 'k', label='Static dissipation', **plotSettings['line'])

	ax.legend(**plotSettings['legend'])

def calculateStiffness(defListRad, forceListStr, classOfData):

	# Convert to force-displacement

	defVectRad = np.asarray(defListRad)
	forceList = [float(key) for key in forceListStr]
	forceVect = np.asarray(forceList)
	y_total = forceVect * float(classOfData.ForceMagnitude) / 1000 #Expressed in kN
	x_total = defVectRad

	nPointsRegre = 3
	stiffArray = np.zeros(len(y_total))
	for i in range(len(y_total)):

		if i >= nPointsRegre:
			x = x_total[i-nPointsRegre:i]
			y = y_total[i-nPointsRegre:i]

			# least-squares solution to a linear matrix equation, y = mx + c
			A = np.vstack([x, np.ones(len(x))]).T
			m, c = np.linalg.lstsq(A, y)[0]

			stiffArray[i] = m

	stiffArray[0:nPointsRegre] = stiffArray[nPointsRegre]

	return np.ndarray.tolist(stiffArray)