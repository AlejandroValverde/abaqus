import matplotlib.pyplot as plt
import numpy as np
import pdb #pdb.set_trace()
import math
import getopt
import os
import sys

from moduleCommon import *

def readCMDoptions(argv, CMDoptionsDict):

    short_opts = "i:p:s:m:"
    long_opts = ["ifile=", "plotOptions=", "saveFigure=", "plotMean="]
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
				print(item.rjust(self.columnWidths[i]), end='')

			else: #Is number
				formatSpec = '%'+str(self.columnWidths[i])+'.3f'
				print(str(formatSpec % item).rjust(self.columnWidths[i]), end='')

			i += 1

		print('\n')
		

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
			ax.grid(**plotSettings['grid'])
			figure.set_size_inches(10, 6, forward=True)
			axDict[keyCurrent] = ax
			figDict[keyCurrent] = figure
			flagDict[keyCurrent] = False

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

	rangeID = 1
	for case in data:

		for (keyCurrent, rangeCurrent) in studyDefDict.items():

			if studyDefDict[keyCurrent][0] <= case.id <= studyDefDict[keyCurrent][1]:

				#Store key
				if not keyCurrent in keysUsed: keysUsed.append(keyCurrent)

				#Plotting operations
				if case.typeAnalysis == 'linear':

					if plotSettings['typeOfPlot'] == 'RF':
						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						plotRF(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'K':
						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						plotK(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'UR1':
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent], **plotSettings['title'])
						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						plotUR1(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'U2_z':
						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent], **plotSettings['title'])
						plotU2_z(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'U2_x':
						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent], **plotSettings['title'])
						plotU2_x(case, plotSettings, keyCurrent, axDict[keyCurrent])

				elif 'nonlinear' in case.typeAnalysis:

					if plotSettings['typeOfPlot'] == 'UR1_frame':

						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						if 'Force' in case.typeLoad:
							axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
						elif 'displacement' in case.typeLoad:
							axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $displ_y$=' + str(case.displ)+'mm', **plotSettings['title'])
						scatterHandles[keyCurrent] = plotUR1_frame(case, plotSettings, keyCurrent, axDict[keyCurrent], counterNperKey, scatterHandles)
						counterNperKey[keyCurrent] += 1

						if not keyCurrent in keysWith_UR1_tau_plot: keysWith_UR1_tau_plot.append(keyCurrent)

					elif plotSettings['typeOfPlot'] == 'UR1_tau':

						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						if 'Force' in case.typeLoad:
							axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
						elif 'displacement' in case.typeLoad:
							axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $displ_y$=' + str(case.displ)+'mm', **plotSettings['title'])

						scatterHandles[keyCurrent] = plotUR1_tau(case, table, plotSettings, keyCurrent, axDict[keyCurrent], counterNperKey, scatterHandles)
						counterNperKey[keyCurrent] += 1

						if not keyCurrent in keysWith_UR1_tau_plot: keysWith_UR1_tau_plot.append(keyCurrent)

					elif plotSettings['typeOfPlot'] == 'plotU2_z_LastTau':

						flagDict, axDict, figDict = figureInitialization(flagDict, axDict, figDict, keyCurrent, plotSettings)
						flagDict[keyCurrent] = True
						if 'Force' in case.typeLoad:
							axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + '='+str(getattr(case, keyCurrent))+', $Q_y$=' + str(case.ForceMagnitude)+'N'+'/last frame', **plotSettings['title'])
						elif 'displacement' in case.typeLoad:
							axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + '='+str(getattr(case, keyCurrent))+', $displ_y$=' + str(case.displ)+'mm'+'/last frame', **plotSettings['title'])

						plotU2_z_LastTau(case, table, plotSettings, keyCurrent, axDict[keyCurrent])

				#Saving for plots that have more than one figure per parameter used
				if CMDoptionsDict['flagSaveFigure'] and plotSettings['typeOfPlot'] in ['plotU2_z_LastTau']:
					globalCreateDir(os.getcwd(), '-figures') #Create directory if it does not already exists
					# pdb.set_trace()
					figDict[keyCurrent].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyCurrent+'_'+str(getattr(case, keyCurrent))+'.pdf'))


	if plotSettings['typeOfPlot'] == 'UR1_tau':

		for key in tuple(keysWith_UR1_tau_plot):

			axDict[key].legend(handles = scatterHandles[key], **plotSettings['legend'])

	#Save figures
	if keysUsed and CMDoptionsDict['flagSaveFigure'] and plotSettings['typeOfPlot'] in ['UR1_tau']: #If at least one plot was crated: if keysUsed
		globalCreateDir(os.getcwd(), '-figures') #Create directory if it does not already exists
		for keyUsed in keysUsed:

			figDict[keyUsed].savefig(os.path.join('figures', plotSettings['typeOfPlot'] + '-' + keyUsed+'.pdf'))


def plotRF(classOfData, plotSettings, attr, ax):

	#Axes labels
	ax.set_xlabel(plotSettings['xLabel'][attr], **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	ax.plot([getattr(classOfData, attr)], [classOfData.rf], '-ok', **plotSettings['line'])

	applyPlottingSettingsToAxesTicks(ax, plotSettings)
	

def plotK(classOfData, plotSettings, attr, ax):

	#Axes labels
	ax.set_xlabel(plotSettings['xLabel'][attr], **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])


	ax.plot([getattr(classOfData, attr)], classOfData.K, '-ok', **plotSettings['line'])

	applyPlottingSettingsToAxesTicks(ax, plotSettings)


def plotUR1(classOfData, plotSettings, attr, ax):


	#Axes labels
	ax.set_xlabel('Distance along the wing box $x/L$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	# pdb.set_trace()
	ax.plot(classOfData.xOverL, classOfData.ur1_xOverL * (180/math.pi), '-o', label = attr+':'+str(getattr(classOfData, attr)), **plotSettings['line'])
	
	applyPlottingSettingsToAxesTicks(ax, plotSettings)

	ax.legend(**plotSettings['legend'])

def plotU2_z(classOfData, plotSettings, attr, ax):

	#Axes labels
	ax.set_xlabel('Distance along the wing box chordwise direction $z/C3$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	ax.plot(classOfData.zOverC3, classOfData.u2_zOverC3, '-o', label = attr+':'+str(getattr(classOfData, attr)), **plotSettings['line'])
	
	applyPlottingSettingsToAxesTicks(ax, plotSettings)

	ax.legend(**plotSettings['legend'])

def plotU2_x(classOfData, plotSettings, attr, ax):

	#Axes labels
	ax.set_xlabel('Distance along the wing box $x/L$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	ax.plot(classOfData.xOverL, classOfData.u2_xOverL, '-o', label = attr+':'+str(getattr(classOfData, attr)), **plotSettings['line'])
	
	applyPlottingSettingsToAxesTicks(ax, plotSettings)

	ax.legend(**plotSettings['legend'])

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
	table.printRow([attr, getattr(classOfData, attr), float(max(classOfData.framesFraction)), classOfData.fineSize, classOfData.courseSize, classOfData.damp, maxU2, posZ_maxU2, posX_maxU2, minU2, posZ_minU2, posX_minU2])


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

	if ((abs(a - b) / b)*100) > 5: #If error > 5%
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


def plotUR1_tau(classOfData, table, plotSettings, attr, ax, counterNperKey, scatterHandles):

	ax.set_xlabel('$Q_{frame} / Q_{total}$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	storeMeans = []
	errorStore = []
	i = 0
	for fraction in classOfData.framesFraction:

		#Find index for current frame, classOfData.dataFrames contains unsorted data
		indexFrame = classOfData.framesCount.index(i)

		dataThisFrame = classOfData.dataFrames[indexFrame]

		#To obtain twist from U2 difference
		indexForMaxX = dataThisFrame.xPosForU2.index(max(dataThisFrame.xPosForU2))

		#Mean operations
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

	#Check error if twist calculated from different parts
	a = (classOfData.linear_u2_zOverC3[-1] - classOfData.linear_u2_zOverC3[0]) / float(classOfData.C3)
	b = classOfData.linear_ur1_xOverL[-1]
	errorLinear = (abs(a - b) / b)*100

	if errorLinear > 5: #If error > 5%
		raise ValueError('ERROR: The error in the calculation of the twist from different parts is more than 5%')	
	meanTwist_linear = np.mean([a, b])
	ax.plot([0.0, 1.0] , [0.0, meanTwist_linear * (180/math.pi)], linestyle = '-.', c = plotSettings['colors'][counterNperKey[attr]], **plotSettings['line'])

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
	

	#Output on CMD
	table.printRow([attr, getattr(classOfData, attr), float(max(classOfData.framesFraction)), classOfData.fineSize, classOfData.courseSize, classOfData.damp, min(storeMeans)*(180/math.pi), max(errorStore), meanTwist_linear * (180/math.pi), errorLinear])
	
	# pdb.set_trace()
	return scatterHandles[attr]

def maxErrorFromMeanFunction(values):

	meanValues = np.mean(values)
	errorList = []	
	for v in values:

		if meanValues < 0.0001: #Treat as zero
			error = 0.0
		else:
			error = (abs(v - meanValues) / meanValues)*100

		if error > 5: #If error > 5%
			raise ValueError('ERROR: The error in the calculation of the twist from different parts is more than 5%. Error is '+str(error)+'%')

		errorList += [error]

	return max(errorList)

