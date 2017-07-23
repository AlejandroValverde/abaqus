import matplotlib.pyplot as plt
import numpy as np
import pdb #pdb.set_trace()
import math
import getopt

from moduleCommon import *

def readCMDoptions(argv):

    try:
    	opts, args = getopt.getopt(argv,"i:o:",["ifile=", "plotOptions="])

    except getopt.GetoptError:
        raise ValueError('ERROR: Not correct input to script')

    for opt, arg in opts:

        if opt in ("-i", "--ifile"):
            postProcFolderName = arg
        elif opt in ("-o", "--plotOptions"):
            plotOptString = arg

    return postProcFolderName, plotOptString

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

		dictOut[nameParater] = [float(lineRange.split(',')[0]), float(lineRange.split(',')[1]), float(lineRange.split(',')[2])]

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

def caseDistintion(data, studyDefDict, plotSettings):

	def figureInitialization(flagDict, axDict, keyCurrent, plotSettings):

		if flagDict[keyCurrent]: #Figure initialization

			figure, ax = plt.subplots(1, 1)
			ax.grid(**plotSettings['grid'])
			figure.set_size_inches(10, 6, forward=True)
			axDict[keyCurrent] = ax
			flagDict[keyCurrent] = False

		return flagDict, axDict

	keysAll = [(key) for key in studyDefDict.keys()]
	flagDict = studyDefDict.fromkeys(keysAll, True)
	axDict = studyDefDict.fromkeys(keysAll, 0)
	
	if plotSettings['typeOfPlot'] == 'UR1_tau' or plotSettings['typeOfPlot'] == 'UR1_frame':
		counterNperKey = studyDefDict.fromkeys(keysAll, 0)
		scatterHandles = studyDefDict.fromkeys(keysAll, [])
		keysWith_UR1_tau_plot = []
	
	rangeID = 1
	for case in data:

		for (keyCurrent, rangeCurrent) in studyDefDict.items():

			if studyDefDict[keyCurrent][0] <= case.id <= studyDefDict[keyCurrent][1]:

				#Plotting operations
				if case.typeAnalysis == 'linear':

					if plotSettings['typeOfPlot'] == 'RF':
						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)
						plotRF(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'K':
						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)
						plotK(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'UR1':
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent], **plotSettings['title'])
						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)
						plotUR1(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'U2_z':
						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent], **plotSettings['title'])
						plotU2_z(case, plotSettings, keyCurrent, axDict[keyCurrent])

					elif plotSettings['typeOfPlot'] == 'U2_x':
						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent], **plotSettings['title'])
						plotU2_x(case, plotSettings, keyCurrent, axDict[keyCurrent])

				elif case.typeAnalysis == 'nonlinear':

					if plotSettings['typeOfPlot'] == 'UR1_frame':

						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)

						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
						scatterHandles[keyCurrent] = plotUR1_frame(case, plotSettings, keyCurrent, axDict[keyCurrent], counterNperKey, scatterHandles)
						counterNperKey[keyCurrent] += 1

						if not keyCurrent in keysWith_UR1_tau_plot: keysWith_UR1_tau_plot.append(keyCurrent)

					elif plotSettings['typeOfPlot'] == 'UR1_tau':

						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)

						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + ' | $Q_y$=' + str(case.ForceMagnitude)+'N', **plotSettings['title'])
						scatterHandles[keyCurrent] = plotUR1_tau(case, plotSettings, keyCurrent, axDict[keyCurrent], counterNperKey, scatterHandles)
						counterNperKey[keyCurrent] += 1

						if not keyCurrent in keysWith_UR1_tau_plot: keysWith_UR1_tau_plot.append(keyCurrent)

					elif plotSettings['typeOfPlot'] == 'plotU2_z_LastTau':

						flagDict, axDict = figureInitialization(flagDict, axDict, keyCurrent, plotSettings)
						flagDict[keyCurrent] = True
						axDict[keyCurrent].set_title(plotSettings['xLabel'][keyCurrent] + '='+str(getattr(case, keyCurrent))+', $Q_y$=' + str(case.ForceMagnitude)+'N'+'/last frame', **plotSettings['title'])
						plotU2_z_LastTau(case, plotSettings, keyCurrent, axDict[keyCurrent])

	if plotSettings['typeOfPlot'] == 'UR1_tau':

		for key in tuple(keysWith_UR1_tau_plot):

			axDict[key].legend(handles = scatterHandles[key], **plotSettings['legend'])


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

def plotU2_z_LastTau(classOfData, plotSettings, attr, ax):

	ax.set_xlabel('$z / C_3$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	indexForMaxFrame = classOfData.framesCount.index(max(classOfData.framesCount))
	j=0
	# pdb.set_trace()
	for x in classOfData.dataFrames[indexForMaxFrame].xPosForU2:

		# pdb.set_trace()
		ax.plot(classOfData.dataFrames[indexForMaxFrame].dataU2OverX[j].zOverC3, classOfData.dataFrames[indexForMaxFrame].dataU2OverX[j].u2_zOverC3, label='x='+str(round(x,2)), **plotSettings['line'])
		j += 1

	ax.legend(**plotSettings['legend'])


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


def plotUR1_tau(classOfData, plotSettings, attr, ax, counterNperKey, scatterHandles):

	ax.set_xlabel('$Q_{frame} / Q_{total}$', **plotSettings['axes_x'])
	ax.set_ylabel(plotSettings['yLabel'], **plotSettings['axes_y'])

	i = 0
	for fraction in classOfData.framesFraction:

		#Find index for current frame, classOfData.dataFrames contains unsorted data
		indexFrame = classOfData.framesCount.index(i)

		dataThisFrame = classOfData.dataFrames[indexFrame]

		#To obtain twist from U2 difference
		indexForMaxX = dataThisFrame.xPosForU2.index(max(dataThisFrame.xPosForU2))

		if plotSettings['meanOption']:

			meanTwist = np.mean([dataThisFrame.ur1_xOverL_up[-1], dataThisFrame.ur1_xOverL_dn[-1], dataThisFrame.twistFromU2[indexForMaxX]])
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

	if ((abs(a - b) / b)*100) > 5: #If error > 5%
		raise ValueError('ERROR: The error in twist for the linear solution is bigger than 5%, having it been obtained from different parts of the beam')	
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

	return scatterHandles[attr]