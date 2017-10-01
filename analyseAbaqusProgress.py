import os
import sys
import pdb #pdb.set_trace()
import math
import getopt

try:
	import matplotlib.pyplot as plt
	from matplotlib.ticker import MaxNLocator
except Exception as e:
	print('->>>> Failed to import matplotlib')
else:
	print('Matplotlib imported successfully')

def readCMDoptionsMainAbaqusParametric(argv, CMDoptionsDict):

	short_opts = "v:"
	long_opts = ["velocity="]
	try:
		opts, args = getopt.getopt(argv,short_opts,long_opts)
	except getopt.GetoptError:
		raise ValueError('ERROR: Not correct input to script')

	# check input
	if len(opts) != len(long_opts):
		raise ValueError('ERROR: Invalid number of inputs')	

	for opt, arg in opts:

		if opt in ("-v", "--velocity"):

			CMDoptionsDict['V'] = float(arg)

	return CMDoptionsDict

class AbaqusProgressFile(object):
	"""docstring for AbaqusProgressFile"""
	def __init__(self, arg):
		super(AbaqusProgressFile, self).__init__()
		self.id = arg

	def loadOutput(self):

		file = open('Output-Job-54_20_1-'+str(self.id)+'.txt', 'r')

		lines = file.readlines()

		# pdb.set_trace()

		nameParater = lines[0]

		nameParater = nameParater.replace('\r\n','')

		nameParater = nameParater.replace('\n','')

		parameters = nameParater.split('\t')

		self.step = float(parameters[2])
		self.twist = float(parameters[4])
		self.changeTwist = float(parameters[5])
		self.RF = float(parameters[6])

class dataFromAbaqusInFolder(object):
	"""docstring for dataFromAbaqusInFolder"""
	def __init__(self, folderName):
		super(dataFromAbaqusInFolder, self).__init__()
		self.folderName = folderName

		print(self.folderName)


	def loadData(self):

		data = []
		steps = []
		id_unsorted = []
		twist_unsorted = []
		RF_unsorted = []
		id_sorted = []
		twist_sorted = []
		RF_sorted = []		
		i = -1
		for file in os.listdir(self.folderName):

			if file.startswith('Output-Job'):

				file = file.replace('.txt','')
				file = file.replace('Output-Job-54_20_1-','')

				fileData = AbaqusProgressFile(int(file))

				fileData.loadOutput()

				if not fileData.step in steps:
					steps.append(fileData.step)
					
					i += 1
					id_unsorted.append([])
					twist_unsorted.append([])
					RF_unsorted.append([])

				id_unsorted[i].append(fileData.id)
				twist_unsorted[i].append(fileData.twist)
				RF_unsorted[i].append(fileData.RF)

				data.append(fileData)

		for id_j in range(len(id_unsorted)):

			c, d, e = (list(t) for t in zip(*sorted(zip(id_unsorted[id_j], twist_unsorted[id_j], RF_unsorted[id_j])))) #Follow the order of id_unsorted[id_j]
			id_sorted.append(c) 
			twist_sorted.append(d) 
			RF_sorted.append(e) 


		self.dataForAllJobs = data
		self.steps = sorted(steps)
		self.steps_unsorted = steps
		self.nsteps = len(steps)
		self.id_unsorted = id_unsorted
		self.id_sorted = id_sorted
		self.twist_unsorted = twist_unsorted
		self.twist_sorted = twist_sorted
		self.RF_unsorted = RF_unsorted
		self.RF_sorted = RF_sorted


################################
#Script

CMDoptionsDict = {}
CMDoptionsDict = readCMDoptionsMainAbaqusParametric(sys.argv[1:], CMDoptionsDict)

#Flight data
V = CMDoptionsDict['V']

#Load data
pwd = os.getcwd()
os.chdir(pwd+'\output_step0-1_V30_s3')
data = dataFromAbaqusInFolder(os.getcwd())
data.loadData()
os.chdir(pwd+'\output_step0-1_V30_s3_lin')
data1 = dataFromAbaqusInFolder(os.getcwd())
data1.loadData()
# pdb.set_trace()

#Plotting options
axes_label_x  = {'size' : 18, 'weight' : 'bold', 'verticalalignment' : 'top', 'horizontalalignment' : 'center'} #'verticalalignment' : 'top'
axes_label_y  = {'size' : 18, 'weight' : 'bold', 'verticalalignment' : 'bottom', 'horizontalalignment' : 'center'} #'verticalalignment' : 'bottom'
text_title_properties = {'weight' : 'bold', 'size' : 14}
axes_ticks = {'labelsize' : 16}
line = {'linewidth' : 2, 'markersize' : 10}
scatter = {'linewidths' : 2}
legend = {'fontsize' : 18, 'loc' : 'best'}
grid = {'alpha' : 0.7}
colors = ['k', 'b', 'y', 'm', 'r', 'c','k', 'b', 'y', 'm', 'r', 'c','k', 'b', 'y', 'm', 'r', 'c','k', 'b', 'y', 'm', 'r', 'c']
markers = ['o', 'v', '^', 's', '*', '+']
linestyles = ['-', '--', '-.', ':'] 

plotSettings = {'axes_x':axes_label_x,'axes_y':axes_label_y, 'title':text_title_properties,
                'axesTicks':axes_ticks, 'line':line, 'legend':legend, 'grid':grid, 'scatter':scatter,
                'colors' : colors, 'markers' : markers, 'linestyles' : linestyles}

#Plotting

###################
attributes = ('twist', 'RF')
ylabels = ('$\phi_{\mathrm{tip}}$ [deg]', '$RF_{\mathrm{root}}$ [N]')
# Sim evolution
for atr, ylabel in zip(attributes, ylabels):
	figure, axss = plt.subplots(nrows = 2, ncols = 5, sharey = True)
	figure.set_size_inches(10*math.pow(figure.get_dpi()/100, -1), 6*math.pow(figure.get_dpi()/100, -1), forward=True)
	# mapIndex = [[0,0], [0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [1, 5]]
	mapIndex = [[0,0], [0, 1], [0, 2], [0, 3], [0, 4], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4]]
	if hasattr(axss, '__iter__'):#If iterable
		i = 0
		breakFlag = False
		axss[0][0].set_ylabel(ylabel, **plotSettings['axes_y'])
		axss[1][0].set_ylabel(ylabel, **plotSettings['axes_y'])
		for axs in axss:
			for ax in axs:
				# pdb.set_trace()
				ax.grid(which='both', **plotSettings['grid'])
				ax.tick_params(axis='both', **plotSettings['axesTicks'])

				ax.set_xlabel('iter', **plotSettings['axes_x'])

				ax.set_title('Step:'+str(data.steps[i]) + ', V='+str(round(math.sqrt(data.steps[i])*V, 2)), **plotSettings['title'])

				i += 1

				if i == data.nsteps:
					breakFlag = True
					break

			if breakFlag:
				break

	else:

		axs.grid(which='both', **plotSettings['grid'])
		axs.tick_params(axis='both', **plotSettings['axesTicks'])
		figure.set_size_inches(10, 6, forward=True)

		axs.set_xlabel('iteration', **plotSettings['axes_x'])
		axs.set_ylabel('$\phi_{\mathrm{tip}}$ [deg]', **plotSettings['axes_y'])

		axs.set_title('Step: '+str(data.steps[0]) + ', V='+str(round(math.sqrt(data.steps[0])*V, 2)), **plotSettings['title'])
	for fileOfData in data.dataForAllJobs:
		if hasattr(axs, '__iter__'):#If iterable
			axss[mapIndex[data.steps.index(fileOfData.step)][0],mapIndex[data.steps.index(fileOfData.step)][1]].plot(int(fileOfData.id), getattr(fileOfData, atr), marker = 'o', c = plotSettings['colors'][0], **plotSettings['line']) #/ max(classOfData.framesCount)
			axss[mapIndex[data.steps.index(fileOfData.step)][0],mapIndex[data.steps.index(fileOfData.step)][1]].xaxis.set_major_locator(MaxNLocator(integer=True))
		else:
			axs.plot(int(fileOfData.id), getattr(fileOfData, atr), marker = 'o', c = plotSettings['colors'][0], **plotSettings['line']) #/ max(classOfData.framesCount)

	####################
	#Twist vs speed

	figure, ax = plt.subplots(nrows = 1, ncols = 1, sharey = True)
	figure.set_size_inches(10*math.pow(figure.get_dpi()/100, -1), 6*math.pow(figure.get_dpi()/100, -1), forward=True)
	ax.set_ylabel(ylabel, **plotSettings['axes_y'])
	ax.set_xlabel('$V$ [m/s]', **plotSettings['axes_x'])
	ax.grid(which='both', **plotSettings['grid'])
	ax.tick_params(axis='both', **plotSettings['axesTicks'])
	# ax.set_title('Step:'+str(data.steps[i]) + ', V='+str(round(math.sqrt(data.steps[i])*V, 2)))
	i = 0
	for step in data.steps_unsorted:
		ax.plot(math.sqrt(step)*V, getattr(data, atr+'_sorted')[i][-1], marker = 'o', c = plotSettings['colors'][0], **plotSettings['line']) #/ max(classOfData.framesCount
		i += 1

	j = 0
	for step in data1.steps_unsorted:
		ax.plot(math.sqrt(step)*V, getattr(data1, atr+'_sorted')[j][-1], marker = '+', c = plotSettings['colors'][0], **plotSettings['line']) #/ max(classOfData.framesCount
		j += 1

	handle1 = plt.Line2D([],[], color=plotSettings['colors'][0], marker='o', linestyle='', label='nonlinear')
	handle2 = plt.Line2D([],[], color=plotSettings['colors'][0], marker='+', linestyle='', label='linear')
	ax.legend(handles = [handle1, handle2], **plotSettings['legend'])

plt.show(block = True)