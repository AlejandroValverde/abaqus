import os
import sys
import getopt

def readCMDoptionsMainAbaqusParametric(argv, CMDoptionsDict):

	short_opts = "o:"
	long_opts = ["options="]
	try:
		opts, args = getopt.getopt(argv,short_opts,long_opts)
	except getopt.GetoptError:
		raise ValueError('ERROR: Not correct input to script')

	# check input
	if len(opts) != len(long_opts):
		raise ValueError('ERROR: Invalid number of inputs')	

	for opt, arg in opts:

		if opt in ("-o", "--options"):

			if arg.lower() in ('script'):
				CMDoptionsDict['type'] = 'script'
			elif arg.lower() in ('nogui'):
				CMDoptionsDict['type'] = 'noGUI'

	return CMDoptionsDict

###############

CMDoptionsDict = {}
CMDoptionsDict = readCMDoptionsMainAbaqusParametric(sys.argv[1:], CMDoptionsDict)

folder_path= os.getcwd()
os.chdir(folder_path + '/code')
print(os.getcwd())

Runs=[54]
Elements=[20]

for RunNb in range(1,len(Runs)+1):
    for ElemNb in range(1,len(Elements)+1):
        # commandLine='abaqus cae noGUI=Parameters'+str(Runs[RunNb-1])+'_'+str(Elements[ElemNb-1])+'.py'
        if CMDoptionsDict['type'] == 'script':
        	commandLine='abaqus cae script=main_sim.py'
        elif CMDoptionsDict['type'] == 'noGUI':
        	commandLine='abaqus cae noGUI=main_sim.py'
        print(commandLine)
        os.system(commandLine)

    


