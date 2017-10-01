########################################################################
# MAIN FILE AEROELASTIC SIMULATION
########################################################################

#session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE) # only necessary when expanding the file in GUI

import math			# necessary for certain mathematical functions
import os			# necessary for creation of folders

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import string
import numpy as np
# from parameters_Dom import N, c_0#, b #Loaded from input file

from functions_Urban_sim import *

## From Parameters54_20
iSIM=54
PositionFrontSpar=2.000000e-01
PositionRearSpar=3.000000e-01
FaserWinkel=85
thickness_bucklingSpar=1.500000e-04

elements_inWeb=20

# DEFINE DIRECTORIES GRIFFIN
#start_dir = '/cluster/home/runkelf/opt/code'
##start_dir = folder_path=os.getcwd()
#folder_path = '/cluster/home/runkelf/opt/'

#start_dir = 'D:/faselu/ftero_opt_local/code'
start_dir = folder_path=os.getcwd()
#folder_path = 'D:/faselu/ftero_opt_local'
[folder_path,folder_path_ende]=os.path.split(start_dir)
xfoil_path = folder_path+'/xfoil6.96/bin/xfoil.exe'

##################################################
jobname = 'Wing-stat'
res_name = '/xy_result.txt'
res_name_spars = '/SparVorne.txt' 
res_name_spars_hinten ='/SparHinten.txt'

daempfung=[50000]#[500,600,700]
jobNumber=1

geklappt = open(folder_path+'/Geklappt'+str(iSIM)+'_'+str(elements_inWeb)+'.txt','w')
geklappt.write('No'+'\n')
geklappt.close()

Nein='No'

#iSIM = i

for trials in range(0,len(daempfung)):
    if os.path.exists(folder_path+'/Damping'+str(iSIM)+'_'+str(elements_inWeb)+'.txt'):
        os.remove(folder_path+'/Damping'+str(iSIM)+'_'+str(elements_inWeb)+'.txt')
    
    damping = open(folder_path+'/Damping'+str(iSIM)+'_'+str(elements_inWeb)+'.txt','a')
    damping.write(str(daempfung[trials]))
    damping.close()
    
    geklappt = open(folder_path+'/Geklappt'+str(iSIM)+'_'+str(elements_inWeb)+'.txt','a')
    lines = [line.rstrip('\n') for line in open(folder_path+'/Geklappt'+str(iSIM)+'_'+str(elements_inWeb)+'.txt')]
    geklappt.close()

    if lines[0]==Nein:
		LPF=1
		os.chdir(folder_path+'/code')
		WingModel = mdb.Model(name = 'Model-SparAngle-'+str(int(jobNumber)))
		# execfile('RunModelling.py')
		execfile('codeForWing_sim.py')
##		# DEFINE RESULT FOLDER
		if not os.path.exists(folder_path+'/Results-'+str(iSIM)+'_'+str(elements_inWeb)):
			os.makedirs(folder_path+'/Results-'+str(iSIM)+'_'+str(elements_inWeb))
		os.chdir(folder_path+'/Results-'+str(iSIM)+'_'+str(elements_inWeb))

		#Clean all the previous files from execution
		for f in os.listdir(os.getcwd()):
		    # if f.startswith(jobNameComplete.replace('nonlinear','linear')) or f.startswith('abaqus.rpy'):
			os.remove(f)
		
		# mdb.models['Model-SparAngle-'+str(int(jobNumber))].sections['Spar-Rear_GLASS_ONLY'].setValues(idealization=NO_IDEALIZATION, integrationRule=SIMPSON, layup=(SectionLayer(thickness=0.001*thickness_bucklingSpar, orientAngle=FaserWinkel, material='Glass-Fabric', plyName='ply1111'), ), preIntegrate=OFF, symmetric=False)
		execfile(start_dir+'/AeroelasticLoopCP_sim.py')





