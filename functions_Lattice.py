from abaqus import *
from abaqusConstants import *
import visualization
import os
import re
import sys
import shutil
import subprocess as sp
import numpy
import math
import regionToolset
import copy
# import numpy as np

##(jobname_title,index,model,initialIncrement,maxInc,varsFieldOut,varsHistoryOut,dampMagnit)  ####
##
##jobName=jobname_title
##stepNo=index
##currentModel=model
##initialIncrement=initialIncrement
##maxIncr=maxInc
##varsFieldOut=varsFieldOut
##varsHistoryOut=varsHistoryOut
##stabMag=dampMagnit

 
def aeroLoopStepsNL(jobName,stepNo,currentModel,initialIncrement,maxIncr,varsFieldOut,varsHistoryOut,stabMag):
    
    currentStep = 'Step-' + str(stepNo)# + '-' + str(jobNo)

    if stepNo == 1:
        previousStep = 'Initial'
    else:
        previousStep = 'Step-' + str(stepNo-1)# + '-' + str(jobNo)
        previousJobName = jobName+str(stepNo-1)#+ '-' + str(jobNo)
        #-----------------------------------------------------
        # define where to restart analysis for higher steps
        currentModel.setValues( restartJob = previousJobName, restartStep = previousStep )
    
    #-----------------------------------------------------
    # create next step and restart
    stepLoop = currentModel.StaticStep(
           description='weakly coupled static aeroelastic analysis with xfoil', 
           continueDampingFactors=False,
           stabilizationMagnitude=stabMag, 
           stabilizationMethod=DAMPING_FACTOR,#DISSIPATED_ENERGY_FRACTION
           initialInc=initialIncrement,#initialIncrement,
           maxInc=maxIncr,
           minInc=10e-10,
           maxNumInc=3000, 
           name=currentStep, 
           nlgeom=ON, 
           previous=previousStep )
            
    #-----------------------------------------------------
    # restart
    stepLoop.Restart(frequency=1, numberIntervals=0, 
        overlay=ON, timeMarks=OFF)

    #-----------------------------------------------------
    # limit output requests
    if stepNo == 1:
        try:
            del currentModel.fieldOutputRequests['F-Output-1']
        except:
            pass
        currentModel.FieldOutputRequest(createStepName=currentStep, 
            name='FieldOutput', variables=varsFieldOut)
        # history output
        try:
            del currentModel.historyOutputRequests['H-Output-1']
        except:
            pass
        currentModel.HistoryOutputRequest(createStepName=currentStep, 
            name='HistoryOutput', variables=varsHistoryOut)			
			
			
#######################################################################
