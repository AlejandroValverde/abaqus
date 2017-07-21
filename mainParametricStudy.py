import os, re

from moduleParametricStudy import *
from moduleCommon import *
import pdb #pdb.set_trace()
import warnings

#### PLOTTING OPTIONS ####

#Results options
meanOptionalFlag = True #Plot mean of the twist (True) or plot twist obtained from different parts of the model (upper and lower flanges; and difference on U2 for upper flange)

#Plotting options
axes_label_x  = {'size' : 14, 'weight' : 'bold', 'verticalalignment' : 'center', 'horizontalalignment' : 'center'} #'verticalalignment' : 'top'
axes_label_y  = {'size' : 14, 'weight' : 'bold', 'verticalalignment' : 'center', 'horizontalalignment' : 'center'} #'verticalalignment' : 'bottom'
text_title_properties = {'weight' : 'bold', 'size' : 18}
axes_ticks = {'size' : 12}
line = {'linewidth' : 1, 'markersize' : 5}
scatter = {'linewidths' : 2}
legend = {'fontsize' : 14, 'loc' : 'best'}
grid = {'alpha' : 0.7}
colors = ['b', 'k', 'y', 'm', 'r', 'c']
markers = ['o', 'v', '^', 's', '*', '+']

#x labels
xLabel={'N' : 'Number of unit cells in transversal direction', 
        'Cbox_t' : 'C-box wall thickness, $t_{C}$ (mm)', 
        'rib_t' : 'Rib thickness, $t_{rib}$ (mm)', 
        'rib_a' : 'Rib dimension frame width, $a$ (mm)', 
        'eOverB' : 'Chiral ligament eccentricity, $e/B$ (%)', 
        'tChiral' : 'Chiral lattice section thickness, $t_{chiral}$ (mm)',
        'forceMagnitude' : 'Applied force magnitude per unit of length (N/mm)',
        'forceXStart' : 'Initial x-coordinate of distributed force (mm)',
        'forceXEnd' : 'Final x-coordinate of distributed force (mm)',
        'forceXn' : 'Number of points where the force is applied',
        'forceZ3' : 'Z-coordinate of distributed force (mm)',
        'innerRibs_n' : 'Number of ribs in the inner part of the wing box',
        'courseSize' : 'Course mesh size',
        'fineSize' : 'Fine mesh size',
        'maxTimeIncrement' : 'Maximum time increment allowed',
        'initialTimeIncrement' : 'Initial time increment',
        'maxNumInc' : 'Maximum number of increments in a step'}

plotSettings = {'xLabel':xLabel,'axes_x':axes_label_x,'axes_y':axes_label_y, 'title':text_title_properties,
                'axesTicks':axes_ticks, 'line':line, 'legend':legend, 'grid':grid, 'scatter':scatter,
                'colors' : colors, 'markers' : markers, 'meanOption' : meanOptionalFlag}

#### INITIALIZE FOLDERS ####

warnings.simplefilter('default', UserWarning) # print the first occurrence of matching warnings for each location where the warning is issued

print('Parametric study processing started...')

#Get working directory
cwd = os.getcwd()

#Move to post-processing directory
globalChangeDir(cwd, '-postProc')
postProcFolder = os.getcwd()

#Read parameter study definition file
studyDefDict = importParametricStudyDeffile('parametricStudyDef.txt')

#### IMPORT DATA ####
data = []
for file in os.listdir(postProcFolder):
    globalChangeDir(cwd, '-postProc')
    if file.endswith('inputAbaqus_nonlinear.txt'):
        
        #Create case study class where store all the results obtained from Abaqus at termination of its computation
        temp = caseStudy(int(file[:-26])) #file[:-26] - Returns the index

        temp.importInputData(file)

        #Import data for case
        globalChangeDir(cwd, '-postProc-'+str(temp.id))
        postProcFolderForCase = os.getcwd()

        #Show warning message when there are not result files
        if not os.listdir(postProcFolderForCase): # equals to: if os.listdir(postProcFolderForCase) == []
            warnings.warn('-> No result files found for iteration '+str(temp.id))

        #LINEAR
        if temp.typeAnalysis == 'linear': #NOT IN USE

            temp.importReactionForce()

            temp.import_data_from_path(str(temp.id)+'-'+kindY+'_'+kindX+'.rpt', 'ur1', 'xOverL')

            temp.import_data_from_path(str(temp.id)+'-'+kindY+'_'+kindX+'.rpt', 'u2', 'zOverC3')

            temp.import_data_from_path(str(temp.id)+'-'+kindY+'_'+kindX+'.rpt', 'u2', 'xOverL')

            temp.calculateK()

        #NONLINEAR
        elif temp.typeAnalysis == 'nonlinear':

            #Obtain information of the fraction of load applied at each frame
            temp.framesFraction = temp.obtainFrameLoadFractionInfo('frameInfo.txt')

            #For each frame
            dataFrames = []
            framesCount = []
            for file2 in os.listdir(postProcFolderForCase):

                if file2.startswith('ur1_up'):

                    #Initialize class for frame
                    tempPerFrame = dataPerFrame(int(file2.replace('.rpt','')[12:])) #int(file2.replace('.rpt','')[12:]) returns the frame ID

                    #Import data from upper flange
                    tempPerFrame.import_data_from_path(file2, 'ur1', 'xOverL_'+file2[4:6])

                    #For lower flange
                    file2_down = file2.replace('up','dn')
                    tempPerFrame.import_data_from_path(file2_down, 'ur1', 'xOverL_'+file2_down[4:6])

                    #For u2 vs x
                    dataU2OverX = []
                    dataxPos = []
                    dataTwistFromU2 = []
                    for file3 in os.listdir(postProcFolderForCase):
                        if file3.startswith('u2_frame'+str(tempPerFrame.frameID)+'_x'):

                            tempPerFramePerX = dataPerFramePerX(float(file3.replace('u2_frame'+str(tempPerFrame.frameID)+'_x','')[:-4]))

                            tempPerFramePerX.import_data_from_path(file3, 'u2', 'zOverC3')

                            twistFromU2 = (tempPerFramePerX.u2_zOverC3[-1] - tempPerFramePerX.u2_zOverC3[0])/float(temp.C3)

                            dataTwistFromU2 += [twistFromU2]

                            dataU2OverX += [tempPerFramePerX]

                            dataxPos += [tempPerFramePerX.xPos]

                    #Save results into class
                    tempPerFrame.dataU2OverX = dataU2OverX

                    tempPerFrame.xPosForU2 = dataxPos

                    tempPerFrame.twistFromU2 = dataTwistFromU2

                    #Save class into global list for this case
                    dataFrames += [tempPerFrame]

                    framesCount += [tempPerFrame.frameID]


            temp.dataFrames = dataFrames
            temp.framesCount = framesCount


            #For the corresponding linear simulation
            temp.import_data_from_path('linear_ur1_xOverL.rpt', 'linear_ur1', 'xOverL')
            temp.import_data_from_path('linear_u2_zOverC3.rpt', 'linear_u2', 'zOverC3')
            temp.import_data_from_path('linear_u2_xOverL.rpt', 'linear_u2', 'xOverL') #Not essential for later calculations

            #Store global class
            data += [temp]

        #Return to original working folder
        globalChangeDir(cwd, '')

#### PLOTTING ####

##Plot reaction force (RF-2) as a function of parameter values
plotSettings['yLabel'] = 'Reaction force, $R_y$ (N)'
plotSettings['typeOfPlot'] = 'RF'
caseDistintion(data, studyDefDict, plotSettings)

#Plot initial stiffness (K) as a function of parameter values
plotSettings['yLabel'] = 'Initial stiffness, $K$ (N/mm)'
plotSettings['typeOfPlot'] = 'K'
caseDistintion(data, studyDefDict, plotSettings)

#Plot UR-1 along the wing box length
plotSettings['typeOfPlot'] = 'UR1'
plotSettings['yLabel'] = 'Angular rotation $UR_1$ (deg)'
caseDistintion(data, studyDefDict, plotSettings)

#Plot vertical displacement U2 along the wing box chordwise direction 
plotSettings['typeOfPlot'] = 'U2_z'
plotSettings['yLabel'] = 'Vertical displacement $U_2$ (mm)'
caseDistintion(data, studyDefDict, plotSettings)

#Plot vertical displacement U2 along the wing box spanwise direction 
plotSettings['typeOfPlot'] = 'U2_x'
plotSettings['yLabel'] = 'Vertical displacement $U_2$ (mm)'
caseDistintion(data, studyDefDict, plotSettings)

# NONLINEAR PLOTS
# plotSettings['typeOfPlot'] = 'UR1_frame'
# plotSettings['yLabel'] = 'Angular rotation (deg)'
# caseDistintion(data, studyDefDict, plotSettings)

plotSettings['typeOfPlot'] = 'UR1_tau'
plotSettings['yLabel'] = 'Angular rotation (deg)'
caseDistintion(data, studyDefDict, plotSettings)

plotSettings['typeOfPlot'] = 'plotU2_z_LastTau'
plotSettings['yLabel'] = 'Vertical displacement $U_2$ (mm)'
caseDistintion(data, studyDefDict, plotSettings)


# plt.show(block = True)

# for i in plt.get_fignums(): #NOT WORKING
#     plt.figure(i)
#     plt.savefig('figure%d.png' % i)

#Return to main working directory
os.chdir(cwd)

plt.show(block = not isUnix())

print('-> Parametric study finished')