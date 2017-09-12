import os, re, sys

from moduleParametricStudy import *
from moduleCommon import *
import pdb #pdb.set_trace()
import warnings

#### PLOTTING OPTIONS ####

#Plotting options
axes_label_x  = {'size' : 18, 'weight' : 'bold', 'verticalalignment' : 'top', 'horizontalalignment' : 'center'} #'verticalalignment' : 'top'
axes_label_y  = {'size' : 18, 'weight' : 'bold', 'verticalalignment' : 'bottom', 'horizontalalignment' : 'center'} #'verticalalignment' : 'bottom'
text_title_properties = {'weight' : 'bold', 'size' : 18}
axes_ticks = {'labelsize' : 14}
line = {'linewidth' : 2, 'markersize' : 5}
scatter = {'linewidths' : 2}
legend = {'fontsize' : 16, 'loc' : 'best'}
grid = {'alpha' : 0.7}
colors = ['k', 'b', 'y', 'm', 'r', 'c','k', 'b', 'y', 'm', 'r', 'c','k', 'b', 'y', 'm', 'r', 'c','k', 'b', 'y', 'm', 'r', 'c']
markers = ['o', 'v', '^', 's', '*', '+']
linestyles = ['-', '--', '-.', ':'] 

#x labels
xLabel={'N' : 'Number of unit cells in spanwise direction',
        'M' : 'Number of unit cells in transversal direction',
        'r' : 'Node radius (mm)',
        'B' : 'Node depth (mm)',
        'L' : 'Ligament half length (mm)',
        'C3' : 'C-box length in the chordwise direction (mm)',
        'Cbox_t' : 'C-box wall thickness, $t_{C}$ (mm)', 
        'rib_t' : 'Rib thickness, $t_{rib}$ (mm)', 
        'rib_a' : 'Rib dimension frame width, $a$ (mm)', 
        'eOverB' : 'Chiral ligament eccentricity, $e/B$ (%)', 
        'tChiral' : 'Chiral lattice section thickness, $t_{chiral}$ (mm)',
        'forceMagnitude' : 'Applied force magnitude per unit of length (N/mm)',
        'forceXStart' : 'Initial x-coordinate of distributed force (mm)',
        'forceXEnd' : 'Final x-coordinate of distributed force (mm)',
        'forceXn' : 'Number of points where the force is applied',
        'forceZPos' : 'Load point position in the chordwise direction',
        'innerRibs_n' : 'Number of ribs in the inner part of the wing box',
        'courseSize' : 'Course mesh size',
        'fineSize' : 'Fine mesh size',
        'maxTimeIncrement' : 'Maximum time increment allowed',
        'initialTimeIncrement' : 'Initial time increment',
        'maxNumInc' : 'Maximum number of increments in a step',
        'damp' : 'Artificial constant damping factor',
        'additionalBC' : 'Boundary condition at the border',
        'ForceMagnitude' : 'Magnitude of the force applied',
        'dofContraint' : 'Degrees of freedom constrained',
        'typeLoad' : 'Type of load introduction method'}

plotSettings = {'xLabel':xLabel,'axes_x':axes_label_x,'axes_y':axes_label_y, 'title':text_title_properties,
                'axesTicks':axes_ticks, 'line':line, 'legend':legend, 'grid':grid, 'scatter':scatter,
                'colors' : colors, 'markers' : markers, 'linestyles' : linestyles}

#### INITIALIZE FOLDERS ####

print('Parametric study processing started...')

#Get working directory
cwd = os.getcwd()

#Read postProc folder name from CMD
CMDoptionsDict = {}
CMDoptionsDict = readCMDoptions(sys.argv[1:], CMDoptionsDict)

#Results options
plotSettings['meanOption'] = CMDoptionsDict['plotMean'] #Plot mean of the twist (True) or plot twist obtained from different parts of the model (upper and lower flanges; and difference on U2 for upper flange)

#Move to post-processing directory
globalChangeDir(cwd, '-'+CMDoptionsDict['postProcFolderName'])
postProcFolder = os.getcwd()

#Read parameter study definition file
studyDefDict = importParametricStudyDeffile('parametricStudyDef.txt')

#### IMPORT DATA ####
data = []
for file in os.listdir(postProcFolder):
    globalChangeDir(cwd, '-'+CMDoptionsDict['postProcFolderName'])
    if file.endswith('inputAbaqus_nonlinear.txt'):
        
        #Create case study class where store all the results obtained from Abaqus at termination of its computation
        temp = caseStudy(int(file[:-26])) #file[:-26] - Returns the index

        temp.importInputData(file)

        #Import data for case
        globalChangeDir(cwd, '-'+CMDoptionsDict['postProcFolderName']+'-'+str(temp.id))
        postProcFolderForCase = os.getcwd()

        #Show warning message when there are not result files
        if not os.listdir(postProcFolderForCase): # equals to: if os.listdir(postProcFolderForCase) == []
            warnings.warn('-> No result files found for iteration '+str(temp.id))

        #NONLINEAR

        #Obtain information of the fraction of load applied at each frame
        temp.framesFraction = temp.obtainFrameLoadFractionInfo('frameInfo.txt')

        #Obtain data of the external work vs static dissipation energy given by Abaqus in (N * mm)
        if temp.damp != '0.0':
            temp.obtainEnergyData('extWork_stab.rpt')

        #Obtain final values of twist
        if CMDoptionsDict['collectTwist']:
            try:
                temp.obtainTwistData('twist.rpt', '')
                temp.obtainTwistData('linear_twist.rpt', 'linear_')
            except FileNotFoundError as e:
                print('-> Data for twist on selected points not found for iteration: '+str(temp.id))
                plotSettings['twistData'] = False
            else:
                print('-> Data for twist on selected points loaded for iteration: '+str(temp.id))
                plotSettings['twistData'] = True
        else:
            plotSettings['twistData'] = False

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
        try:
            temp.import_data_from_path('linear_ur1_xOverL.rpt', 'linear_ur1', 'xOverL')
            temp.import_data_from_path('linear_u2_zOverC3.rpt', 'linear_u2', 'zOverC3')
            temp.import_data_from_path('linear_u2_xOverL.rpt', 'linear_u2', 'xOverL') #Not essential for later calculations
        except FileNotFoundError as e:
            print('-> Linear results not found for iteration: '+str(temp.id))
            plotSettings['plotLinear'] = False
        else:
            plotSettings['plotLinear'] = True
        
        #Store global class
        data += [temp]

        #Return to original working folder
        globalChangeDir(cwd, '')

#Return to original working folder
globalChangeDir(cwd, '')

print('-> Data loaded...')

#Print summary table
table_sum = tableOutput('Simulation summary', ['parameter', 'value', 'max Q_fr/Q_to', 'frames', 'f_mesh', 'c_mesh', 'damp', 'forceMagnitude', 'wingBoxLength', 'wingBoxHeight'])
for case in data:
    for (keyCurrent, rangeCurrent) in studyDefDict.items():
        if studyDefDict[keyCurrent][0] <= case.id <= studyDefDict[keyCurrent][1]:
            table_sum.printRow([keyCurrent, getattr(case, keyCurrent), float(max(case.framesFraction)), int(max(case.framesCount)), case.fineSize, case.courseSize, case.damp, case.ForceMagnitude, float(case.wingBoxLength), 2*88.3176086632785*(float(case.M)-1)])

#### PLOTTING ####

# NONLINEAR PLOTS
# plotSettings['typeOfPlot'] = 'UR1_frame'
# plotSettings['yLabel'] = 'Angular rotation (deg)'
# caseDistintion(data, studyDefDict, plotSettings, CMDoptionsDict)
if 'energy' in CMDoptionsDict['plotOptString']:
    plotSettings['typeOfPlot'] = 'energy'
    plotSettings['yLabel'] = 'Energy (kJ)'
    caseDistintion(data, studyDefDict, plotSettings, CMDoptionsDict, [])

if 'UR1_tau' in CMDoptionsDict['plotOptString']:
    plotSettings['typeOfPlot'] = 'UR1_tau'
    plotSettings['yLabel'] = '$\phi_{\mathrm{tip}} (\mathrm{deg})$'
    table = tableOutput('Angular rotation at tip UR1 (deg)', ['parameter', 'value', 'max UR1', 'error UR1 (%)', 'UR1, linear', 'error UR1, linear (%)'])
    caseDistintion(data, studyDefDict, plotSettings, CMDoptionsDict, table)

if 'U2_z' in CMDoptionsDict['plotOptString']:
    plotSettings['typeOfPlot'] = 'plotU2_z_LastTau'
    plotSettings['yLabel'] = 'Vertical displacement $U_2$ (mm)'
    table = tableOutput('Vertical displacement U2 (mm)', ['parameter', 'value', 'max U2', 'zOverC3_maxU2', 'xOverL_maxU2','min U2', 'zOverC3_minU2', 'xOverL_minU2'])
    caseDistintion(data, studyDefDict, plotSettings, CMDoptionsDict, table)


# plt.show(block = True)

# for i in plt.get_fignums(): #NOT WORKING
#     plt.figure(i)
#     plt.savefig('figure%d.png' % i)

#Return to main working directory
os.chdir(cwd)

if not isUnix() and CMDoptionsDict['showFigures']:
    plt.show(block = True)

print('---> Parametric study finished')