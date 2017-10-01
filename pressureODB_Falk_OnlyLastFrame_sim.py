#######################################################################
#
#   CALCULATE LIFT AND DRAG BY APPLIED PRESSURE
#
#######################################################################

import visualization
import math

#index = index  # step-Nr
#alpha = 7

primaryJobName = 'Wing-stat_'+str(jobNumber)+'-' + str(index)#+' - Copy'
##primaryJobName = 'Wing-stat_'+str(jobNumber) +'-'+str(int(1000*thickness_bucklingSpar[digge]))+'-'+str(int(FaserWinkel[faserw]))+ '-' + str(index)

odbName = primaryJobName + '.odb'
odb = visualization.openOdb(odbName)

stepName = 'Step-' + str(index)
#assembly = mdb.Model(name='Model-SparAngle-5').rootAssembly


iterations=odb.steps[stepName].frames[-1]
##print(frameValueOld)
##frameValueOld=frameValueOld+0.049
twistCur = [0.0]*N
for i in range(N):
    EiD=iterations.frameId
    twistCur[i] = coordsDeformedAerofoilFalk(odb,path_xfoilFiles,primaryJobName,jobNumber,index,chord,i+1,EiD)
##            twistCur[i] = coordsDeformedAerofoilFalkCAMBER(odb,path_xfoilFiles,primaryJobName,jobNumber,index,chord,i+1,EiD,xCL, yCL)

pp = []
pI = []
pP = []
#cP = []
#cI = []

Elem_NO_TIP='A_ALL_ELEMENTS_NO_TIP_CORNER'
setElem_NO_TIP= odb.rootAssembly.elementSets[Elem_NO_TIP]

pField = iterations.fieldOutputs['PDLOAD']
mStress = iterations.fieldOutputs['MSTRS']
mStressTsaiW= iterations.fieldOutputs['TSAIW']
mStrain= iterations.fieldOutputs['MSTRN']

MAXSTR_NO_TIP= mStress.getSubset(region=setElem_NO_TIP)
TSAIW_NO_TIP= mStressTsaiW.getSubset(region=setElem_NO_TIP)
STRAIN_NO_TIP= mStrain.getSubset(region=setElem_NO_TIP)

mStr_NO_TIP=0
mTsaiW_NO_TIP=0
meps_NO_TIP=0
for it in xrange(1,len(MAXSTR_NO_TIP.values)):
    if MAXSTR_NO_TIP.values[it].data > mStr_NO_TIP:
        mStr_NO_TIP=MAXSTR_NO_TIP.values[it].data
    if TSAIW_NO_TIP.values[it].data > mTsaiW_NO_TIP:
        mTsaiW_NO_TIP=TSAIW_NO_TIP.values[it].data
    if STRAIN_NO_TIP.values[it].data > meps_NO_TIP:
        meps_NO_TIP=STRAIN_NO_TIP.values[it].data





# Reaction Force auslesen

RFs = iterations.fieldOutputs['RF']
BCs = odb.rootAssembly.nodeSets['ALL_BC']
RF_BC = RFs.getSubset(region=BCs)

ReferencePoint_BC = odb.rootAssembly.nodeSets['REFERENCEPOINT']
RF_ReferencePoint_BC = RFs.getSubset(region=ReferencePoint_BC)
# exit()

RF_x_all=0
RF_y_all=0
RF_z_all=0
    
for j in RF_BC.values:
    RF_x=j.data[0]
    RF_y=j.data[1]
    RF_x_all=RF_x_all+RF_x
    RF_y_all=RF_y_all+RF_y

#Total RF
TotalRF_ReferencePoint = math.sqrt( math.pow(RF_ReferencePoint_BC.values[0].data[0], 2) + math.pow(RF_ReferencePoint_BC.values[0].data[1], 2) + math.pow(RF_ReferencePoint_BC.values[0].data[2], 2))


mStr=0
mTsaiW=0
meps=0
for it in xrange(1,len(mStress.values)):
    if mStress.values[it].data > mStr:
        mStr=mStress.values[it].data
    if mStressTsaiW.values[it].data > mTsaiW:
        mTsaiW=mStressTsaiW.values[it].data
    if mStrain.values[it].data > meps:
        meps=mStrain.values[it].data

pValues = pField.values
for node in pValues:
        pp.append([node.elementLabel,node.data])
        pI.append(node.elementLabel)
        pP.append(node.data)
#cField = odb.steps[stepName].frames[-1].fieldOutputs['COORD']
#cValues = cField.values
#for node in cValues:
#	cI.append(node.elementLabel)
#	cP.append(node.data[2])

Fx = []
Fz = []
Mx = []
Mz = []

Mxtot=0
FzTot=0

fo_fail = open(folder_path+'/rootBendingMoment/Failure-Job-'+str(iSIM) +'_'+str(elements_inWeb) +'_'+str(jobNumber)+'-' + str(index)+'.txt','a')
fo_fail.write(str(mStr_NO_TIP) + '\t' +str(mTsaiW_NO_TIP)+'\t' + str(meps_NO_TIP)+ '\t' + str(mStr)+ '\t' +str(mTsaiW)+ '\t' +str(meps)+'\n')   
fo_fail.close()

##        fo = open(folder_path+'/rootBendingMoment/Output-Job-'+str(jobNumber) +'-'+str(int(1000*PositionFrontSpar)) + '-' + str(1000*int(AngleFrontSpar)) + '-' + str(int(1000*PositionRearSpar)) + '-' + str(int(1000*AngleRearSpar)) +'-'+str(int(1000*thickness_bucklingSpar[digge]))+'-'+str(int(FaserWinkel[faserw]))+ '-' + str(index)+'.txt','a')
fo = open(folder_path+'/rootBendingMoment/Output-Job-'+str(iSIM) +'_'+str(elements_inWeb) +'_'+str(jobNumber)+'-' + str(index)+'.txt','a')
fo.write(str(iterations.frameId) + '\t' +str(iterations.frameValue)+'\t' + str(t_step_readout)+ '\t' +str(twist[-1])+ '\t' +str(twistCur[-1])+'\t' + str(TwistDifference)+'\t' + str(TotalRF_ReferencePoint) + '\t' + str(RF_y_all)+ '\t' +str(FzTot)+ '\t' +str(Mxtot)+ '\t' + str(mStr)+ '\t' +str(mTsaiW)+ '\t' +str(meps)+'\n')   
fo.close()

##f = open('D:/runkelf/New3Bending3/rootBendingMoment/pODB.txt','a')
##f.write('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\t' +'Fx = ' + str(FxTot)+ '\t'+ 'Fz = ' + str(FzTot)+ '\n')
###f.write('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\t' +'Fx = ' + str(FxPerp)+ '\t'+ 'Fz = ' + str(FzPerp)+ '\n')
###f.write('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\t' +'cD = ' + str(cdP)+ '\t'+ 'cL = ' + str(clP)+ '\n')
##f.write('Job = ' + str(jobNumber)+ '\t' + 'Step = ' + str(index) + '\t' +'Mx = ' + str(Mxtot)+ '\t'+ 'Mz = ' + str(Mztot)+ '\n')
##f.close()

odb.close()
