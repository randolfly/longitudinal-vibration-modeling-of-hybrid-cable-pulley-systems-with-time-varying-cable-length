# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2023 replay file
# Internal Version: 2022_09_29-02.11.55 183150
# Run by randolf on Thu May 15 14:38:59 2025
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=290.270629882812, 
    height=200.41667175293)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
#: Executing "onCaeStartup()" in the site directory ...
openMdb(
    pathName='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/osc_beam.cae')
#: The model database "D:\doc\randolf.lib\phd\papers\04-绳索动力学模型论文\数值仿真\julia\longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length\appendix\comment1\abaqus\fix_rod\osc_beam.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
#--- Recover file: 'osc_beam.rec' ---
# -*- coding: mbcs -*- 
from part import *
from material import *
from section import *
from optimization import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', 
    frequency=1, name='H-Output-3', variables=('IRF1', 'IRF2', 'IRF3', 'IRM1', 
    'IRM2', 'IRM3'))
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
mdb.jobs['Job-1']._Message(STARTED, {'phase': BATCHPRE_PHASE, 
    'clientHost': 'randolf', 'handle': 0, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    'message': 'OUTPUT REQUEST IRF1 IS NOT AVAILABLE FOR THIS TYPE OF ANALYSIS', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    'message': 'OUTPUT REQUEST IRF2 IS NOT AVAILABLE FOR THIS TYPE OF ANALYSIS', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    'message': 'OUTPUT REQUEST IRF3 IS NOT AVAILABLE FOR THIS TYPE OF ANALYSIS', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    'message': 'OUTPUT REQUEST IRM1 IS NOT AVAILABLE FOR THIS TYPE OF ANALYSIS', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    'message': 'OUTPUT REQUEST IRM2 IS NOT AVAILABLE FOR THIS TYPE OF ANALYSIS', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    'message': 'OUTPUT REQUEST IRM3 IS NOT AVAILABLE FOR THIS TYPE OF ANALYSIS', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(ODB_FILE, {'phase': BATCHPRE_PHASE, 
    'file': 'D:\\doc\\randolf.lib\\phd\\papers\\04-\xc9\xfe\xcb\xf7\xb6\xaf\xc1\xa6\xd1\xa7\xc4\xa3\xd0\xcd\xc2\xdb\xce\xc4\\\xca\xfd\xd6\xb5\xb7\xc2\xd5\xe6\\osc_test\\abaqus\\fix_rod_with_endpoint_mass\\Job-1.odb', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(COMPLETED, {'phase': BATCHPRE_PHASE, 
    'message': 'Analysis phase complete', 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STARTED, {'phase': STANDARD_PHASE, 
    'clientHost': 'randolf', 'handle': 28704, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STEP, {'phase': STANDARD_PHASE, 'stepId': 1, 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(MEMORY_ESTIMATE, {'phase': STANDARD_PHASE, 
    'jobName': 'Job-1', 'memory': 31.0})
mdb.jobs['Job-1']._Message(PHYSICAL_MEMORY, {'phase': STANDARD_PHASE, 
    'physical_memory': 24254.0, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(MINIMUM_MEMORY, {'minimum_memory': 22.0, 
    'phase': STANDARD_PHASE, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 0, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(MEMORY_ESTIMATE, {'phase': STANDARD_PHASE, 
    'jobName': 'Job-1', 'memory': 31.0})
mdb.jobs['Job-1']._Message(PHYSICAL_MEMORY, {'phase': STANDARD_PHASE, 
    'physical_memory': 24254.0, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(MINIMUM_MEMORY, {'minimum_memory': 22.0, 
    'phase': STANDARD_PHASE, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1e-05, 'attempts': 1, 
    'timeIncrement': 1e-05, 'increment': 1, 'stepTime': 1e-05, 'step': 1, 
    'jobName': 'Job-1', 'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 
    'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 1, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2e-05, 'attempts': 1, 
    'timeIncrement': 1e-05, 'increment': 2, 'stepTime': 2e-05, 'step': 1, 
    'jobName': 'Job-1', 'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 
    'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 2, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 3e-05, 'attempts': 1, 
    'timeIncrement': 1e-05, 'increment': 3, 'stepTime': 3e-05, 'step': 1, 
    'jobName': 'Job-1', 'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 
    'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 3, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 4.25e-05, 'attempts': 1, 
    'timeIncrement': 1.25e-05, 'increment': 4, 'stepTime': 4.25e-05, 'step': 1, 
    'jobName': 'Job-1', 'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 
    'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 4, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 5.8125e-05, 'attempts': 1, 
    'timeIncrement': 1.5625e-05, 'increment': 5, 'stepTime': 5.8125e-05, 
    'step': 1, 'jobName': 'Job-1', 'severe': 0, 'iterations': 1, 
    'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 5, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 7.765625e-05, 'attempts': 1, 
    'timeIncrement': 1.953125e-05, 'increment': 6, 'stepTime': 7.765625e-05, 
    'step': 1, 'jobName': 'Job-1', 'severe': 0, 'iterations': 1, 
    'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 6, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0001020703125, 
    'attempts': 1, 'timeIncrement': 2.44140625e-05, 'increment': 7, 
    'stepTime': 0.0001020703125, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 7, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000132587890625, 
    'attempts': 1, 'timeIncrement': 3.0517578125e-05, 'increment': 8, 
    'stepTime': 0.000132587890625, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 8, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00017073486328125, 
    'attempts': 1, 'timeIncrement': 3.814697265625e-05, 'increment': 9, 
    'stepTime': 0.00017073486328125, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 9, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000218418579101563, 
    'attempts': 1, 'timeIncrement': 4.76837158203125e-05, 'increment': 10, 
    'stepTime': 0.000218418579101563, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 10, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000278023223876953, 
    'attempts': 1, 'timeIncrement': 5.96046447753906e-05, 'increment': 11, 
    'stepTime': 0.000278023223876953, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 11, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000352529029846191, 
    'attempts': 1, 'timeIncrement': 7.45058059692383e-05, 'increment': 12, 
    'stepTime': 0.000352529029846191, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 12, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000445661287307739, 
    'attempts': 1, 'timeIncrement': 9.31322574615479e-05, 'increment': 13, 
    'stepTime': 0.000445661287307739, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 13, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000562076609134674, 
    'attempts': 1, 'timeIncrement': 0.000116415321826935, 'increment': 14, 
    'stepTime': 0.000562076609134674, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 14, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000707595761418343, 
    'attempts': 1, 'timeIncrement': 0.000145519152283669, 'increment': 15, 
    'stepTime': 0.000707595761418343, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 15, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.000889494701772928, 
    'attempts': 1, 'timeIncrement': 0.000181898940354586, 'increment': 16, 
    'stepTime': 0.000889494701772928, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 16, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00111686837721616, 
    'attempts': 1, 'timeIncrement': 0.000227373675443232, 'increment': 17, 
    'stepTime': 0.00111686837721616, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 17, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0014010854715202, 
    'attempts': 1, 'timeIncrement': 0.00028421709430404, 'increment': 18, 
    'stepTime': 0.0014010854715202, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 18, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00175635683940025, 
    'attempts': 1, 'timeIncrement': 0.00035527136788005, 'increment': 19, 
    'stepTime': 0.00175635683940025, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 19, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00220044604925031, 
    'attempts': 1, 'timeIncrement': 0.000444089209850063, 'increment': 20, 
    'stepTime': 0.00220044604925031, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 20, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00275555756156289, 
    'attempts': 1, 'timeIncrement': 0.000555111512312578, 'increment': 21, 
    'stepTime': 0.00275555756156289, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 21, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00344944695195361, 
    'attempts': 1, 'timeIncrement': 0.000693889390390723, 'increment': 22, 
    'stepTime': 0.00344944695195361, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 22, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00431680868994202, 
    'attempts': 1, 'timeIncrement': 0.000867361737988404, 'increment': 23, 
    'stepTime': 0.00431680868994202, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 23, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.00540101086242752, 
    'attempts': 1, 'timeIncrement': 0.0010842021724855, 'increment': 24, 
    'stepTime': 0.00540101086242752, 'step': 1, 'jobName': 'Job-1', 
    'severe': 0, 'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 24, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0067562635780344, 
    'attempts': 1, 'timeIncrement': 0.00135525271560688, 'increment': 25, 
    'stepTime': 0.0067562635780344, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 25, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.008450329472543, 
    'attempts': 1, 'timeIncrement': 0.0016940658945086, 'increment': 26, 
    'stepTime': 0.008450329472543, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 26, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0105679118406788, 
    'attempts': 1, 'timeIncrement': 0.00211758236813575, 'increment': 27, 
    'stepTime': 0.0105679118406788, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 27, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0132148898008484, 
    'attempts': 1, 'timeIncrement': 0.00264697796016969, 'increment': 28, 
    'stepTime': 0.0132148898008484, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 28, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0165236122510606, 
    'attempts': 1, 'timeIncrement': 0.00330872245021211, 'increment': 29, 
    'stepTime': 0.0165236122510606, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 29, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0206595153138257, 
    'attempts': 1, 'timeIncrement': 0.00413590306276514, 'increment': 30, 
    'stepTime': 0.0206595153138257, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 30, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0258293941422821, 
    'attempts': 1, 'timeIncrement': 0.00516987882845642, 'increment': 31, 
    'stepTime': 0.0258293941422821, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 31, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0322917426778526, 
    'attempts': 1, 'timeIncrement': 0.00646234853557053, 'increment': 32, 
    'stepTime': 0.0322917426778526, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 32, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0403696783473158, 
    'attempts': 1, 'timeIncrement': 0.00807793566946316, 'increment': 33, 
    'stepTime': 0.0403696783473158, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 33, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0504670979341448, 
    'attempts': 1, 'timeIncrement': 0.0100974195868289, 'increment': 34, 
    'stepTime': 0.0504670979341448, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 34, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0630888724176809, 
    'attempts': 1, 'timeIncrement': 0.0126217744835362, 'increment': 35, 
    'stepTime': 0.0630888724176809, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 35, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0788660905221012, 
    'attempts': 1, 'timeIncrement': 0.0157772181044202, 'increment': 36, 
    'stepTime': 0.0788660905221012, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 36, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.0985876131526265, 
    'attempts': 1, 'timeIncrement': 0.0197215226305253, 'increment': 37, 
    'stepTime': 0.0985876131526265, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 37, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.123239516440783, 
    'attempts': 1, 'timeIncrement': 0.0246519032881566, 'increment': 38, 
    'stepTime': 0.123239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 38, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.153239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 39, 
    'stepTime': 0.153239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 39, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.183239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 40, 
    'stepTime': 0.183239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 40, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.213239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 41, 
    'stepTime': 0.213239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 41, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.243239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 42, 
    'stepTime': 0.243239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 42, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.273239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 43, 
    'stepTime': 0.273239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 43, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.303239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 44, 
    'stepTime': 0.303239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 44, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.333239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 45, 
    'stepTime': 0.333239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 45, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.363239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 46, 
    'stepTime': 0.363239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 46, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.393239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 47, 
    'stepTime': 0.393239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 47, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.423239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 48, 
    'stepTime': 0.423239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 48, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.453239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 49, 
    'stepTime': 0.453239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 49, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.483239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 50, 
    'stepTime': 0.483239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 50, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.513239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 51, 
    'stepTime': 0.513239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 51, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.543239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 52, 
    'stepTime': 0.543239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 52, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.573239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 53, 
    'stepTime': 0.573239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 53, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.603239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 54, 
    'stepTime': 0.603239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 54, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.633239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 55, 
    'stepTime': 0.633239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 55, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.663239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 56, 
    'stepTime': 0.663239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 56, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.693239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 57, 
    'stepTime': 0.693239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 57, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.723239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 58, 
    'stepTime': 0.723239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 58, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.753239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 59, 
    'stepTime': 0.753239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 59, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.783239516440783, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 60, 
    'stepTime': 0.783239516440783, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 60, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.813239516440784, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 61, 
    'stepTime': 0.813239516440784, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 61, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.843239516440784, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 62, 
    'stepTime': 0.843239516440784, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 62, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.873239516440784, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 63, 
    'stepTime': 0.873239516440784, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 63, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.903239516440784, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 64, 
    'stepTime': 0.903239516440784, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 64, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.933239516440784, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 65, 
    'stepTime': 0.933239516440784, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 65, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.963239516440784, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 66, 
    'stepTime': 0.963239516440784, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 66, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 0.993239516440784, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 67, 
    'stepTime': 0.993239516440784, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 67, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.02323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 68, 
    'stepTime': 1.02323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 68, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.05323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 69, 
    'stepTime': 1.05323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 69, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.08323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 70, 
    'stepTime': 1.08323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 70, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.11323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 71, 
    'stepTime': 1.11323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 71, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.14323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 72, 
    'stepTime': 1.14323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 72, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.17323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 73, 
    'stepTime': 1.17323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 73, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.20323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 74, 
    'stepTime': 1.20323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 74, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.23323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 75, 
    'stepTime': 1.23323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 75, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.26323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 76, 
    'stepTime': 1.26323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 76, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.29323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 77, 
    'stepTime': 1.29323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 77, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.32323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 78, 
    'stepTime': 1.32323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 78, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.35323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 79, 
    'stepTime': 1.35323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 79, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.38323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 80, 
    'stepTime': 1.38323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 80, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.41323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 81, 
    'stepTime': 1.41323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 81, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.44323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 82, 
    'stepTime': 1.44323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 82, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.47323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 83, 
    'stepTime': 1.47323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 83, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.50323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 84, 
    'stepTime': 1.50323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 84, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.53323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 85, 
    'stepTime': 1.53323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 85, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.56323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 86, 
    'stepTime': 1.56323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 86, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.59323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 87, 
    'stepTime': 1.59323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 87, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.62323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 88, 
    'stepTime': 1.62323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 88, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.65323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 89, 
    'stepTime': 1.65323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 89, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.68323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 90, 
    'stepTime': 1.68323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 90, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.71323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 91, 
    'stepTime': 1.71323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 91, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.74323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 92, 
    'stepTime': 1.74323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 92, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.77323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 93, 
    'stepTime': 1.77323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 93, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.80323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 94, 
    'stepTime': 1.80323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 94, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.83323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 95, 
    'stepTime': 1.83323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 95, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.86323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 96, 
    'stepTime': 1.86323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 96, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.89323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 97, 
    'stepTime': 1.89323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 97, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.92323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 98, 
    'stepTime': 1.92323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 98, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.95323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 99, 
    'stepTime': 1.95323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 99, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 1.98323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 100, 
    'stepTime': 1.98323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 100, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.01323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 101, 
    'stepTime': 2.01323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 101, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.04323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 102, 
    'stepTime': 2.04323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 102, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.07323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 103, 
    'stepTime': 2.07323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 103, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.10323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 104, 
    'stepTime': 2.10323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 104, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.13323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 105, 
    'stepTime': 2.13323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 105, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.16323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 106, 
    'stepTime': 2.16323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 106, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.19323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 107, 
    'stepTime': 2.19323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 107, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.22323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 108, 
    'stepTime': 2.22323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 108, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.25323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 109, 
    'stepTime': 2.25323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 109, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.28323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 110, 
    'stepTime': 2.28323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 110, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.31323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 111, 
    'stepTime': 2.31323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 111, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.34323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 112, 
    'stepTime': 2.34323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 112, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.37323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 113, 
    'stepTime': 2.37323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 113, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.40323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 114, 
    'stepTime': 2.40323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 114, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.43323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 115, 
    'stepTime': 2.43323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 115, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.46323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 116, 
    'stepTime': 2.46323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 116, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.49323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 117, 
    'stepTime': 2.49323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 117, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.52323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 118, 
    'stepTime': 2.52323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 118, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.55323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 119, 
    'stepTime': 2.55323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 119, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.58323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 120, 
    'stepTime': 2.58323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 120, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.61323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 121, 
    'stepTime': 2.61323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 121, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.64323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 122, 
    'stepTime': 2.64323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 122, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.67323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 123, 
    'stepTime': 2.67323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 123, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.70323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 124, 
    'stepTime': 2.70323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 124, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.73323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 125, 
    'stepTime': 2.73323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 125, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.76323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 126, 
    'stepTime': 2.76323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 126, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.79323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 127, 
    'stepTime': 2.79323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 127, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.82323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 128, 
    'stepTime': 2.82323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 128, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.85323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 129, 
    'stepTime': 2.85323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 129, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.88323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 130, 
    'stepTime': 2.88323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 130, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.91323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 131, 
    'stepTime': 2.91323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 131, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.94323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 132, 
    'stepTime': 2.94323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 132, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 2.97323951644078, 
    'attempts': 1, 'timeIncrement': 0.03, 'increment': 133, 
    'stepTime': 2.97323951644078, 'step': 1, 'jobName': 'Job-1', 'severe': 0, 
    'iterations': 1, 'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 133, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(WARNING, {'phase': STANDARD_PHASE, 
    'message': 'There is zero MOMENT everywhere in the model based on the default criterion. please check the value of the average MOMENT during the current iteration to verify that the MOMENT is small enough to be treated as zero. if not, please use the solution controls to reset the criterion for zero MOMENT.', 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    'frame': 134, 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(STATUS, {'totalTime': 3.0, 'attempts': 1, 
    'timeIncrement': 0.0267604835592219, 'increment': 134, 'stepTime': 3.0, 
    'step': 1, 'jobName': 'Job-1', 'severe': 0, 'iterations': 1, 
    'phase': STANDARD_PHASE, 'equilibrium': 1})
mdb.jobs['Job-1']._Message(END_STEP, {'phase': STANDARD_PHASE, 'stepId': 1, 
    'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(COMPLETED, {'phase': STANDARD_PHASE, 
    'message': 'Analysis phase complete', 'jobName': 'Job-1'})
mdb.jobs['Job-1']._Message(JOB_COMPLETED, {'time': 'Mon Apr  7 16:36:24 2025', 
    'jobName': 'Job-1'})
#--- End of Recover file ------
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.mdbData['Model-1'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.64956, 
    farPlane=4.66751, width=1.44794, height=0.629502, viewOffsetX=0.00961619, 
    viewOffsetY=0.0123532)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
o1 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          3
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].setValues(displayedObject=None)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb'])
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'V', 'A', 'RF', 'CF', 
    'CSTRESS', 'CDISP'))
mdb.models['Model-1'].historyOutputRequests['H-Output-2'].setValues(variables=(
    'E11', 'E22', 'E33', 'E12', 'E13', 'E23', 'EP', 'U1', 'U2', 'U3', 'UR1', 
    'UR2', 'UR3', 'COOR1', 'COOR2', 'COOR3'))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: The job input file "Job-1.inp" has been submitted for analysis.
o3 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
#: Model: D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          3
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
#: Job Job-1: Analysis Input File Processor completed successfully.
#: Job Job-1: Abaqus/Standard completed successfully.
#: Job Job-1 completed successfully. 
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb'])
o3 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
#: Model: D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          3
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ), 
    elementSets=("PART-1-1._I1", ))
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
#: Warning: Requested operation will result in the creation of a very large number of xyDataObjects. Performance can be affected. Please reduce the number of specified entities using Display Group operations before re-performing this operation.
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E12'), )), ), 
    elementSets=("PART-1-1._I1", ))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
#: Warning: Requested operation will result in the creation of a very large number of xyDataObjects. Performance can be affected. Please reduce the number of specified entities using Display Group operations before re-performing this operation.
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ), 
    elementSets=("PART-1-1._I2", ))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
#: Warning: Requested operation will result in the creation of a very large number of xyDataObjects. Performance can be affected. Please reduce the number of specified entities using Display Group operations before re-performing this operation.
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].setValues(displayedObject=None)
xyp = session.xyPlots['XYPlot-1']
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects['_E:E11 SP:1 PI: PART-1-1 E: 301 IP: 1']
c1 = session.Curve(xyData=xy1)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Artificial strain energy: ALLAE for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Contact constraint discontinuity work: ALLCCDW for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Contact constraint elastic energy: ALLCCE for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Contact constraint elastic normal energy: ALLCCEN for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Coordinates: COOR1 at Node 2 in NSET FORCE_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Coordinates: COOR2 at Node 2 in NSET FORCE_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Coordinates: COOR3 at Node 2 in NSET FORCE_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Creep dissipation energy: ALLCD for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Damage dissipation energy: ALLDMD for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Dynamic time integration energy: ALLDTI for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Electrostatic energy: ALLEE for Whole Model', steps=(
    'Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Energy lost to quiet boundaries: ALLQB for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='External work: ALLWK for Whole Model', steps=('Step-1', 
    ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Frictional dissipation: ALLFD for Whole Model', steps=(
    'Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Internal energy: ALLIE for Whole Model', steps=(
    'Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Plastic dissipation: ALLPD for Whole Model', steps=(
    'Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Rotational displacement: UR1 at Node 2 in NSET FORCE_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial displacement: U1 at Node 2 in NSET FORCE_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 2 in NSET FORCE_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial displacement: U3 at Node 2 in NSET FORCE_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Static dissipation (stabilization): ALLSD for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Strain energy: ALLSE for Whole Model', steps=('Step-1', 
    ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Total energy of the output set: ETOTAL for Whole Model', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Viscous dissipation: ALLVD for Whole Model', steps=(
    'Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.mdbData['Model-1'])
session.viewports['Viewport: 1'].odbDisplay.setFrame(step='Step-1')
o7 = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=o7)
x = session.xyPlots['XYPlot-1']
session.viewports['Viewport: 1'].setValues(displayedObject=x)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects['_temp_2']
c1 = session.Curve(xyData=xy1)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects['_temp_2']
c1 = session.Curve(xyData=xy1)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
del session.xyDataObjects['_temp_2']
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.models['Model-1'].historyOutputRequests['H-Output-2'].setValues(variables=(
    'E11', 'E22', 'E33', 'E12', 'E13', 'E23', 'EP', 'VE11', 'VE22', 'VE33', 
    'VE12', 'VE13', 'VE23', 'PE11', 'PE22', 'PE33', 'PE12', 'PE13', 'PE23', 
    'PEP', 'VEEQ', 'PEEQ', 'PEEQT', 'PEMAG', 'PEQC1', 'PEQC2', 'PEQC3', 
    'PEQC4', 'EE11', 'EE22', 'EE33', 'EE12', 'EE13', 'EE23', 'EEP', 'IE11', 
    'IE22', 'IE33', 'IE12', 'IE13', 'IE23', 'IEP', 'THE11', 'THE22', 'THE33', 
    'THE12', 'THE13', 'THE23', 'THEP', 'NE11', 'NE22', 'NE33', 'NE12', 'NE13', 
    'NE23', 'NEP', 'LE11', 'LE22', 'LE33', 'LE12', 'LE13', 'LE23', 'LEP', 
    'TE11', 'TE22', 'TE33', 'TE12', 'TE13', 'TE23', 'TEEQ', 'TEVOL', 'EEQUT', 
    'ER11', 'ER22', 'ER33', 'ER12', 'ER13', 'ER23', 'ERP', 'SE1', 'SE2', 'SE3', 
    'SK1', 'SK2', 'SK3', 'BICURV', 'SPE1', 'SPE2', 'SPE3', 'SPE4', 'SEPE1', 
    'SEPE2', 'SEPE3', 'SEPE4', 'SEE1', 'SKE1', 'SKE2', 'SKE3', 'SEP1', 'SKP1', 
    'SKP2', 'SKP3', 'SALPHA1', 'SALPHA2', 'SALPHA3', 'SALPHA4', 'NBEEQ', 
    'NBPEEQ', 'GKEEQ', 'GKPEEQ', 'U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3', 
    'COOR1', 'COOR2', 'COOR3'))
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'PE', 'PEEQ', 'PEMAG', 'EE', 'LE', 'U', 'V', 'A', 'RF', 'CF', 
    'CSTRESS', 'CDISP'))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: The job input file "Job-1.inp" has been submitted for analysis.
o3 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
#: Model: D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          3
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
#: Job Job-1: Analysis Input File Processor completed successfully.
#: Job Job-1: Abaqus/Standard completed successfully.
#: Job Job-1 completed successfully. 
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb'])
#: Warning: The output database 'D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb' disk file has changed.
#: 
#: The current plot operation has been canceled, re-open the file to view the results
o3 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
#: Model: D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          3
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
o1 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ), 
    elementPick=(('PART-1-1', 1, ('[#0:7 #100000 ]', )), ), )
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects['_E:E11 SP:1 PI: PART-1-1 E: 245 IP: 1']
c1 = session.Curve(xyData=xy1)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    SYMBOLS_ON_DEF, ))
#* The basic plot option 'render beam profiles' is supported only in deformed 
#* or undeformed or contour plot modes.
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.88166, 
    farPlane=2.15814, width=1.00556, height=0.437175, viewOffsetX=0.000596225, 
    viewOffsetY=0.00730695)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ), 
    elementPick=(('PART-1-1', 5, ('[#0:2 #4 #0 #80000000 #0:3 #200 #0:2', 
    ' #10000 #0:3 #40000 ]', )), ), )
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['Part-1']
p.DatumPointByCoordinate(coords=(0.1, 0.0, 0.0))
p = mdb.models['Model-1'].parts['Part-1']
p.DatumPointByCoordinate(coords=(0.2, 0.0, 0.0))
p = mdb.models['Model-1'].parts['Part-1']
p.DatumPointByCoordinate(coords=(0.3, 0.0, 0.0))
p = mdb.models['Model-1'].parts['Part-1']
p.DatumPointByCoordinate(coords=(0.4, 0.0, 0.0))
p = mdb.models['Model-1'].parts['Part-1']
p.DatumPointByCoordinate(coords=(0.5, 0.0, 0.0))
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF, adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
p = mdb.models['Model-1'].parts['Part-1']
n = p.nodes
nodes = n.getSequenceFromMask(mask=('[#0:3 #80 #0:2 #400 #0:2 #8000 #0:2', 
    ' #80000 #0:2 #800000 ]', ), )
p.Set(nodes=nodes, name='plot_set')
#: The set 'plot_set' has been created (5 nodes).
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
regionDef=mdb.models['Model-1'].rootAssembly.allInstances['Part-1-1'].sets['plot_set']
mdb.models['Model-1'].FieldOutputRequest(name='F-Output-2', 
    createStepName='Step-1', variables=('E', 'VE', 'PE', 'VEEQ', 'PEEQ', 
    'PEEQT', 'PEEQMAX', 'PEMAG', 'PEQC', 'EE', 'IE', 'THE', 'NE', 'LE', 'TE', 
    'TEEQ', 'TEVOL', 'EEQUT', 'ER', 'SE', 'SPE', 'SEPE', 'SEE', 'SEP', 
    'SALPHA', 'NBEEQ', 'NBPEEQ', 'GKEEQ', 'GKPEEQ', 'U', 'UT', 'UR', 'V', 'VT', 
    'VR', 'A', 'AT', 'AR', 'RBANG', 'RBROT', 'RF', 'RT', 'RM', 'CF', 'SF', 
    'SQEQ', 'TF', 'VF', 'ESF1', 'NFORC', 'NFORCSO', 'RBFOR', 'BF', 'CORIOMAG', 
    'ROTAMAG', 'CENTMAG', 'CENTRIFMAG', 'GRAV', 'P', 'HP', 'TRSHR', 'TRNOR', 
    'TRVEC'), frequency=1, region=regionDef, sectionPoints=DEFAULT, 
    rebar=EXCLUDE)
regionDef=mdb.models['Model-1'].rootAssembly.allInstances['Part-1-1'].sets['plot_set']
mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-4', 
    createStepName='Step-1', variables=('S11', 'S22', 'S33', 'S12', 'S13', 
    'S23', 'SP', 'TRESC', 'PRESS', 'INV3', 'MISES', 'TSHR13', 'TSHR23', 
    'CTSHR13', 'CTSHR23', 'ALPHA11', 'ALPHA22', 'ALPHA33', 'ALPHA12', 
    'ALPHA13', 'ALPHA23', 'ALPHAP', 'TRIAX', 'LODE', 'VS11', 'VS22', 'VS33', 
    'VS12', 'VS13', 'VS23', 'PS11', 'PS22', 'PS33', 'PS12', 'PS13', 'PS23', 
    'CS11', 'ALPHA1', 'ALPHA2', 'ALPHA3', 'ALPHA4', 'ALPHA5', 'ALPHA6', 
    'ALPHA7', 'ALPHA8', 'ALPHA9', 'ALPHA10', 'SSAVG1', 'SSAVG2', 'SSAVG3', 
    'SSAVG4', 'SSAVG5', 'SSAVG6', 'SEQUT', 'YIELDPOT', 'NBSEQ', 'GKSEQ', 'E11', 
    'E22', 'E33', 'E12', 'E13', 'E23', 'EP', 'VE11', 'VE22', 'VE33', 'VE12', 
    'VE13', 'VE23', 'PE11', 'PE22', 'PE33', 'PE12', 'PE13', 'PE23', 'PEP', 
    'VEEQ', 'PEEQ', 'PEEQT', 'PEMAG', 'PEQC1', 'PEQC2', 'PEQC3', 'PEQC4', 
    'EE11', 'EE22', 'EE33', 'EE12', 'EE13', 'EE23', 'EEP', 'IE11', 'IE22', 
    'IE33', 'IE12', 'IE13', 'IE23', 'IEP', 'THE11', 'THE22', 'THE33', 'THE12', 
    'THE13', 'THE23', 'THEP', 'NE11', 'NE22', 'NE33', 'NE12', 'NE13', 'NE23', 
    'NEP', 'LE11', 'LE22', 'LE33', 'LE12', 'LE13', 'LE23', 'LEP', 'TE11', 
    'TE22', 'TE33', 'TE12', 'TE13', 'TE23', 'TEEQ', 'TEVOL', 'EEQUT', 'ER11', 
    'ER22', 'ER33', 'ER12', 'ER13', 'ER23', 'ERP', 'SE1', 'SE2', 'SE3', 'SK1', 
    'SK2', 'SK3', 'BICURV', 'SPE1', 'SPE2', 'SPE3', 'SPE4', 'SEPE1', 'SEPE2', 
    'SEPE3', 'SEPE4', 'SEE1', 'SKE1', 'SKE2', 'SKE3', 'SEP1', 'SKP1', 'SKP2', 
    'SKP3', 'SALPHA1', 'SALPHA2', 'SALPHA3', 'SALPHA4', 'NBEEQ', 'NBPEEQ', 
    'GKEEQ', 'GKPEEQ', 'U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3', 'UT', 'UR', 'V1', 
    'V2', 'V3', 'VR1', 'VR2', 'VR3', 'VT', 'VR', 'A1', 'A2', 'A3', 'AR1', 
    'AR2', 'AR3', 'AT', 'AR', 'GU', 'GV', 'GA', 'WARP', 'RBANG', 'RBROT', 
    'IRA1', 'IRA2', 'IRA3', 'IRAR1', 'IRAR2', 'IRAR3', 'RF1', 'RF2', 'RF3', 
    'RM1', 'RM2', 'RM3', 'RT', 'RM', 'RWM', 'CF1', 'CF2', 'CF3', 'CM1', 'CM2', 
    'CM3', 'CW', 'SF1', 'SF2', 'SF3', 'SM1', 'SM2', 'SM3', 'BIMOM', 'SQEQ', 
    'TF1', 'TF2', 'TF3', 'TM1', 'TM2', 'TM3', 'VF1', 'VF2', 'VF3', 'VM1', 
    'VM2', 'VM3', 'ESF1', 'NFORC', 'NFORCSO', 'RBFOR', 'IRF1', 'IRF2', 'IRF3', 
    'IRM1', 'IRM2', 'IRM3'), frequency=1, region=regionDef, 
    sectionPoints=DEFAULT, rebar=EXCLUDE)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: The job input file "Job-1.inp" has been submitted for analysis.
#: Job Job-1: Analysis Input File Processor completed successfully.
#: Job Job-1: Abaqus/Standard completed successfully.
#: Job Job-1 completed successfully. 
o3 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
#: Model: D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          4
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
o1 = session.openOdb(
    name='D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Point loads: CF1 at Node 104 in NSET PLOT_SET', steps=(
    'Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Reaction force: RF1 at Node 104 in NSET PLOT_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial acceleration: A1 at Node 104 in NSET PLOT_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial displacement: U1 at Node 404 in NSET PLOT_SET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ), 
    elementLabels=(('PART-1-1', ('101', '201', )), ))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ), 
    elementLabels=(('PART-1-1', ('101', '201', '301', '401', '501', )), ))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('E', INTEGRATION_POINT, ((COMPONENT, 'E11'), )), ), 
    elementLabels=(('PART-1-1', ('100', '200', '300', '400', '500', )), ))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
x = session.xyPlots['XYPlot-1']
session.viewports['Viewport: 1'].setValues(displayedObject=x)
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
del session.xyDataObjects['_E:E11 SP:1 PI: PART-1-1 E: 100 IP: 1']
del session.xyDataObjects['_E:E11 SP:1 PI: PART-1-1 E: 200 IP: 1']
del session.xyDataObjects['_E:E11 SP:1 PI: PART-1-1 E: 300 IP: 1']
del session.xyDataObjects['_E:E11 SP:1 PI: PART-1-1 E: 400 IP: 1']
del session.xyDataObjects['_E:E11 SP:1 PI: PART-1-1 E: 500 IP: 1']
del session.xyDataObjects['_E:E11 SP:5 PI: PART-1-1 E: 100 IP: 1']
del session.xyDataObjects['_E:E11 SP:5 PI: PART-1-1 E: 200 IP: 1']
del session.xyDataObjects['_E:E11 SP:5 PI: PART-1-1 E: 300 IP: 1']
del session.xyDataObjects['_E:E11 SP:5 PI: PART-1-1 E: 400 IP: 1']
del session.xyDataObjects['_E:E11 SP:5 PI: PART-1-1 E: 500 IP: 1']
del session.xyDataObjects['_E:E11 SP:21 PI: PART-1-1 E: 100 IP: 1']
del session.xyDataObjects['_E:E11 SP:21 PI: PART-1-1 E: 200 IP: 1']
del session.xyDataObjects['_E:E11 SP:21 PI: PART-1-1 E: 300 IP: 1']
del session.xyDataObjects['_E:E11 SP:21 PI: PART-1-1 E: 400 IP: 1']
del session.xyDataObjects['_E:E11 SP:21 PI: PART-1-1 E: 500 IP: 1']
del session.xyDataObjects['_E:E11 SP:25 PI: PART-1-1 E: 100 IP: 1']
del session.xyDataObjects['_E:E11 SP:25 PI: PART-1-1 E: 200 IP: 1']
del session.xyDataObjects['_E:E11 SP:25 PI: PART-1-1 E: 300 IP: 1']
del session.xyDataObjects['_E:E11 SP:25 PI: PART-1-1 E: 400 IP: 1']
del session.xyDataObjects['_E:E11 SP:25 PI: PART-1-1 E: 500 IP: 1']
odb = session.odbs['D:/doc/randolf.lib/phd/papers/04-绳索动力学模型论文/数值仿真/julia/longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length/appendix/comment1/abaqus/fix_rod/Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
    variable=(('EE', INTEGRATION_POINT, ((COMPONENT, 'EE11'), )), ), 
    elementLabels=(('PART-1-1', ('100', '200', '300', '400', '500', )), ))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
del session.xyDataObjects['_EE:EE11 SP:5 PI: PART-1-1 E: 100 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:5 PI: PART-1-1 E: 200 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:5 PI: PART-1-1 E: 300 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:5 PI: PART-1-1 E: 400 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:5 PI: PART-1-1 E: 500 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:21 PI: PART-1-1 E: 100 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:21 PI: PART-1-1 E: 200 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:21 PI: PART-1-1 E: 300 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:21 PI: PART-1-1 E: 400 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:21 PI: PART-1-1 E: 500 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:25 PI: PART-1-1 E: 100 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:25 PI: PART-1-1 E: 200 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:25 PI: PART-1-1 E: 300 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:25 PI: PART-1-1 E: 400 IP: 1']
del session.xyDataObjects['_EE:EE11 SP:25 PI: PART-1-1 E: 500 IP: 1']
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects['_EE:EE11 SP:1 PI: PART-1-1 E: 100 IP: 1']
c1 = session.Curve(xyData=xy1)
xy2 = session.xyDataObjects['_EE:EE11 SP:1 PI: PART-1-1 E: 200 IP: 1']
c2 = session.Curve(xyData=xy2)
xy3 = session.xyDataObjects['_EE:EE11 SP:1 PI: PART-1-1 E: 300 IP: 1']
c3 = session.Curve(xyData=xy3)
xy4 = session.xyDataObjects['_EE:EE11 SP:1 PI: PART-1-1 E: 400 IP: 1']
c4 = session.Curve(xyData=xy4)
xy5 = session.xyDataObjects['_EE:EE11 SP:1 PI: PART-1-1 E: 500 IP: 1']
c5 = session.Curve(xyData=xy5)
chart.setValues(curvesToPlot=(c1, c2, c3, c4, c5, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
mdb.save()
#: The model database has been saved to "D:\doc\randolf.lib\phd\papers\04-绳索动力学模型论文\数值仿真\julia\longitudinal-vibration-modeling-of-hybrid-cable-pulley-systems-with-time-varying-cable-length\appendix\comment1\abaqus\fix_rod\osc_beam.cae".
