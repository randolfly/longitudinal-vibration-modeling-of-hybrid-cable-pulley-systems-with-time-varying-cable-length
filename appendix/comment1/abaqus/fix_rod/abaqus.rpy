# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2023 replay file
# Internal Version: 2022_09_29-02.11.55 183150
# Run by randolf on Sun Nov 10 11:19:18 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=245.008895874023, 
    height=152.899337768555)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
#: Executing "onCaeStartup()" in the site directory ...
openMdb(pathName='D:/tmp/desktop/osc_test/abaqus/fix_rod/osc_beam.cae')
#: The model database "D:\tmp\desktop\osc_test\abaqus\fix_rod\osc_beam.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].setValues(displayedObject=None)
odb = session.mdbData['Model-1']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.51853, 
    farPlane=4.79855, width=1.89399, height=0.867366, viewOffsetX=0.00235913, 
    viewOffsetY=-0.0110244)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.6184, 
    farPlane=4.731, width=1.9691, height=0.901763, cameraPosition=(2.06596, 
    1.13837, 3.12374), cameraUpVector=(-0.268954, 0.785205, -0.557779), 
    cameraTarget=(0.524684, -0.0398475, 0.0219367), viewOffsetX=0.00245269, 
    viewOffsetY=-0.0114616)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
o1 = session.openOdb(name='D:/tmp/desktop/osc_test/abaqus/fix_rod/Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: D:/tmp/desktop/osc_test/abaqus/fix_rod/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       1
#: Number of Node Sets:          3
#: Number of Steps:              1
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.83158, 
    farPlane=2.20822, width=0.976271, height=0.447091, viewOffsetX=0.00164913, 
    viewOffsetY=0.110535)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.81904, 
    farPlane=2.22076, width=1.16735, height=0.534599, viewOffsetX=0.0138313, 
    viewOffsetY=0.121308)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=134 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=134 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=134 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=134 )
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=TIME_HISTORY)
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
session.animationOptions.setValues(frameRate=26)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.7878, 
    farPlane=2.35866, width=1.14731, height=0.525418, cameraPosition=(0.774655, 
    0.780728, 1.901), cameraUpVector=(-0.189745, 0.915043, -0.355941), 
    cameraTarget=(0.525957, 0.00944773, 0.0507895), viewOffsetX=0.0135938, 
    viewOffsetY=0.119225)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.74687, 
    farPlane=2.39959, width=1.53375, height=0.702394, viewOffsetX=0.010494, 
    viewOffsetY=0.166013)
session.viewports['Viewport: 1'].view.setValues(nearPlane=1.73807, 
    farPlane=2.40838, width=1.52603, height=0.698859, cameraPosition=(0.759053, 
    0.780174, 1.90333), cameraUpVector=(-0.111321, 0.922502, -0.369592), 
    cameraTarget=(0.510355, 0.00889367, 0.0531176), viewOffsetX=0.0104412, 
    viewOffsetY=0.165177)
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    contourStyle=CONTINUOUS, maxValue=1.16886, minValue=0.0182123)
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    spectrum='Blue to red')
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    spectrum='div: cool2warm')
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    spectrum='Black to white')
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    spectrum='Rainbow')
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    spectrum='seq: batlow')
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    spectrum='seq: turbo')
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(
    contourMethod=TESSELLATED)
#* The 'tessellated' contour method is not supported in combination with 
#* 'render beam profiles'.
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.showFirstFrame()
session.viewports[session.currentViewportName].animationController.showFrame(
    frame=57)
