import logging
import os

import vtk
import qt
import ctk

import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin


#
# VirtualEndoscopy
#

class VirtualEndoscopy(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Virtual Endoscopy"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["VirtualEndoscopy"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["Yi Gao(Shenzhen University), Haofan Huang (Shenzhen University)"]  # TODO: replace with "Firstname Lastname (Organization)"
    # TODO: update with short description of the module and a link to online module documentation
    self.parent.helpText = """
This extension performs the virtual endoscopy, including: centerline extraction, lumen segmentation, and virtual fly-through.

This is an example of scripted loadable module bundled in an extension.
See more information in <a href="https://github.com/organization/projectname#VirtualEndoscopy">module documentation</a>.
"""
    # TODO: replace with organization, grant and thanks
    self.parent.acknowledgementText = """
This work is partially suppored by the Key-Area Research and Development Program of Guangdong Province grant 2021B0101420005, the Key Technology Development Program of Shenzhen grant JSGG20210713091811036, the Department of Education of Guangdong Province grant 2017KZDXM072, the National Natural Science Foundation of China grant 61601302, the Shenzhen Key Laboratory Foundation grant ZDSYS20200811143757022, the Shenzhen Peacock Plan grant KQTD2016053112051497, and the SZU Top Ranking Project grant 86000000210.

This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""

    # Additional initialization step after application startup is complete
    slicer.app.connect("startupCompleted()", registerSampleData)


#
# Register sample data sets in Sample Data module
#

def registerSampleData():
  """
  Add data sets to Sample Data module.
  """
  # It is always recommended to provide sample data for users to make it easy to try the module,
  # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

  import SampleData
  iconsPath = os.path.join(os.path.dirname(__file__), 'Resources/Icons')

  # To ensure that the source code repository remains small (can be downloaded and installed quickly)
  # it is recommended to store data sets that are larger than a few MB in a Github release.

  # VirtualEndoscopy1
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='VirtualEndoscopy',
    sampleName='VirtualEndoscopy1',
    # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
    # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
    thumbnailFileName=os.path.join(iconsPath, 'VirtualEndoscopy1.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
    fileNames='VirtualEndoscopy1.nrrd',
    # Checksum to ensure file integrity. Can be computed by this command:
    #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
    checksums = 'SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95',
    # This node name will be used when the data set is loaded
    nodeNames='VirtualEndoscopy1'
  )

  # VirtualEndoscopy2
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='VirtualEndoscopy',
    sampleName='VirtualEndoscopy2',
    thumbnailFileName=os.path.join(iconsPath, 'VirtualEndoscopy2.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
    fileNames='VirtualEndoscopy2.nrrd',
    checksums = 'SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97',
    # This node name will be used when the data set is loaded
    nodeNames='VirtualEndoscopy2'
  )


#
# VirtualEndoscopyWidget
#

class VirtualEndoscopyWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)  # needed for parameter node observation
    self.logic = None
    self._parameterNode = None
    self._updatingGUIFromParameterNode = False




  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...
    self.setupWholePanel()


    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()


  def onComputeLumenButton(self):
    originalVolumeNode = self.inputOriginalImageSelector1.currentNode()
    inputAxisLabelNode = self.inputAxisLabelImageSelector2.currentNode()
    outputLumenLabelVolumeNode = self.outputLumenLabelImageSelector.currentNode()


    # run the filter
    ijkToRAS = vtk.vtkMatrix4x4()
    originalVolumeNode.GetIJKToRASMatrix(ijkToRAS)
    outputLumenLabelVolumeNode.SetIJKToRASMatrix(ijkToRAS)
    #outputLumenLabelVolumeNode.SetName("vesselLumenLabelVolume")

    parameters = {}
    parameters['inputVolume'] = originalVolumeNode.GetID()
    parameters['inputAxisLabelVolume'] = inputAxisLabelNode.GetID()
    parameters['lowerThreshold'] = self.thresholdSliderWidget.value
    parameters['outputLumenMaskVolume'] = outputLumenLabelVolumeNode.GetID()

    slicer.cli.run( slicer.modules.segmentlumenfromaxis, None, parameters, wait_for_completion=True )


  def setupWholePanel(self):
    wholeProcessCollapsibleButton = ctk.ctkCollapsibleButton()
    wholeProcessCollapsibleButton.text = "The Whole Process Pipeline"
    self.layout.addWidget(wholeProcessCollapsibleButton)

    # Layout within the dummy collapsible button
    wholeProcessFormLayout = qt.QFormLayout(wholeProcessCollapsibleButton)

    #
    # input volume selector
    #
    self.inputOriginalImageSelector = slicer.qMRMLNodeComboBox()
    self.inputOriginalImageSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputOriginalImageSelector.selectNodeUponCreation = True
    self.inputOriginalImageSelector.addEnabled = False
    self.inputOriginalImageSelector.removeEnabled = False
    self.inputOriginalImageSelector.renameEnabled = False
    self.inputOriginalImageSelector.noneEnabled = False
    self.inputOriginalImageSelector.showHidden = False
    self.inputOriginalImageSelector.showChildNodeTypes = False
    self.inputOriginalImageSelector.setMRMLScene( slicer.mrmlScene )
    self.inputOriginalImageSelector.setToolTip( "Pick the input to the algorithm." )
    wholeProcessFormLayout.addRow("Input Volume: ", self.inputOriginalImageSelector)

    # Input fiducials node selector
    self.treeFiducialsNodeSelector = slicer.qMRMLNodeComboBox()
    self.treeFiducialsNodeSelector.objectName = 'treeFiducialsNodeSelector'
    self.treeFiducialsNodeSelector.toolTip = "Select a fiducial list to define control points for the axis."
    self.treeFiducialsNodeSelector.nodeTypes = ['vtkMRMLMarkupsFiducialNode', 'vtkMRMLAnnotationHierarchyNode', 'vtkMRMLFiducialListNode']
    #self.treeFiducialsNodeSelector.nodeTypes = ['vtkMRMLMarkupsFiducialNode']
    self.treeFiducialsNodeSelector.noneEnabled = False
    self.treeFiducialsNodeSelector.renameEnabled = False
    self.treeFiducialsNodeSelector.addEnabled = False
    self.treeFiducialsNodeSelector.removeEnabled = False
    wholeProcessFormLayout.addRow("Axis Tree Fiducials:", self.treeFiducialsNodeSelector)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)', self.treeFiducialsNodeSelector, 'setMRMLScene(vtkMRMLScene*)') #not sure what this line is for???

    #
    # output lumen label image selector
    #
    self.outputLumenLabelImageSelector = slicer.qMRMLNodeComboBox()
    self.outputLumenLabelImageSelector.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.outputLumenLabelImageSelector.objectName = "lumenLabelVolume"
    self.outputLumenLabelImageSelector.selectNodeUponCreation = True
    self.outputLumenLabelImageSelector.renameEnabled = True
    self.outputLumenLabelImageSelector.addEnabled = True
    self.outputLumenLabelImageSelector.removeEnabled = False
    self.outputLumenLabelImageSelector.noneEnabled = False
    self.outputLumenLabelImageSelector.showHidden = False
    self.outputLumenLabelImageSelector.showChildNodeTypes = False
    self.outputLumenLabelImageSelector.setMRMLScene( slicer.mrmlScene )
    self.outputLumenLabelImageSelector.setToolTip( "Pick the output to the algorithm." )
    wholeProcessFormLayout.addRow("Lumen Label Volume: ", self.outputLumenLabelImageSelector)


    self.vesselBrighterCheckBox = qt.QCheckBox()
    self.vesselBrighterCheckBox.checked = False
    self.vesselBrighterCheckBox.setToolTip("Is vessel brighter than the surroundings?")
    wholeProcessFormLayout.addRow("Vessel brighter than the surroundings?", self.vesselBrighterCheckBox)

    #
    # Calcification threshold
    #
    self.calcificationSlicerWidget = ctk.ctkSliderWidget()
    self.calcificationSlicerWidget.singleStep = 1
    self.calcificationSlicerWidget.minimum = -1000
    self.calcificationSlicerWidget.maximum = 10000
    self.calcificationSlicerWidget.value = 800
    self.calcificationSlicerWidget.setToolTip("Calcification value.")
    wholeProcessFormLayout.addRow("Calcification", self.calcificationSlicerWidget)


    self.thresholdSliderWidget = ctk.ctkSliderWidget()
    self.thresholdSliderWidget.tracking = True
    self.thresholdSliderWidget.singleStep = 1
    self.thresholdSliderWidget.minimum = -1000
    self.thresholdSliderWidget.maximum = 1000
    self.thresholdSliderWidget.decimals = 0
    self.thresholdSliderWidget.value = 0
    self.thresholdSliderWidget.setToolTip("Threshold Value")
    self.thresholdSliderWidget.enabled = True
    wholeProcessFormLayout.addRow("Treshold", self.thresholdSliderWidget)

    self.superResolutionCheckBox = qt.QCheckBox()
    self.superResolutionCheckBox.checked = False
    self.superResolutionCheckBox.setToolTip("Perform super-resolution segmentation?")
    wholeProcessFormLayout.addRow("Perform super-resolution segmentation?", self.superResolutionCheckBox)

    #
    # Apply Button
    #
    self.computeWholeProcess = qt.QPushButton("Run it")
    self.computeWholeProcess.toolTip = "Run the whole process."
    self.computeWholeProcess.enabled = True
    wholeProcessFormLayout.addRow(self.computeWholeProcess)

    # connections
    self.computeWholeProcess.connect('clicked(bool)', self.onWholeProcessButton)

  def enableOrDisableComputeVesslenessButton(self):
    self.computeVesslenessButton.enabled = self.inputOriginalImageSelector.currentNodeID and self.outputVesslenessImageSelector.currentNodeID


  def onWholeProcessButton(self):
    #--------------------------------------------------------------------------------
    # Run the ComputeVesselness CLI to get the vesselness
    parameters = {}
    parameters['inputVolume'] = self.inputOriginalImageSelector.currentNode().GetID()
    parameters['vesselIsBrighter'] = self.vesselBrighterCheckBox.checked
    vesselNessImageNode = slicer.vtkMRMLScalarVolumeNode()
    vesselNessImageNode.SetName("vesselNessImage")
    slicer.mrmlScene.AddNode( vesselNessImageNode )
    parameters['outputVesselnessVolume'] = vesselNessImageNode.GetID()
    parameters['sigma'] = 1.0
    parameters['alpha1'] = 2.0
    parameters['alpha2'] = 2.0
    parameters['calcificationThreshold'] = self.calcificationSlicerWidget.value
    slicer.cli.run( slicer.modules.computevesselness, None, parameters, wait_for_completion=True )
    # Run the ComputeVesselness CLI to get the vesselness
    #================================================================================


    #--------------------------------------------------------------------------------
    # Compute axis tree
    parameters = {}
    parameters['inputVesselnessVolume'] = vesselNessImageNode.GetID()
    parameters['fiducialsAlongCA'] = self.treeFiducialsNodeSelector.currentNode()

    axisTreeLabelImageNode = slicer.vtkMRMLLabelMapVolumeNode()
    axisTreeLabelImageNode.SetName("axisTreeLabelImage")
    slicer.mrmlScene.AddNode( axisTreeLabelImageNode )
    parameters['outputAxisMaskVolume'] = axisTreeLabelImageNode.GetID()

    ################################################################################
    # This CLI module "computetree" treat the 1st fiducial point as
    #the root and trace every other fiducial points back to the root
    #slicer.cli.run( slicer.modules.computetree, None, parameters, wait_for_completion=True )

    ################################################################################
    # This CLI module "ComputeAxisFromVesselness" treat the (i-1)-th
    #fiducial point as the root and trace from the i-th point to the
    #(i-1)-th one. Then from (i+1) to i, and on. So this module gives
    #a SINGEL axis.
    #slicer.cli.run( slicer.modules.ComputeAxisFromVesselness, None, parameters, wait_for_completion=True )
    slicer.cli.run( slicer.modules.computeaxisfromvesselness, None, parameters, wait_for_completion=True )
    #================================================================================

    #--------------------------------------------------------------------------------
    # Run the SegmentLumenFromAxis CLI to get the vesselness
    parameters = {}
    parameters['inputVolume'] = self.inputOriginalImageSelector.currentNode().GetID()
    parameters['inputAxisLabelVolume'] = axisTreeLabelImageNode.GetID()
    parameters['lowerThreshold'] = self.thresholdSliderWidget.value
    parameters['calcificationThreshold'] = self.calcificationSlicerWidget.value
    parameters['outputLumenMaskVolume'] = self.outputLumenLabelImageSelector.currentNode()
    parameters['superResolution'] = self.superResolutionCheckBox.checked

    slicer.cli.run( slicer.modules.segmentlumenfromaxis, None, parameters, wait_for_completion=True )
    # Run the SegmentLumenFromAxis CLI to get the vesselness
    #================================================================================

    # d = vesselNessImageNode.GetDisplayNode()
    # d.SetVisibility(0) # do not show the vesselness image



    
    # #--------------------------------------------------------------------------------
    # # Run the model maker to genreate surface
    # parameters = {}
    # parameters['InputVolume'] = self.outputLumenLabelImageSelector.currentNode().GetID()
    # lumenModelHierarchy = slicer.vtkMRMLModelHierarchyNode()
    # lumenModelHierarchy.SetName("theLumenModelHierarchy")
    # slicer.mrmlScene.AddNode( lumenModelHierarchy )
    # parameters['ModelSceneFile'] = lumenModelHierarchy.GetID()
    # slicer.cli.run( slicer.modules.modelmaker, None, parameters, wait_for_completion=True )
    # # Run the model maker to genreate surface
    # #================================================================================



    # #--------------------------------------------------------------------------------
    # # Backface To Frontface 
    # e = slicer.util.getNode('Model_1_jake')
    # f = e.GetDisplayNode()
    # f.SetBackfaceCulling(0)
    # f.SetFrontfaceCulling(1)
    # # Backface To Frontface 
    # #================================================================================



    # #--------------------------------------------------------------------------------
    # # To get the central coordinates of pancreatic duct and insert fiducials on them
    # inputvalues = slicer.util.getNode('axisTreeLabelImage')
    # # print('ehhhhh.')
    # # print(inputvalues)
    # inputdata = slicer.util.arrayFromVolume(inputvalues)
    # coordinates = np.where(inputdata != 0)
    # I = coordinates[2][3::3]
    # # print('I:')
    # # print(I)
    # J = coordinates[1][3::3]
    # # print('J:')
    # # print(J)
    # K = coordinates[0][3::3]
    # # print('K:')
    # # print(K)
    # RasCoordinates = []
    # N = 0
    
    # # IJK TO RAS
    # ijkToRasMatrix = vtk.vtkMatrix4x4()
    # inputvalues.GetIJKToRASMatrix(ijkToRasMatrix)	
    # for N in range(len(I)):
    #   c_Ijk = [I[N],J[N],K[N],1]
    #   c_Ras = ijkToRasMatrix.MultiplyFloatPoint(c_Ijk)
    #   RasCoordinates.append(c_Ras)
    # # print('RasCoordinates:')
    # # print(RasCoordinates)

    # # Change the bottom fiducial point to the second point
    # fidNode = slicer.util.getNode('F')
    # ras = [0,0,0]
    # fidNode.GetNthFiducialPosition(1,ras)
    # # print(ras)
    # RasCoordinates.append(ras)
    # fidNode.SetNthFiducialPosition(1,RasCoordinates[0][0],RasCoordinates[0][1],RasCoordinates[0][2])
    # RasCoordinates.pop(0)
    # for c in RasCoordinates:
    #   slicer.modules.markups.logic().AddFiducial(c[0],c[1],c[2])
    # # To get the central coordinates of pancreatic duct and insert fiducials on them
    # #================================================================================
   




    # TODO
    # parameters = {}
    # parameters["InputVolume"] = vesselNessImage.GetID()
    # grayModel = slicer.vtkMRMLModelNode()
    # slicer.mrmlScene.AddNode( grayModel )
    # parameters["OutputGeometry"] = grayModel.GetID()
    # parameters['Threshold'] = self.vesselSlicerWidget.value
    # grayMaker = slicer.modules.grayscalemodelmaker
    # slicer.cli.runSync(grayMaker, None, parameters)
    # d = grayModel.GetDisplayNode()
    # d.SetVisibility(0) # do not show the gray model


    # parameters = {}
    # parameters["ModelSceneFile"] = grayModel.GetID()
    # parameters["fiducialsAlongCA"] = self.axisFiducialsNodeSelectorNew.currentNode().GetID()
    # parameters["axisPolylineName"] = self.outputGeometrySelector.currentNode().GetID()
    # slicer.cli.run( slicer.modules.computepathonsurface, None, parameters, wait_for_completion=True )

  def onNewButton(self):

    parameters = {}
    parameters['inputVolume'] = self.inputOriginalImageSelector.currentNode().GetID()
    parameters['vesselIsBrighter'] = self.vesselBrighterCheckBox.checked
    vesselNessImage = slicer.vtkMRMLScalarVolumeNode()
    vesselNessImage.SetName("vesselNessImage")
    slicer.mrmlScene.AddNode( vesselNessImage )
    parameters['outputVesselnessVolume'] = vesselNessImage.GetID()
    parameters['sigma'] = 1.0
    parameters['alpha1'] = 2.0
    parameters['alpha2'] = 2.0
    parameters['calcificationThreshold'] = self.calcificationSlicerWidget.value
    slicer.cli.run( slicer.modules.computevesselness, None, parameters, wait_for_completion=True )

    parameters = {}
    parameters["InputVolume"] = vesselNessImage.GetID()
    grayModel = slicer.vtkMRMLModelNode()
    slicer.mrmlScene.AddNode( grayModel )
    parameters["OutputGeometry"] = grayModel.GetID()
    parameters['Threshold'] = self.vesselSlicerWidget.value
    grayMaker = slicer.modules.grayscalemodelmaker
    slicer.cli.runSync(grayMaker, None, parameters)
    d = grayModel.GetDisplayNode()
    d.SetVisibility(0) # do not show the gray model


    parameters = {}
    parameters["ModelSceneFile"] = grayModel.GetID()
    parameters["fiducialsAlongCA"] = self.axisFiducialsNodeSelectorNew.currentNode().GetID()
    parameters["axisPolylineName"] = self.outputGeometrySelector.currentNode().GetID()
    slicer.cli.run( slicer.modules.computepathonsurface, None, parameters, wait_for_completion=True )


  def onComputeBranchButton(self):
    vesselnessVolumeNode = self.inputVesslenessImageSelector2.currentNode()
    inputAxisLabelNode = self.inputAxisLabelImageSelector1.currentNode()
    mainAndBranchAxisLabelImageNode = self.outputBranchLabelImageSelector.currentNode()
    #self.mainAndBranchAxisLabelImageNode.SetName("axisWithBranchLabelVolume-label")
    if not (vesselnessVolumeNode and mainAndBranchAxisLabelImageNode):
      qt.QMessageBox.critical(
        slicer.util.mainWindow(),
        'Compute', 'Input and output volumes are required for computing branch axis')
      return
    # run the filter
    ijkToRAS = vtk.vtkMatrix4x4()
    vesselnessVolumeNode.GetIJKToRASMatrix(ijkToRAS)
    mainAndBranchAxisLabelImageNode.SetIJKToRASMatrix(ijkToRAS)

    branchFiducialsNode = self.branchFiducialsNodeSelector.currentNode();

    parameters = {}
    parameters['inputVesselnessVolume'] = vesselnessVolumeNode.GetID()
    parameters['fiducialsOnBranches'] = branchFiducialsNode.GetID()
    parameters['inputAxisLabelVolume'] = inputAxisLabelNode.GetID()
    parameters['outputAxisAndBranchMaskVolume'] = mainAndBranchAxisLabelImageNode.GetID()

    slicer.cli.run( slicer.modules.computebranchaxis, None, parameters, wait_for_completion=True )

  def onPreProcessImageButton(self):
    originalVolumeNode = self.inputOriginalImageSelector.currentNode()
    vesselnessVolumeNode = self.outputVesslenessImageSelector.currentNode()
    #vesselnessVolumeNode.SetName("vesselnessVolume")
    if not (originalVolumeNode and vesselnessVolumeNode):
      qt.QMessageBox.critical(
          slicer.util.mainWindow(),
          'Compute', 'Input and output volumes are required for computing vessleness')
      return
    # run the filter
    ijkToRAS = vtk.vtkMatrix4x4()
    originalVolumeNode.GetIJKToRASMatrix(ijkToRAS)
    vesselnessVolumeNode.SetIJKToRASMatrix(ijkToRAS)

    parameters = {}
    parameters['inputVolume'] = originalVolumeNode.GetID()
    parameters['vesselIsBrighter'] = self.vesselBrighterCheckBox.checked
    parameters['outputVesselnessVolume'] = vesselnessVolumeNode.GetID()
    parameters['sigma'] = 0.5
    parameters['alpha1'] = 0.5
    parameters['alpha2'] = 2.0
    parameters['calcificationThreshold'] = self.calcificationSlicerWidget.value
    # parameters['sigma'] = self.sigmaSlicerWidget.value
    # parameters['alpha1'] = self.alpha1SlicerWidget.value
    # parameters['alpha2'] = self.alpha2SlicerWidget.value

    slicer.cli.run( slicer.modules.computevesselness, None, parameters, wait_for_completion=True )


    # fill computed vesselness image to other UI selectors
    self.inputVesslenessImageSelector1.setCurrentNode(self.outputVesslenessImageSelector.currentNode());
    self.inputVesslenessImageSelector2.setCurrentNode(self.outputVesslenessImageSelector.currentNode());


    appLogic = slicer.app.applicationLogic()
    selectionNode = appLogic.GetSelectionNode()
    selectionNode.SetReferenceActiveVolumeID(self.inputOriginalImageSelector.currentNode().GetID())
    #selectionNode.SetReferenceSecondaryVolumeID(fg)
    appLogic.PropagateVolumeSelection()

#    self.computeAxisButton.enabled = True



  def onComputeAxisButton(self):
    inputVesselnessVolumeNode = self.inputVesslenessImageSelector1.currentNode()
    axisLabelVolumeNode = self.outputAxisLabelImageSelector.currentNode()
    #self.axisLabelVolumeNode.SetName("axisLabelVolume-label")
    if not (inputVesselnessVolumeNode and axisLabelVolumeNode):
      qt.QMessageBox.critical(
        slicer.util.mainWindow(),
        'Compute', 'Input and output volumes are required for computing axis')
      return
    # run the filter
    ijkToRAS = vtk.vtkMatrix4x4()
    inputVesselnessVolumeNode.GetIJKToRASMatrix(ijkToRAS)
    axisLabelVolumeNode.SetIJKToRASMatrix(ijkToRAS)

    axisFiducialsNode = self.axisFiducialsNodeSelector.currentNode();

    parameters = {}
    parameters['inputVesselnessVolume'] = inputVesselnessVolumeNode.GetID()
    parameters['fiducialsAlongCA'] = axisFiducialsNode.GetID()
    parameters['outputAxisMaskVolume'] = axisLabelVolumeNode.GetID()


    slicer.cli.run( slicer.modules.computeaxisfromvesselness, None, parameters, wait_for_completion=True )


  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.removeObservers()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()

  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)

  def onSceneEndClose(self, caller, event):
    """
    Called just after the scene is closed.
    """
    # If this module is shown while the scene is closed then recreate a new parameter node immediately
    if self.parent.isEntered:
      self.initializeParameterNode()

  def initializeParameterNode(self):
    """
    Ensure parameter node exists and observed.
    """
    # Parameter node stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.

    self.setParameterNode(self.logic.getParameterNode())

    # Select default input nodes if nothing is selected yet to save a few clicks for the user
    if not self._parameterNode.GetNodeReference("InputVolume"):
      firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
      if firstVolumeNode:
        self._parameterNode.SetNodeReferenceID("InputVolume", firstVolumeNode.GetID())

  def setParameterNode(self, inputParameterNode):
    """
    Set and observe parameter node.
    Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
    """

    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    # Update node selectors and sliders
    self.ui.inputSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputVolume"))
    self.ui.outputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputVolume"))
    self.ui.invertedOutputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputVolumeInverse"))
    self.ui.imageThresholdSliderWidget.value = float(self._parameterNode.GetParameter("Threshold"))
    self.ui.invertOutputCheckBox.checked = (self._parameterNode.GetParameter("Invert") == "true")

    # Update buttons states and tooltips
    if self._parameterNode.GetNodeReference("InputVolume") and self._parameterNode.GetNodeReference("OutputVolume"):
      self.ui.applyButton.toolTip = "Compute output volume"
      self.ui.applyButton.enabled = True
    else:
      self.ui.applyButton.toolTip = "Select input and output volume nodes"
      self.ui.applyButton.enabled = False

    # All the GUI updates are done
    self._updatingGUIFromParameterNode = False

  def updateParameterNodeFromGUI(self, caller=None, event=None):
    """
    This method is called when the user makes any change in the GUI.
    The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID("InputVolume", self.ui.inputSelector.currentNodeID)
    self._parameterNode.SetNodeReferenceID("OutputVolume", self.ui.outputSelector.currentNodeID)
    self._parameterNode.SetParameter("Threshold", str(self.ui.imageThresholdSliderWidget.value))
    self._parameterNode.SetParameter("Invert", "true" if self.ui.invertOutputCheckBox.checked else "false")
    self._parameterNode.SetNodeReferenceID("OutputVolumeInverse", self.ui.invertedOutputSelector.currentNodeID)

    self._parameterNode.EndModify(wasModified)

  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

      # Compute output
      self.logic.process(self.ui.inputSelector.currentNode(), self.ui.outputSelector.currentNode(),
        self.ui.imageThresholdSliderWidget.value, self.ui.invertOutputCheckBox.checked)

      # Compute inverted output (if needed)
      if self.ui.invertedOutputSelector.currentNode():
        # If additional output volume is selected then result with inverted threshold is written there
        self.logic.process(self.ui.inputSelector.currentNode(), self.ui.invertedOutputSelector.currentNode(),
          self.ui.imageThresholdSliderWidget.value, not self.ui.invertOutputCheckBox.checked, showResult=False)


#
# VirtualEndoscopyLogic
#

class VirtualEndoscopyLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    """
    Called when the logic class is instantiated. Can be used for initializing member variables.
    """
    ScriptedLoadableModuleLogic.__init__(self)



  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def run(self, inputVolume, outputVolume, imageThreshold, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
    cliParams = {'InputVolume': inputVolume.GetID(), 'OutputVolume': outputVolume.GetID(), 'ThresholdValue' : imageThreshold, 'ThresholdType' : 'Above'}
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

    # Capture screenshot
    if enableScreenshots:
      self.takeScreenshot('VirtualEndoscopyTest-Start','MyScreenshot',-1)

    logging.info('Processing completed')

    return True


#
# VirtualEndoscopyTest
#

class VirtualEndoscopyTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_VirtualEndoscopy1()

  def test_VirtualEndoscopy1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    # Get/create input data

    import SampleData
    registerSampleData()
    inputVolume = SampleData.downloadSample('VirtualEndoscopy1')
    self.delayDisplay('Loaded test data set')

    inputScalarRange = inputVolume.GetImageData().GetScalarRange()
    self.assertEqual(inputScalarRange[0], 0)
    self.assertEqual(inputScalarRange[1], 695)

    outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    threshold = 100

    # Test the module logic

    logic = VirtualEndoscopyLogic()

    # Test algorithm with non-inverted threshold
    logic.process(inputVolume, outputVolume, threshold, True)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.process(inputVolume, outputVolume, threshold, False)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')


