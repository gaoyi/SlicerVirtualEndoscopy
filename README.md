# SlicerVirtualEndoscopy

## Introduction

SlicerVirtualEndoscopy is an extension for simulating Virutla Endoscopy from volumetric image data in 3D Slicer.

Virtual endoscopy extracts the lumen from volumetric images and creates a fly-through visualization similar to the real endoscopy examinations.

Slicer already has an Endoscopy module, based on Delphine Nain’s work: https://dspace.mit.edu/handle/1721.1/87240

The current Endoscopy module requires the already extracted surface as input. The proposed module, however, starts from volumetric data and can perform vessel segmentation etc.

## Build and Installation
Following the standard exntesion build directions to build this one.

## Modules
This extension currently has two modules:
* VirtualEndoscopy: This is the main module of this extension. It calls the following modules to enable the entire computation pipeline.
* ComputeVesselness: From a scalar volumetric image, compute the "vesselness image", in which the pixel more likely to be inside a vessel is brighter than thouse aren't.
* ComputeAxisFromVesselness: From a few fiducial points user gives, trace in the vesselness image and form a axis of the vessel, as a binary image.
* SegmentLumenFromAxis: From the binary label image representing the axis of the vessel, segment the lumen of the vessel, and form a binary image for the lumen.

After the above lumen segmentation step, the surface of the lumen could be extracted and the Endoscopy module in Slicer could be used to simulate the fly through.

## Usage

### Example data
Example abdominal CT data from the TCIA pancrease image dataset:
`Roth, H.R., et al., Data from pancreas-CT. The cancer imaging archive, 2016. 32.`

### Module panel. The module panel of the VIrtual Endoscopy looks like:
![image](https://user-images.githubusercontent.com/920557/176148714-c95592a1-cf31-4770-8e24-d8997d2f657d.png)

Adding few fiducial points along the vessel (or tubular structure, the pancreatic duct in this case):
![image](https://user-images.githubusercontent.com/920557/176149073-04abdbd0-64e9-4d26-8134-095a21f37dba.png)

Run the module and it will extract the vessel/duct axis as well as the lumen of the vessel/duct:
![image](https://user-images.githubusercontent.com/920557/176149685-b1e8f01c-1178-4a1e-90bc-dc859e63c05e.png)


## Citation

If you find this extension helpful please cite this paper:

Haofan Huang, Xiaxia Yu, Mu Tian, Weizhen He, Shawn Xiang Li, Zhengrong Liang, and Yi Gao, "Open-source algorithm and software for CT-based virtual pancreatoscopy and other applications" in Visual Computing for Industry, Biomedicine, and Art, 2022, In print


