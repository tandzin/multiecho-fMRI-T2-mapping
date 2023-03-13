# multiecho-fMRI-T2-mapping
General information
The tool consists of 2 components:
1.	temporal voxel-by-voxel noise filtering of 3 gradient echoes acquired as 3-D multi-volume fMRI data
2.	voxel-by-voxel estimation of T2* at all acquired repetition times
These 2 steps are not independent: Step 2, the T2* estimation, can not be done before the input data have been processed by Step 1, the temporal voxel-by-voxel noise filtering.
The input data must contain three .nii files named *Echo1*.nii, *Echo2*.nii and *Echo3*.nii files containing fMRI echo series, as well as three .mat files named *Echo1*.mat, *Echo2*.mat  and *Echo3*.mat containing DICOM information related to *Echo1*.nii, *Echo2*.nii and *Echo3*.nii.
Location of the datasets: the datasets used so far are placed on the CBIA gryf data server in the directory /data/common/PROJECTS/_Current/CEITEC_T2star/T2star_data. 
An example dataset  that can be immediately used by the CBIA as input to temporal voxel-by-voxel noise filtering is located at /data/common/PROJECTS/_Current/CEITEC_T2star/T2star_data/Example
Installation of the tool: simply place the MFL directory from the MFL.zip file to a location of your choice.
Startup of the tool:
•	launch MATLAB >=2016b
•	change to your MFL directory
The tool was developed and tested using MATLAB 2016b under Windows 7 and WIndows 10.
