# T2*-estimation MATLAB Tool

This repository contains MATLAB files necessary to compute dynamic fMRI maps as described in a paper to appear in ... The user guide is in the file 
"MATLAB T2star Estimation Tool github guidelines.pdf"

The tool estimates T2* maps from three series of gradient echoes acquired in sessions in multiple data volumes with a constant repetition time and three constant echo times. 
The tool consists of 2 components:
1.	temporal voxel-by-voxel noise filtering of 3 gradient echoes acquired as 3D-t multi-volume fMRI data
2.	voxel-by-voxel estimation of T2* at all acquired repetition times

These 2 steps are not independent: Step 2, the T2* estimation, can not be done before the input data have been processed by Step 1, the temporal voxel-by-voxel noise filtering.

The input data must contain

1.	three .nii files named *Echo1*.nii, *Echo2*.nii and *Echo3*.nii files containing fMRI multi-echo series, as well as 
2.	three .mat files named *Echo1*.mat, *Echo2*.mat  and *Echo3*.mat containing DICOM headers  belonging  to *Echo1*.nii, *Echo2*.nii and *Echo3*.nii.

The .mat files *Echo1*.mat, *Echo2*.mat  and *Echo3*.mat containing DICOM headers must be loaded into MATLAB.  In the MATLAB workspace, the structure "dicom_header" must appear  and must contain the following values

1.	dicom_header.RepetitionTime
2.	dicom_header.EchoTime
