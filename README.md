# multiecho-fMRI-T2-mapping
General information
The tool consists of 2 components:
1.	temporal voxel-by-voxel noise filtering of 3 gradient echoes acquired as 3-D multi-volume fMRI data
2.	voxel-by-voxel estimation of T2* at all acquired repetition times
These 2 steps are not independent: Step 2, the T2* estimation, can not be done before the input data have been processed by Step 1, the temporal voxel-by-voxel noise filtering.
The input data must contain three .nii files named *Echo1*.nii, *Echo2*.nii and *Echo3*.nii files containing fMRI echo series, as well as three .mat files named *Echo1*.mat, *Echo2*.mat  and *Echo3*.mat containing DICOM information related to *Echo1*.nii, *Echo2*.nii and *Echo3*.nii.
Location of the datasets: the datasets used so far are placed on the CBIA gryf data server in the directory /data/common/PROJECTS/_Current/CEITEC_T2star/T2star_data. 
An example dataset  that can be immediately used by the CBIA as input to temporal voxel-by-voxel noise filtering is located at /data/common/PROJECTS/_Current/CEITEC_T2star/T2star_data/Example
Installation of the tool: simply place the T2star github directory from the T2star github.zip file to a location of your choice.
Startup of the tool:
•	launch MATLAB >=2016b
•	change to your MFL directory
The tool was developed and tested using MATLAB 2016b under Windows 7 and WIndows 10.

Running the MATLAB implementation of the denoising algorithm.
The denoising algorithm is encapsulated in the script named T2star_TV_MFL.m which calls several other scripts and one mex-function. 
On completion of each echo’s denoising, the algorithm creates a summary plot which helps pinpoint parameters that should be adjusted to obtain optimum convergence, as explained further on:
![image](https://user-images.githubusercontent.com/25533528/224768176-e8235303-06a9-49e7-8094-40d05431ff7b.png)

The denoising script is launched by the command
>> T2star_TV_github
after which a self-explanatory dialog follows. _

Algorithm notifications are usually introduced by “CBIA>>” to distinguish them from system messages.
CBIA>>please make sure spm12 toolbox is on your MATLAB path
Prompts for user actions are mostly terminated by a colon:
CBIA>>please enter full directory path containing
      xxx_Echo_N.nii and
      yyy_Echo_N_dicom_header.mat,
      such as I:\Data\T2star\4055B\FCNI1_mag\cwu:
while informative messages about the algorithm progress begin with an ellipsis (...):
CBIA>>...unfiltered echoes will be read from
         I:\Data\T2star\4055B\FCNI1_mag\cwu
CBIA>>...reading 4-D echoes Echo_1, Echo_2, Echo_3 
         and Echo_1_dicom_header, Echo_2_dicom_header, Echo_3_dicom_header 
CBIA>>...reading TE(1),TE(2),TE(3) and TR from dicom headers 
CBIA>>...concatenating 4-D echoes Echo_1, Echo_2, Echo_3 
         to 5-D datacube Echoes with echo number as the 5-th dimension
--------------------------------------------------------------------
Denoising is performed on 3-D volumes of all 3 echo trains for all fMRI-acquired time instants. 
A z-slice is selected here for 2-D visual representation only, and does not affect denoising:
CBIA>>pick a z-slice, or press RETURN for default: [23]:
--------------------------------------------------------------------
To speed up computations by parallel processing on multi-core CPUs, the number of threads of the user’s CPU can be entered:
CBIA>>enter number of CPU threads, or press RETURN for 8 threads:
--------------------------------------------------------------------
 
Setting the algorithm parameters:

CBIA>>enter log2(mu), or press RETURN for default [-10]:
The parameter mu=2^log2(mu) defines weighting of the L2-norm (mean square error) between the measured and the denoised echo time-course: the smaller mu, the smoother, in terms of total variation TV, will the denoised signal be. mu affects the speed and the stability of the TV optimization process.
![image](https://user-images.githubusercontent.com/25533528/224766207-f321e455-71bf-41b7-956e-86d3aae8978b.png)
--------------------------------------------------------------------
CBIA>>enter log2(beta), or press RETURN for default [-4]:
The parameter beta=2^log2(beta)affects the speed with which the true total variation TV will approach its approximation w_norm: the larger beta, the faster will TV and w_norm meet. beta affects the speed and the stability of the TV optimization process.
![image](https://user-images.githubusercontent.com/25533528/224766125-b2dabeb0-f510-4bf9-833e-6a633f70cc45.png)
--------------------------------------------------------------------
CBIA>>enter nr_shrinks, or press RETURN for default [10]:
The parameter nr_shrinks defines the number of shrinkage iterations, i.e.those in which total variation TV is reduced. nr_shrinks affects the accuracy of the TV optimization process.
CBIA>>enter nr_grads, or press RETURN for default [5]:
The parameter nr_grads defines the number of steepest descent operations, i.e.those in which square error data_term between the filtered and the noisy echo is reduced. nr_grads affects the accuracy and the stability of the data_term optimization process.
The product nr_shrinks*nr_grads should be large enough for the data_term to reach the plateau:
![image](https://user-images.githubusercontent.com/25533528/224765965-0bdf0dd3-2618-48bf-9ba3-b99579e58133.png)
Also, the product nr_shrinks*nr_grads should be large enough for the primal and the dual functional to meet:
![image](https://user-images.githubusercontent.com/25533528/224765890-5ea4b66a-d3b3-42d1-9ec9-5120c1314b54.png)
--------------------------------------------------------------------
CBIA>>...TV filtering results will be stored in
         I:\Data\T2star\4055B\FCNI1_mag\cwu\TV_filtered_echoes\-16   -6   10   10
--------------------------------------------------------------------
The algorithm progress is reported:
CBIA>>...temporal TV filtering of Echo_1
shrink_count=1
shrink_count=2
shrink_count=3
shrink_count=4
shrink_count=5
shrink_count=6
shrink_count=7
shrink_count=8
shrink_count=9
shrink_count=10
ElapsedTime= 176.0191
CBIA>>...temporal TV filtering of Echo_2
shrink_count=1
shrink_count=2
shrink_count=3
shrink_count=4
shrink_count=5
shrink_count=6
shrink_count=7
shrink_count=8
shrink_count=9
shrink_count=10
ElapsedTime= 176.7533
CBIA>>...temporal TV filtering of Echo_3
shrink_count=1
shrink_count=2
shrink_count=3
shrink_count=4
shrink_count=5
shrink_count=6
shrink_count=7
shrink_count=8
shrink_count=9
shrink_count=10
ElapsedTime= 178.1268
CBIA>>...saving TV filtering results to
         I:\Data\T2star\4055B\FCNI1_mag\cwu\TV_filtered_echoes\-16   -6   10   10
--------------------------------------------------------------------
after which the termination notification appears:
CBIA>>...temporal TV filtering of the 5D datacube has been completed.
--------------------------------------------------------------------
Finally the user is prompted to peruse the denoising results:
         To visually compare unfiltered and TV-filtered echoes at voxels of your choice,
         pick repeatedly those you are interested in.
![image](https://user-images.githubusercontent.com/25533528/224765732-7d89def7-37df-4133-a2b5-ff3cc007993d.png)
pick_another_pixel:[return/0]
pick_another_pixel:[return/0] 0
>>
>>2.  voxel-by-voxel estimation of T2* at all acquired repetition times.
2.1  Estimation of T2* by weighted exponential matching of the three TV-filtered gradient echoes.
The MATLAB tool is based on the original implementation of the algorithm described in Michálek et al.:“Fast and accurate compensation of signal offset for T2 mapping“, https://link.springer.com/article/10.1007/s10334-019-00737-3, which processed only single 2D slices of the input MRI data. For use with 3D fMRI data, the original code has been enhanced to include also the z-dimension (volume) and t-dimension (time). Resulting T2* estimates and noise masks are output as 4-D .nii files. 
2.2  Running the T2*estimation MATLAB script.
The T2*estimation script is launched by the command
>> T2star_series_WLS_MFL
after which a self-explanatory dialog consisting of informative messages (introduced with ...) or prompts for user action (terminated with:) follows:
CBIA>>...please make sure spm12 toolbox is on your MATLAB path
CBIA>>please enter full directory path containing TV-filtered echoes
      TV_Echo_1.nii, TV_Echo_2.nii and TV_Echo_3.nii
      such as I:\Data\T2star\4055B\FCNI1_mag\cwu\TV_filtered_echoes\-16   -6   10   10:
CBIA>>...TV-filtered echoes will be read from
         I:\Data\T2star\4055B\FCNI1_mag\cwu\TV_filtered_echoes\-16   -6   10   10
CBIA>>...reading 4-D TV-filtered echoes TV_Echo_1, TV_Echo_2, TV_Echo_3 
CBIA>>...concatenating 4-D TV-filtered echoes TV_Echo_1, TV_Echo_2, TV_Echo_3 
         to 5-D datacube Echoes with echo number as the 5-th dimension
CBIA>>...unfiltered echoes will be read from
         I:\Data\T2star\4055B\FCNI1_mag\cwu
CBIA>>...reading 4-D unfiltered echoes Echo_1, Echo_2, Echo_3 
         and Echo_1_dicom_header, Echo_2_dicom_header, Echo_3_dicom_header 
CBIA>>...concatenating 4-D unfiltered echoes Echo_1, Echo_2, Echo_3 
         to 5-D datacube Ref_echoes with echo number as the 5-th dimension
CBIA>>...reading TE(1),TE(2),TE(3) and TR from dicom headers 
CBIA>>maximum slice:46
--------------------------------------------------------------------
T2* estimation is performed on whole 3-D volumes for all fMRI acquired time instants. 
z-slice is selected here for 2-D visual representation only, and does not affect denoising:
CBIA>>pick a z-slice for 2D visualizations, or press RETURN for default: [23]:
slc=23
--------------------------------------------------------------------
Before starting the T2* estimation, a gray-value threshold separating background pixels from valid fMRI information is generated automatically. The threshold generation relies on the assumption that the histogram of Echo_1 array (i.e. the brightest echo) exhibits two distinct peaks. The user is presented a plot with the histograms as well as the threshold-generated background mask, and is given the possibility to change the threshold:
CBIA>>press RETURN to accept the threshold, or enter a different one: [12692.5]:
![image](https://user-images.githubusercontent.com/25533528/224765542-4a437e36-d625-42a8-a938-453009efebde.png)
--------------------------------------------------------------------
CBIA>>...interleaving 4-D TV-filtered echoes TV_Echo_1, TV_Echo_2, TV_Echo_3 
         according to their absolute acquisition time 
         to 4-D datacube interleaved_echoes
--------------------------------------------------------------------
CBIA>>...calculating T2star values for all volumes of input data 
--------------------------------------------------------------------
The total execution time needed to estimate T2* for all voxels at all measured fMRI volumes is reported
ElapsedTime=34.8382
--------------------------------------------------------------------
and the directory where the 4-D results are saved as nifti files is displayed:
CBIA>>...saving noise mask and unmasked/masked T2star values as 4-D .nii files to 
         I:\Data\T2star\4055B\FCNI1_mag\cwu\TV_filtered_echoes\-16   -6   10   10\T2star_series_WLS_MFL
--------------------------------------------------------------------
For the user-selected slice, some raw/limited/masked  3-D T2* maps and a noise mask are generated, visualized  and saved as multipage tiffs in a slice-related directory. These can later be viewed, e.g., using ImageJ:
CBIA>>...extracting 3-D data at your selected slice number 23 for visualization
masked_pixel_ratio = 0.41266
CBIA>>...saving raw/limited/masked 3-D T2stars and noise mask at your selected slice number in
I:\Data\T2star\4055B\FCNI1_mag\cwu\TV_filtered_echoes\-16   -6   10   10\T2star_series_WLS_MFL\slc23
--------------------------------------------------------------------
And finally the user is given the possibility  to to inspect the T2* time courses at voxels of their choice:
CBIA>>...To visually inspect echo and T2star time courses at voxels of your choice,
         pick repeatedly those you are interested in.
![image](https://user-images.githubusercontent.com/25533528/224765213-2efb72ee-821c-47e1-bee6-2c651f9c3ddc.png)
![image](https://user-images.githubusercontent.com/25533528/224765326-0b5472af-d623-4c31-907e-b52e04889a61.png)
         pick_another_pixel:[press RETURN for YES/0 for NO] 
         pick_another_pixel:[press RETURN for YES/0 for NO] 
         pick_another_pixel:[press RETURN for YES/0 for NO] 0



