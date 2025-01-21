# TextureTensorAnalysis
Epithelial tissue deformation analysis by texture tensors

## Description

This repository provides scripts for quantifying epithelial tissue deformation using texture tensors and validating kinematic relationships for tissue morphogenesis. For detailed methodology, refer to [1].

The scripts are based on the texture tensor analysis method originally developed by Guirao et al. [2] and further refined in [1]. This method calculates the deformation rate tensor of epithelial tissues and decomposes it into components caused by different morphogenetic cell events, such as cell shape change, rearrangement, cell division and death/apoptosis.

The repository includes scripts for validating kinematic relationships for epithelial morphogenesis as proposed in [3].

Refer to [1, 2] for detailed explanations.

## Requirement

* matlab
* distrib_computing_toolbox
* image_toolbox
* sensor_fusion_and_tracking

## Usage
* **Step 0: Generate Vector Field Data**  
Use `PIVlab`, an open-source MATLAB tool [4], to generate vector field data by Particle Image Velocimetry (PIV). In the output, you have the following variables: (x, y, u_filtered, v_filtered).
As a sample data, PIVlab is applied to the image file `input/Color-Cregion-WT3.tif`.

* **Step 1: Generate Cell IDs and Tracking Data**  
Run `Main_step1.m` to create cell IDs and their tracking data. This script processes the skeletonized image `input/Skeleton-Cregion-WT3.tif` and the PIV data `input/PIVlab.mat`. Upon completion, it outputs `pos_link_step1.mat`.  
Note: If tracking errors occur, perform manual corrections.

* **Step 2: Calculate Time Series of Strains**  
Run `Main_step2.m` to calculate time series data of texture tensors $\hat{\bf M}^{(1)}$, deformation gradient tensors $\hat{\bf F}$, and strain tensors corresponding to various cellular events (denoted by $\hat{\bf S}$, $\hat{\bf R}$, $\hat{\bf D}$, and $\hat{\bf A}$ in Ref. [1]). This step uses the output from Step 1 (pos_link_step1.mat).
Specify the region of interest (ROI) for the analysis.
Make sure that the image region does not extend beyond the boundaries in the first and last frames. See the sample `input/xyROI.mat` for setting the ROI.
In the sample, the side length of a ROI is set to $L = 100$.
When the calculation is finished, `pos_link_step2.mat` will be output.

* **Step 3: Visualization of Data**  
Load `pos_link_step2.mat` and generate figures for visualization. The resulting figures will be saved in the `Figures` folder.

## References

[1] Namba T, Sugimura K, Ishihara S. In vivo validation of kinematic relationships for epithelial morphogenesis. XXXX. 2025;XX. doi:......

[2] Guirao B, et al. Unified quantitative characterization of epithelial tissue development. Elife, Vol. 4, p. e08519,2015. doi:10.7554/eLife.08519

[3] Ishihara S, Marcq P, Sugimura K. From cells to tissue: A continuum model of epithelial mechanics. Phys. Rev. E96, 022418, 2017. doi:10.1103/PhysRevE.96.022418

[4] Stamhuis E, Thielicke W. PIVlab – Towards User-friendly, Affordable and Accurate Digital Particle Image Velocimetry in MATLAB. Journal of Open Research Software. 2014;2(1). doi:10.5334/jors.bl.

