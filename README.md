# radial-ssfp-2d

_Reconstruction of tiny golden-angle radial 2D balanced steady state free precession MRI with motion-correction and cardiac metric optimized gating._

## Framework

1. **Real-Time Image Reconstruction** - Low temporal resolution frames for motion analysis
2. **Motion Estimation** - Translational motion estimation from real-time data  
3. **Motion Correction** - Application of estimated motion to raw k-space data
4. **Cardiac Synchronization** - Metric-optimized gating for cardiac phase detection
5. **Cine Reconstruction** - High-quality motion-corrected cine image reconstruction

## Requirements

Current implementation requires:
- **MATLAB** - including Signal Processing, Image Processing, Statistics and Machine Learning Toolbox, and Curve toolboxes
- **MCNUFFT Class and NUFFT Toolbox** - from [GRASP ](https://cai2r.net/resources/grasp-matlab-code/) (Li Feng & Ricardo Otazo) and [MIRT](https://github.com/JeffFessler/mirt) (Jeff Fessler), respectively
- **Mutual Information Toolbox** - from [Matlab Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/14888-mutual-information-computation) (Hanchuan Peng)
- **Siemens Raw Data Loading Functions** - mapVBVD.m, available from Siemens user forum and elsewhere
- **Coil Sensitivity Map Functions** - from [ISMRM course](https://github.com/hansenms/ismrm_sunrise_matlab) (Michael Hansen & Philip Beatty)

## Usage

See 'script_radial_ssfp_2d_recon_moco_cine.m'.
