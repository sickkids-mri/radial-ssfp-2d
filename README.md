# radial-ssfp-2d

_reconstruction of motion-corrected tiny golden-angle radial 2D balanced steady state free precession MRI_

## Framework

Framework for reconstructing images of the fetal heart in the presence of motion using (a) a golden angle radial acquisition to reconstruct real-time images to (b) automatically assess and reject through-plane motion and (c) estimate in-plane motion and fetal heart rate. Then (d) k-space is corrected and sorted, and (e) cine images are reconstructed. 

![](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs12968-017-0346-6/MediaObjects/12968_2017_346_Fig1_HTML.gif?as=webp)  
–source: [Roy _et al._ JCMR 2017 Fig. 1](https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0346-6/figures/1). 

## Citation

Please cite the following in publications that use this work

Metric Optimized Gating (MOG) for Cine Imaging of the Fetal Heart:

    Roy CW, Seed M, van Amerom JFP, Al Nafisi B, Grosse-Wortmann L, Yoo S-J, Macgowan CK. __Dynamic imaging of the fetal heart using metric optimized gating.__ Magn Reson Med. 2013;70(6):1598–607. doi: [10.1002/mrm.24614](https://doi.org/10.1002/mrm.24614).

Motion-Compensated Reconstruction:

    Roy CW, Seed M, Kingdom JC, Macgowan CK. __Motion compensated cine CMR of the fetal heart using radial undersampling and compressed sensing.__ J. Cardiovasc. Magn. Reson. 2017;19:29. doi: [10.1186/s12968-017-0346-6](https://doi.org/10.1186/s12968-017-0346-6).

Multi-Slice Evaluation: 

    Roy CW, Marini D, Lloyd DFA, Mawad W, Yoo S-J, Schrauben EM, Jaeggi E, Seed M, Macgowan CK. __Preliminary Experience Using Motion Compensated CINE Magnetic Resonance Imaging to Visualise Fetal Congenital Heart Disease.__ Circ. Cardiovasc. Imaging 2018;11. doi: [10.1161/CIRCIMAGING.118.007745](https://doi.org/10.1161/CIRCIMAGING.118.007745).

## Requirements

Assumes data in Siemens .dat VB/VD MRI raw data format, collected using 2D Tiny Golden Angle Radial bSSFP (available via [Siemens C<sup>2</sup>P Exchange](https://webclient.ca.api.teamplay.siemens.com/#/c2p)).  

Current implementation requires:
- **MATLAB** - including Signal Processing, Image Processing, Statistics and Machine Learning, and Curve Fitting toolboxes
- **MCNUFFT Class and NUFFT Toolbox** - from [GRASP ](https://cai2r.net/resources/grasp-matlab-code/) (Li Feng & Ricardo Otazo) and [MIRT](https://github.com/JeffFessler/mirt) (Jeff Fessler), respectively
- **Mutual Information Toolbox** - from [Matlab Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/14888-mutual-information-computation) (Hanchuan Peng)
- **Siemens Raw Data Loading Functions** - mapVBVD.m, available from Siemens and elsewhere
- **Coil Sensitivity Map Functions** - from [ISMRM course](https://github.com/hansenms/ismrm_sunrise_matlab) (Michael Hansen & Philip Beatty)

Contact authors for assistance with non-Mathworks toolboxes. 

## Installation 

Clone repo to reconstruction system.

## Usage

See 'script_radial_ssfp_2d_recon_moco_cine.m'.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License, Terms and Conditions

See the [LICENSE](LICENSE) file for license rights and limitations (MIT).

Use is also subject to Terms & Conditions associated with the 2D Tiny Golden Angle Radial bSSFP MRI acquisition sequence available on the [Siemens C<sup>2</sup>P Exchange](https://webclient.ca.api.teamplay.siemens.com/#/c2p)).  
