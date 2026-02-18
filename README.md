# radial-ssfp-2d

_Reconstruction of tiny golden-angle radial 2D balanced steady state free precession MRI with motion-correction and cardiac metric optimized gating._

## Framework

1. **Real-Time Image Reconstruction** - Low temporal resolution frames for motion analysis
2. **Motion Estimation** - Translational motion estimation from real-time data  
3. **Motion Correction** - Application of estimated motion to raw k-space data
4. **Cardiac Synchronization** - Metric-optimized gating for cardiac phase detection
5. **Cine Reconstruction** - High-quality motion-corrected cine image reconstruction