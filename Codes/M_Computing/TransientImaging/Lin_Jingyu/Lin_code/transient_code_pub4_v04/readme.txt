Fourier Domain Transient Imaging, version 0.4
================================================

Refer to our papers for technical details: 

Jingyu Lin, Yebin Liu, M. Hullin, Qionghai Dai. 
Fourier Analysis on Transient Imaging by Multifrequency Time-of-Flight Camera. 
IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2014, pp. 3230-3237.

Jingyu Lin, Yebin Liu, Jinli Suo, Qionghai Dai. 
Fourier Domain Transient Imaging. 
Transactions on Pattern Analysis and Machine Intelligence, 2014, submitted.


Copyright (C) 2013-2014 Jinyu Lin, linjy02@hotmail.com

================================================
Scheme
------------------------------------------------
1. Calibration
Invoke commands in 'calib' to capture data and fit the correlation function.
Correlation functions are in 'calib/cfs'.

2. Acquisition 
Invoke commands in 'data' to capture data.
Data are saved in 'data/rawdata'.

3. Processing
Invoke commands in 'cmpu_transient' to reconstruct transient images.
Invoke commands in 'cmpu_depth' to compute depth maps.
Invoke 'rec_ubcdata' in 'cmpu_transient_ubc' to compute transient images from UBC data


================================================
Path information
------------------------------------------------
calib - commands for camera calibration, acquired data in 'calib_xxx', and fitted correlation functions in 'cfs'.

data - commands for data acquisition, and acquired data in 'rawdata'.
data_ubc - path for ubc data.

cmpu_transient - commands for the reconstruction of transient images.
cmpu_transient_ubc - reconstruct transient images from ubc data.

cmpu_depth - commands for the computation of depth maps.

functions - supporting functions

TOF_tools - exe files for camera acquisition.
