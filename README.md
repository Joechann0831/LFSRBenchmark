# Light Field Super-Resolution: A Benchmark
A collection of codes and datasets for light field super-resolution methods evaluated in the following paper </br>
Zhen Cheng, Zhiwei Xiong, Chang Chen, Dong Liu. [Light Field Super-Resolution: A Benchmark](http://openaccess.thecvf.com/content_CVPRW_2019/html/NTIRE/Cheng_Light_Field_Super-Resolution_A_Benchmark_CVPRW_2019_paper.html). In CVPRW 2019. </br>

## Datasets
### HCI
It is a dataset synthesized by graphic software [1]. Our used data is an early version of this dataset, which cannot be found in their website now. We make it available again at </br>
http://pan.bitahub.com/index.php?mod=shares&sid=eTJ2bFFQR3BzTm5FTGc5OHJBWW1zWnc4eHpMYWdaeXE5WVdwSmc

### EPFL
It is a real-world dataset captured by Lytro Illum camera [2]. Due to the vignetting effect, we make a further rectification after calibration with the [Light Field Toolbox for MATLAB](http://dgd.vision/Tools/LFToolbox/). Detail information about the rectification can be found in our supplementary document. We make this rectified dataset available at </br>
http://pan.bitahub.com/index.php?mod=shares&sid=eTJ2bFFQR3BzTm5FTGc5OHJBWW1zWnc4eHpMYWdaeXE5WVdwSmc

## Codes

We select four representative light field SR methods from three categories. Among them, GB [3] and RR [4] are implemented with the official codes, while PRO [5], LFCNN [6] and VDSR [7] are re-implemented by ourselves.

## GB

GB is a graph-based framework for light field super-resolution. The official code can be found at [GB-official-code](https://github.com/rossimattia/light-field-super-resolution).

### Requirements

- MATLAB

### Run

You can simply run the code in the MATLAB command window with the following commands:

```matlab
cd GB
graph_based_SR_on_EPFL
```

### Tips for running
You can increase the parameter *poolSize* in the script graph_based_SR.m to accelerate the running. If your CPU resource is limited, you can decrease the pool size as well.

For more details about the code, please refer to the [official website](https://github.com/rossimattia/light-field-super-resolution).

## RR

RR is a learning-based light field SR method with Principal Component Analysis (PCA) and Ridge Regression (RR). The official code can be found at [RR-official-code](https://github.com/rrfarr/LF-Editing).

### Requirements

- MATLAB

### Run

You can simply run the code in the MATLAB command window with the following commands:

```matlab
cd RR
pca_rr_bm_for_EPFL
```

### Tips for running
The above script can also be executed with the Run button in MATLAB edit mode. For more details, please refer to the comments inside the codes as well as the README file in the [official website](https://github.com/rrfarr/LF-Editing).

## PRO

PRO is a projection-based light field super-resolution method. The method is originated from [5] and implemented by ourselves.

### Requirements

- MATLAB

### Run

You can simply run the code in the MATLAB command window with the following commands:

```matlab
cd PRO
projection_based_SR_on_EPFL_scale2
```

### Tips for running
We recommend to run this script on Windows. In Linux OS such as Ubuntu, the script will get stuck due to some unknown configuration problems. The above script is used for scale 2, if you want to upsample the light field with scale 3, you should change the function *sr_projection_scale2* to *sr_projection_scale3*. Note that the input parameters of these two functions are not exactly the same, please pay some attention.

For the details of the implementation, please refer to the code or our supplementary document.

## LFCNN

LFCNN is the first CNN-based method which is published in ICCVW2015. The original LFCNN uses SRCNN structure as their backbone model. We upgrade the shallow SRCNN structure to the deep VDSR structure, which promotes LFCNN's performance for a fair comparison with VDSR [7].

### Requirements and dependencies

- Caffe 1.0
- CUDA and Cudnn suited for Caffe 1.0
- MATLAB with pre-compiled matcaffe

### Train

Please use the code under folder *data_generation* for test and training data generation at first. Then please change the directory for snapshots, logs, data path, and the executable caffe tool in corresponding files in the directory *training_configs*.

### Test

The test code is *LFCNN/test_code/test_LFCNN_to_whole_LF.m*, please change the folder path of matcaffe to your path and run the script in MATLAB command window.

### Tips for running

Note that we used two relatively small datasets for training and testing, so we used Cross-Validation strategy. With more training data, LFCNN can achieve better performance. Note that the given model is trained for testing *scenes 1,2 and 3* in EPFL dataset under the assumption of Bicubic downsampling and a scale factor of 2.

## VDSR

VDSR is a representative CNN-based single image SR method without using any angular information. We train the network using the same training set as in [7]. The official code is implemented with MatConvNet, we re-implement it with Caffe 1.0.

### Requirements and dependencies

- Caffe 1.0
- CUDA and Cudnn suited for Caffe 1.0
- MATLAB with pre-compiled matcaffe

### Train

Please use the code under folder *data_generation* for test and training data generation at first. Then please change the directory for snapshots, logs, data path and the executable caffe tool in corresponding files in the directory *training_configs*.

### Test

The test code is *VDSR/test_code/test_VDSR_whole_LF.m*, please change the folder path of matcaffe to your path and run the script in MATLAB command window. The given model is trained under the assumption of Bicubic downsampling and a scale factor of 2.

## Reference

[1] S. Wanner, S. Meister, and B. Goldlucke. Datasets and benchmarks for densely sampled 4d light fields. In International Symposium on Vision Modeling and Visualization, 2013.

[2] M. Rerabek and T. Ebrahimi. New light field image dataset. In International Conference on Quality of Multimedia Experience (QoMEX), 2016.

[3] M. Rossi and P. Frossard. Graph-based light field super-resolution. In MMSP, 2017.

[4] R. A. Farrugia, C. Galea, and C. Guillemot. Super resolution of light field images using linear subspace projection of patch-volume. IEEE Journal of Selected Topics in Signal Processing, 11(7):1058-1071, 2017.

[5] C.-K. Liang and R. Ramamoorthi. A light transport framework for lenslet light field cameras. ACM Transactions on Graphics, 34(2):16:1-16:19, 2015.

[6] Y. Yoon, H. G. Jeon, D. Yoo, J. Y. Lee, and I. S. Kweon. Learning a deep convolutional network for light-field image super-resolution. In ICCVW, 2015.

[7] J. Kim, J. Kwon Lee, and K. Mu Lee. Accurate image super-resolution using very deep convolutional networks. In CVPR, 2016.

The citation of our benchmark paper:

```latex
@InProceedings{Cheng_2019_CVPR_Workshops,
author = {Cheng, Zhen and Xiong, Zhiwei and Chen, Chang and Liu, Dong},
title = {Light Field Super-Resolution: A Benchmark},
booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR) Workshops},
month = {June},
year = {2019}
}
```





