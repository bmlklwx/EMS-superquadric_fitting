# EMS: A probabilistic approach for superquadric recovery

This repo provides the source code for the CVPR2022 paper:

> [**Robust and Accurate Superquadric Recovery: a Probabilistic Approach**](https://arxiv.org/abs/2111.14517 "ArXiv version of the paper.")  
> Weixiao Liu, Yuwei Wu, [Sipu Ruan](https://ruansp.github.io/), [Gregory S. Chirikjian](https://cde.nus.edu.sg/me/staff/chirikjian-gregory-s/)

We propose an algorithm to recover a superquadric surface/primitive from a given point cloud, with good robustness, accuracy and efficiency.
The superquadric recovered from a point cloud provides a concise, volumetric, and geometrically meaningful interpretation of objects and environment. It can work as a low-level volumetric representation, from which higher level tasks, *e.g.*, motion planning, collision detection and robot-environment interaction, can be built up.

<img src="/figures/Superquadrics.png" alt="superquadrics" width="600"/>

## Abstract

Interpreting objects with basic geometric primitives has long been studied in computer vision. 
Among geometric primitives, superquadrics are well known for their simple implicit expressions and capability of representing a wide range of shapes with few parameters. 
However, as the first and foremost step, recovering superquadrics accurately and robustly from 3D data still remains challenging. 
The existing methods are subject to local optima and are sensitive to noise and outliers in real-world scenarios, resulting in frequent failure in capturing geometric shapes. 
In this paper, we propose the first probabilistic method to recover superquadrics from point clouds. 
Our method builds a Gaussian-uniform mixture model (GUM) on the parametric surface of a superquadric, which explicitly models the generation of outliers and noise.
The superquadric recovery is formulated as a Maximum Likelihood Estimation (MLE) problem. 
We propose an algorithm, Expectation, Maximization, and Switching (EMS), to solve this problem, where: (1) outliers are predicted from the posterior perspective; (2) the superquadric parameter is optimized by the trust-region reflective algorithm; and (3) local optima are avoided by globally searching and switching among parameters encoding similar superquadrics. 
We show that our method can be extended to the multi-superquadrics recovery for complex objects. 
The proposed method outperforms the state-of-the-art in terms of accuracy, efficiency, and robustness on both synthetic and real-world datasets.

## Implementations

This repo provides two implementations in both [MATALB](/MATLAB) and [Python](/Python). 
The algrithm was first implemented in MATALB and later rewritten in Python for a broader audience.
Detailed guidelines for dependency and installation please refer to the linked READMEs (for [MATLAB](/MATLAB) and for [Python](/Python)).
~~The demo for the multiple superquadric recovery is only available in MATLAB.~~
The demo for the multiple superquadric recovery is now available in Python as well.
Please refer to [multiquadric_test.py](/Python/tests/multiquadric_test.py).
Thanks @[stanzwinkels](https://github.com/stanzwinkels) for his contribution.
The Python implementation shows better efficiency (about 3-5 times faster), with the help of NUMBA (a JIT compiler). 
C++ version is planned.

## Superquadrics Sampling, visiualization and other utilities functions

This repo also provide several useful utility functions related to superquadrics.

In the [MATLAB](/MATLAB) implementation, the utilities are located in [/MATLAB/src/utilities](/MATLAB/src/utilities), where
* [superquadricsFitting.m](/MATLAB/src/utilities/superquadricsFitting.m) summarized the baseline least square superquadric recovery methods, based on different objective functions (implicit function, radial distance ...).
* [numerical_fitting.m](/MATLAB/src/utilities/numerical_fitting.m) is an implementation of the numerical stable recovery method proposed in N. Vaskevicius and A. Birk, "Revisiting Superquadric Fitting: A Numerically Stable Formulation," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 41, no. 1, pp. 220-233, 1 Jan. 2019.
* [showSuperquadrics.m](/MATLAB/src/utilities/showSuperquadrics.m) is for visualization of superquadrics. Note that one can choose to visualize a tapered superquadric in its option.
* [sphericalProduct_sampling.m](/MATLAB/src/utilities/sphericalProduct_sampling.m) is an algorithm to sample points almost uniformly spaced on the surface of a given superquadric.

In the [Python](/Python) implementation, the functions above are summarized in [/Python/src/EMS/utilities.py](/Python/src/EMS/utilities.py).


## Related Works
To be added ...
