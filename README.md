# EMS: A probabilistic approach for superquadric recovery

This repo provides the source code for the CVPR2022 paper:

> [**Robust and Accurate Superquadric Recovery: a Probabilistic Approach**](https://arxiv.org/abs/2111.14517 > "ArXiv version of the paper.")  
> Weixiao Liu, Yuwei Wu, [Sipu Ruan](https://ruansp.github.io/), [Gregory S. Chirikjian](https://cde.nus.edu.sg/me/staff/chirikjian-gregory-s/)

We propose an algorithm to recover a superquadric surface/primitive from a given point cloud, with good robustness, accuracy and efficiency.

<img src="/figures/Superquadrics.png" alt="superquadrics" width="600"/>

## Abstration

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

## Dependency and Installation
