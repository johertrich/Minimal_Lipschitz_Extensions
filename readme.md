# Minimal Lipschitz and ∞-harmonic extensions of vector-valued functions on finite graphs

This code belongs to the papers [1] and [2]. Please cite the paper if you use this code. 

They are available at  
https://doi.org/10.1093/imaiai/iaz033 and https://doi.org/10.1007/978-3-030-22368-7_15  
Further an Arxiv-preprint of [1] can be found at  
https://arxiv.org/abs/1903.04873

The code in this repository reproduces the examples from the paper [1]. Note that we use in the directory `MexNonlocal` some functions from the implementation of [3]. For questions and bug reports, please contact Johannes Hertrich (j.hertrich(at)math.tu-berlin.de).

## Requirements and Usage

All examples were tested using Matlab 2019a. Each example reproduces the images for one of the
figures in the paper. Some of the examples have a very long runtime! The code is highly
experimental and far from optimized (especially the building methods for the nonlocal graphs),
so use it with care.
Probably, the mex files have to be recompiled for executing Example_Singapur.m. Necessary
parameter changes for reproducing all images are indicated at the top of each script.

## Reference

[1] M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.  
Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.  
Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.

[2] J. Hertrich, M. Bačák, S. Neumayer, G. Steidl.  
Minimal Lipschitz extensions for vector-valued functions on finite graphs.  
M. Burger, J. Lellmann and J. Modersitzki (eds.)  
Scale Space and Variational Methods in Computer Vision.  
Lecture Notes in Computer Science, 11603, 183-195, 2019.  

[3] F. Laus, F. Pierre, and G. Steidl.  
Nonlocal myriad filters for Cauchy noise removal.  
Journal of Mathematical Imaging and Vision 60.8 (2018): 1324-1354.
