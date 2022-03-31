## Guidelines for the MATLAB Implementation

This is the guideline and structual explanation of the MATLAB implementation of the EMS algorithm.

### Installation

The code is written and tested on MATALB R2021a.
Should be compatible with releases newer than R2017.
The EMS algorithm is self-contained.
To run hierarchical-EMS and demo scripts, the Computer Vision Toolbox is required.
For compaitibility problems, you are more than welcomed to put up an issue.

### File Structure

- The source code is located in `/src`, where `EMS.m` is the core algorithm for single superquadric recovery, and `Hierarchical_EMS.m` is the algorithm for multisuperquadric recovery.
- Utility functions for demonstrations (*e.g.* superquadric rendering, sampling, baseline algorithms) are located in `/src/utilities`.
- Please add `/src` to path before excuting the test scripts.
- Test scripts are located in `/example_scripts`.
- Demo point cloud files (`.ply`) are provided in `/example_scripts/data`.


### Run Demo

#### Demo for single superquadric recovery

There are two test scripts you can play with: `/example_scripts/example_single_superquadric_script.m` and `/example_scripts/random_example_script`.
- With `example_single_superquadric_script.m`, you will load a `.ply` in `/example_scripts/data/single_superquadric` and recover the superquadric representation. The script will compare and visualize the results obtained by the baseline methods (Radial-LSQ, Numerical Stable) with ours.
- With `/example_scripts/random_example_script`, a point cloud of a random superquadric will be generated, and then recovered. In the script, you can specify the partial, outlier and noise levels.

#### Demo for multiple superquadric recovery

For multiple superquadric recovery, you can run the demo script `/example_scripts/multi_superquadric_script.m`.
You can load point clouds from `/example_scripts/data/multi_superquadrics`.
