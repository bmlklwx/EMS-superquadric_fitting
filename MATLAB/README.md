## Guidelines for the MATLAB Implementation

This is the guideline and structual explanation of the MATLAB implementation of the EMS algorithm.

### Installation

The code is written and tested on MATALB R2021a.
Should be compatible with releases newer than R2017.
For other compaitibility problems, you are more than welcomed to put up an issue.

### File Structure

1. The source code is located in `/src`, where `EMS.m` is the core algorithm for single superquadric recovery, and `Hierarchical_EMS.m` is the algorithm for multisuperquadric recovery.
2. Utility functions for demonstrations (*e.g.* superquadric rendering, sampling, baseline algorithms) are located in `/src/utilities`.
3. Test scripts are located in `/example_scripts`.
4. A few point cloud files (`.ply`) are provided in `/example_scripts/data`.


### Run Demo

#### Demo for single superquadric recovery

#### Demo for multiple superquadric recovery
