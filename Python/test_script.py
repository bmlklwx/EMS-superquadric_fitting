import numpy as np
from mayavi import mlab
from EMS.utilities import read_ply, showPoints
from EMS.EMS_recovery import EMS_recovery

import timeit

if __name__ == '__main__':
    # sq = superquadric([1, 1], [1, 2, 1], [np.pi/3, np.pi/2, 0], [0, 1, 0])
    # sq.showSuperquadric(arclength = 0.1)
    # print(sq.euler)
    # print(sq.RotM)
    # print(sq.quat)
    
    # path_to_file = '/home/saintsbury/src/EMS-probabilistic_superquadric_fitting/MATLAB/example_scripts/data/noisy_pointCloud_example_1.ply'
    # path_to_file = '/home/saintsbury/src/EMS-probabilistic_superquadric_fitting/MATLAB/example_scripts/data/partial_pointCloud_example_2.ply'
    path_to_file = '/home/saintsbury/src/EMS-probabilistic_superquadric_fitting/MATLAB/example_scripts/data/000.ply'
    point = read_ply(path_to_file)
    # showPoints(point)

    start = timeit.default_timer()
    sq_recovered, p = EMS_recovery(point)
    stop = timeit.default_timer()
    print('Time (jit compile included): ', stop - start, 's')

    start = timeit.default_timer()
    sq_recovered, p = EMS_recovery(point, OutlierRatio=0.2)
    stop = timeit.default_timer()
    print('Time (runtime): ', stop - start, 's')
    
    print('shape =', sq_recovered.shape)
    print('scale =', sq_recovered.scale)
    print('euler =', sq_recovered.euler)
    print('translation =', sq_recovered.translation)

    fig = mlab.figure(size=(400, 400), bgcolor=(1, 1, 1))
    sq_recovered.showSuperquadric(arclength = 0.2)
    showPoints(point)
    mlab.show()