from mayavi import mlab
import argparse
import sys

from EMS.utilities import read_ply, showPoints
from EMS.EMS_recovery import EMS_recovery

import timeit

def main(argv):

    parser = argparse.ArgumentParser(
        description='Probabilistic Recovery of a superquadric surface from a point cloud file *.ply.')

    parser.add_argument(
        'path_to_data',
        # default = '~/EMS-probabilistic_superquadric_fitting/MATLAB/example_scripts/data/noisy_pointCloud_example_1.ply',
        help='Path to the directory containing the point cloud file *.ply.'
    )
    parser.add_argument(
        '--visualize',
        action = 'store_true',
        help='Visualize the recoverd superquadric and the input point cloud.'
    )
    parser.add_argument(
        '--runtime',
        action = 'store_true',
        help='Show the runtime.'
    )
    parser.add_argument(
        '--result',
        action = 'store_true',       
        help='Print the recovered superquadric parameter.'
    )

    args = parser.parse_args(argv)
    
    print('----------------------------------------------------')
    print('Loading point cloud from: ', args.path_to_data, '...')
    point = read_ply(args.path_to_data)
    print('Point cloud loaded.')
    print('----------------------------------------------------')

    # first run to eliminate jit compiling time
    sq_recovered, p = EMS_recovery(point)

    start = timeit.default_timer()
    sq_recovered, p = EMS_recovery(point, OutlierRatio=0.2)
    stop = timeit.default_timer()
    print('Superquadric Recovered.')
    if args.runtime is True:
        print('Runtime: ', (stop - start) * 1000, 'ms')
    print('----------------------------------------------------')
    
    if args.result is True:
        print('shape =', sq_recovered.shape)
        print('scale =', sq_recovered.scale)
        print('euler =', sq_recovered.euler)
        print('translation =', sq_recovered.translation)
        print('----------------------------------------------------')
    
    if args.visualize is True:
        fig = mlab.figure(size=(400, 400), bgcolor=(1, 1, 1))
        sq_recovered.showSuperquadric(arclength = 0.2)
        showPoints(point)
        mlab.show()


if __name__ == "__main__":
    main(sys.argv[1:])