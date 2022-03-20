import numpy as np
from numba import njit
from scipy.optimize import least_squares

from EMS.superquadrics import rotations, superquadric

def EMS_recovery(
        point, OutlierRatio=0.1, MaxIterationEM=20,
        ToleranceEM=1e-3, RelativeToleranceEM=1e-1,
        MaxOptiIterations=3, Sigma=0, MaxiSwitch=2,
        AdaptiveUpperBound=False, Rescale=True):
    # The function conducting probabilistic superquadric recovery.
    # Input: point - point cloud np array of N * 3
    #

    # ---------------------------------------INITIALIZATIONS--------------------------------------------
    # translate the points to the center of mass
    point = np.array(point, dtype=float)
    t0 = np.mean(point, 0)
    point = point - t0

    # rescale
    if Rescale is True:
        max_length = np.max(point)
        scale = max_length / 10
        point = point / scale

    # eigen analysis for rotation initialization
    EigVec = EigenAnalysis(point)
    R0 = rotations()
    R0.RotM = np.array([-EigVec[:, 0], -EigVec[:, 2],
                       np.cross(EigVec[:, 0], EigVec[:, 2])]).T
    euler0 = R0.euler

    # scale initialization
    point_rot0 = point @ R0.RotM
    s0 = np.median(np.abs(point_rot0), 0)

    # initialize configuration
    x0 = np.array([1.0, 1.0, s0[0], s0[1], s0[2],
                  euler0[0], euler0[1], euler0[2], 0, 0, 0])

    # set lower and upper bounds for the superquadrics
    upper = 4 * np.max(np.abs(point))
    lb = np.array([0, 0, 0.001, 0.001, 0.001, -2 * np.pi, -2 *
                  np.pi, -2 * np.pi, -upper, -upper, -upper])
    ub = np.array([2.0, 2.0, upper, upper, upper, 2 * np.pi,
                  2 * np.pi, 2 * np.pi, upper, upper, upper])

    # calculate bounding volume of ourlier space
    V = BoundVolume(point_rot0)

    # set prior outlier density
    p0 = 1 / V

    # initialize variance
    if Sigma == 0:
        sigma2 = V ** (1 / 3) / 10
    else:
        sigma2 = Sigma

    # initialize EMS
    x = x0
    cost = 0.0
    num_switch = np.int(0)
    p = np.ones(point.shape[0])

    # ---------------------------------------EMS ALGORITHM--------------------------------------------
    for iterEM in range(MaxIterationEM):
        # evaluating distance from points to superquadric
        dist = Distance(point, x)

        # inferring the postierior outlier probability (E-step)
        if OutlierRatio != 0:
            p = OutlierProb(dist, sigma2, OutlierRatio, p0)

        # calculate adaptive upper bound
        if AdaptiveUpperBound is True:
            R_cur = Euler2RotM(x[5: 8])
            point_cur = point @ R_cur - x[8: 11] @ R_cur
            ub_a = 1.1 * np.max(np.abs(point_cur), 0)
            ub[2: 5] = ub_a
            ub[8: 11] = ub_a
            lb[8: 11] = -ub_a

        # Optimize the superquadric configuration (M-step)
        optfunc = least_squares(CostFunc, x, bounds=(
            lb, ub), max_nfev=MaxOptiIterations, args=(point, p, sigma2))
        x_n = optfunc.x
        cost_n = 2 * optfunc.cost

        # update sigma
        sigma2_n = cost_n / (3 * np.sum(p))       

        # evaluate raletive decreasing of cost
        relative_cost = (cost - cost_n) / cost_n

        # check optimality for termination
        if (cost_n < ToleranceEM and iterEM > 0) or \
                (relative_cost < RelativeToleranceEM and num_switch >= MaxiSwitch and iterEM > 4):
            x = x_n
            break

        # check for entering similarity switch
        if relative_cost < RelativeToleranceEM and iterEM > 0:
            # entering similarity switch (S-step)
            # initialize swith success flag
            switch_success = False

            # search for similarity candidates
            x_candidate = SimilarityCandidates(x)

            # evaluating switch (S-step)
            x, cost, sigma2, switch_success = Switch(
                x_candidate, point, p, AdaptiveUpperBound, ub, lb, MaxOptiIterations, \
                sigma2, sigma2_n, cost, cost_n, x_n, switch_success
            )

            num_switch = num_switch + 1

        else:
            # update parameter and prepare for the next EM iteration
            cost = cost_n
            sigma2 = sigma2_n
            x = x_n
    
    if Rescale is True:
        x[2 : 5] = x[2 : 5] * scale
        x[8 : 11] = x[8 : 11] * scale
    
    x[8 : 11] = x[8 : 11] + t0

    sq = superquadric(x[0 : 2], x[2 : 5], x[5 : 8], x[8 : 11])

    return sq, p

# ---------------------------------------UTILITIES-------------------------------------------
@njit(cache=True)
def SimilarityCandidates(x):
    # axis mismatch similarity
    axis_0 = Euler2RotM(x[5: 8])
    axis_1 = axis_0[:, np.array([1, 2, 0])]
    axis_2 = axis_0[:, np.array([2, 0, 1])]
    eul_1 = RotM2Euler(axis_1)
    eul_2 = RotM2Euler(axis_2)
    x_axis = np.array(
        [[x[1], x[0], x[3], x[4], x[2], eul_1[0], eul_1[1], eul_1[2], x[8], x[9], x[10]],
         [x[1], x[0], x[4], x[2], x[3], eul_2[0], eul_2[1], eul_2[2], x[8], x[9], x[10]]]
    )   

    # duality similarities
    scale_ratio = x[np.array([3, 4, 2])] / x[2 : 5]
    scale_idx = np.argwhere(np.logical_and(scale_ratio > 0.6, scale_ratio < 1.4))
    x_rot = np.zeros((scale_idx.shape[0], 11))
    
    for idx in range(scale_idx.shape[0]):
        if scale_idx[idx, 0] == 0:
            eul_rot = RotM2Euler(axis_0 @ Euler2RotM(np.array([np.pi / 4, 0.0, 0.0])))
            if x[1] <= 1:
                x_rot[idx, :] = np.array(
                    [x[0], 2 - x[1], 
                    ((1 - np.sqrt(2)) * x[1] + np.sqrt(2)) * min(x[2], x[3]),
                    ((1 - np.sqrt(2)) * x[1] + np.sqrt(2)) * min(x[2], x[3]),
                    x[4], eul_rot[0], eul_rot[1], eul_rot[2],
                    x[8], x[9], x[10]]
                )
            else:
                x_rot[idx, :] = np.array(
                    [x[0], 2 - x[1], 
                    ((np.sqrt(2)/2 - 1) * x[1] + 2 - np.sqrt(2) / 2) * min(x[2], x[3]),
                    ((np.sqrt(2)/2 - 1) * x[1] + 2 - np.sqrt(2) / 2) * min(x[2], x[3]),
                    x[4], eul_rot[0], eul_rot[1], eul_rot[2],
                    x[8], x[9], x[10]]
                )

        elif scale_idx[idx, 0] == 1:
            eul_rot = RotM2Euler(axis_1 @ Euler2RotM(np.array([np.pi / 4, 0.0, 0.0])))
            if x[0] <= 1:
                x_rot[idx, :] = np.array(
                    [x[1], 2 - x[0], 
                    ((1 - np.sqrt(2)) * x[0] + np.sqrt(2)) * min(x[3], x[4]),
                    ((1 - np.sqrt(2)) * x[0] + np.sqrt(2)) * min(x[3], x[4]),
                    x[2], eul_rot[0], eul_rot[1], eul_rot[2],
                    x[8], x[9], x[10]]
                )
            else:
                x_rot[idx, :] = np.array(
                    [x[1], 2 - x[0], 
                    ((np.sqrt(2)/2 - 1) * x[0] + 2 - np.sqrt(2)/2) * min(x[3], x[4]),
                    ((np.sqrt(2)/2 - 1) * x[0] + 2 - np.sqrt(2)/2) * min(x[3], x[4]),
                    x[2], eul_rot[0], eul_rot[1], eul_rot[2],
                    x[8], x[9], x[10]]
                )

        elif scale_idx[idx, 0] == 2:
            eul_rot = RotM2Euler(axis_2 @ Euler2RotM(np.array([np.pi / 4, 0.0, 0.0])))
            if x[0] <= 1:
                x_rot[idx, :] = np.array(
                    [x[1], 2 - x[0], 
                    ((1 - np.sqrt(2)) * x[0] + np.sqrt(2)) * min(x[4], x[2]),
                    ((1 - np.sqrt(2)) * x[0] + np.sqrt(2)) * min(x[4], x[2]),
                    x[2], eul_rot[0], eul_rot[1], eul_rot[2],
                    x[8], x[9], x[10]]
                )
            else:
                x_rot[idx, :] = np.array(
                    [x[1], 2 - x[0], 
                    ((np.sqrt(2)/2 - 1) * x[0] + 2 - np.sqrt(2)/2) * min(x[4], x[2]),
                    ((np.sqrt(2)/2 - 1) * x[0] + 2 - np.sqrt(2)/2) * min(x[4], x[2]),
                    x[2], eul_rot[0], eul_rot[1], eul_rot[2],
                    x[8], x[9], x[10]]
                )
    
    x_candidate = np.zeros((2 + x_rot.shape[0], 11))
    x_candidate[0 : 2] = x_axis
    if scale_idx.shape[0] > 0:
        x_candidate[2 : 2 + scale_idx.shape[0]] = x_rot
    
    return x_candidate

def Switch(
    x_candidate, point, p, AdaptiveUpperBound, ub, lb, MaxOptiIterations, \
        sigma2, sigma2_n, cost, cost_n, x_n, switch_success
):
    cost_candidate = SwitchCost(x_candidate, point, p)
    idx_nan = np.argwhere(
        np.logical_and(~np.isnan(cost_candidate), ~np.isinf(cost_candidate))
    ).reshape(1, -1)[0]

    cost_candidate = cost_candidate[idx_nan]
    idx = np.argsort(cost_candidate)

    for i in idx:
        if AdaptiveUpperBound is True:
            R_cur = Euler2RotM(x_candidate[i, 5: 8])
            point_cur = point @ R_cur - x_candidate[i, 8: 11] @ R_cur
            ub_a = 1.1 * np.max(np.abs(point_cur), 0)
            ub[2: 5] = ub_a
            ub[8: 11] = ub_a
            lb[8: 11] = -ub_a

        x_candidate[i] = np.minimum(x_candidate[i], ub)
        x_candidate[i] = np.maximum(x_candidate[i], lb)

        optfunc = least_squares(CostFunc, x_candidate[i], bounds=(
                lb, ub), max_nfev=MaxOptiIterations, args=(point, p, sigma2))
        x_switch = optfunc.x
        cost_switch = 2 * optfunc.cost

        if cost_switch < min(cost_n, cost):
            x = x_switch
            cost = cost_switch

            sigma2 = cost_switch / (3 * sum(p))
            switch_success = True
            break
    
    if switch_success == False:
        cost = cost_n
        sigma2 = sigma2_n
        x = x_n

    return x, cost, sigma2, switch_success
    
@njit(cache=True)
def SwitchCost(x_candidate, point, p):
    val = np.zeros(x_candidate.shape[0])
    for i in range(x_candidate.shape[0]):
        val[i] = np.sum(p * (Distance(point, x_candidate[i]) ** 2))
    return val

@njit(cache=True)
def EigenAnalysis(point):
    CovM = point.T @ point / point.shape[0]
    EVal, EVec = np.linalg.eig(CovM)
    idx = np.flip(np.argsort(EVal))
    return EVec[:, idx]

@njit(cache=True)
def BoundVolume(point):
    V = (np.max(point[:, 0]) - np.min(point[:, 0])) * \
        (np.max(point[:, 1]) - np.min(point[:, 1])) * \
        (np.max(point[:, 2]) - np.min(point[:, 2]))
    return V

@njit(cache=True)
def Distance(point, x):
    # approximate the distance from a point to its nearest point on the superquadric surface
    # extract transformation from superquadric parameters
    R = Euler2RotM(x[5: 8])
    t = x[8: 11]

    # transform to the canonical frame
    point_c = point @ R - t @ R

    # calculating radial distance
    # r_norm = np.linalg.norm(point_c, axis=1)
    r_norm = np.sqrt(np.sum(point_c ** 2, 1))
    
    dist = r_norm * np.abs((
        (((point_c[:, 0] / x[2]) ** 2) ** (1 / x[1]) +
         ((point_c[:, 1] / x[3]) ** 2) ** (1 / x[1])) ** (x[1] / x[0]) +
        ((point_c[:, 2] / x[4]) ** 2) ** (1 / x[0])) ** (-x[0] / 2) - 1
    )
    return dist

@njit(cache=True)
def CostFunc(x, point, p, sigma2):
    if sigma2 > 1e-10:
        value = p ** 0.5 * Distance(point, x)
    else:
        value = np.abs((p * Distance(point, x) ** 2 + 2 *
                       sigma2 * np.log(SurfaceArea(x)))) ** 0.5
    return value

@njit(cache=True)
def OutlierProb(dist, sigma2, w, p0):
    c = (2 * np.pi * sigma2) ** (- 3 / 2)
    const = (w * p0) / (c * (1 - w))
    p = np.exp(-1 / (2 * sigma2) * dist ** 2)
    p = p / (const + p)
    return p

@njit(cache=True)
def SurfaceArea(x):
    a00 = 8 * (x[2] * x[3] + x[3] * x[4] + x[2] * x[4])
    a02 = 8 * (x[2] ** 2 + x[3] ** 2) ** 0.5 * x[4] + 4 * x[2] * x[3]
    a20 = 4 * (x[2] * (x[3] ** 2 + x[4] ** 2) ** 0.5 +
               x[3] * (x[2] ** 2 + x[4] ** 2) ** 0.5)
    a = (x[2] ** 2 + x[3] ** 2) ** 0.5
    b = (x[3] ** 2 + x[4] ** 2) ** 0.5
    c = (x[2] ** 2 + x[4] ** 2) ** 0.5
    s = (a + b + c) / 2
    a22 = 8 * (s * (s - a) * (s - b) * (s - c)) ** (1/2)
    area = np.array([[1 - x[0] / 2, x[0] / 2]]) @ np.array([[a00, a02],
                                                            [a20, a22]]) @ np.array([[1 - x[1] / 2], [x[1] / 2]])
    return area[0, 0]

@njit(cache=True)
def Euler2RotM(euler):
    # from euler angles to rotation matrix (ZYX_intrinsic)

    RotZ = np.array(
        [[np.cos(euler[0]), -np.sin(euler[0]), 0.0],
         [np.sin(euler[0]), np.cos(euler[0]), 0.0],
         [0.0, 0.0, 1.0]]
    )

    RotY = np.array(
        [[np.cos(euler[1]), 0.0, np.sin(euler[1])],
         [0.0, 1.0, 0.0],
         [-np.sin(euler[1]), 0.0, np.cos(euler[1])]]
    )

    RotX = np.array(
        [[1.0, 0.0, 0.0],
         [0.0, np.cos(euler[2]), -np.sin(euler[2])],
         [0.0, np.sin(euler[2]), np.cos(euler[2])]]
    )

    return RotZ @ RotY @ RotX

@njit(cache=True)
def RotM2Euler(R):

    s = np.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])
    singular = s < 1e-6

    if not singular:
        x = np.arctan2(R[2, 1], R[2, 2])
        y = np.arctan2(-R[2, 0], s)
        z = np.arctan2(R[1, 0], R[0, 0])
    else:
        x = np.arctan2(-R[1, 2], R[1, 1])
        y = np.arctan2(-R[2, 0], s)
        z = 0

    return np.array([z, y, x])
