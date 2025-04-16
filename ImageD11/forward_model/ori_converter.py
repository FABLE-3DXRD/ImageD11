# it is adapted from hyperspherical-coverings/quat_utils.py and github transforms3d/euler.py and quaternions.py
# cross verified with GrainRecon_v2 matlab functions
# Haixing Fang, haixing.fang@esrf.fr
# Oct 23, 2023
# Updates on Jan 19th, 2025: implement numba.njit

import numpy as np
import math
from numba import njit

# # axis sequences for Euler angles
# _NEXT_AXIS = [1, 2, 0, 1]

# # map axes strings to/from tuples of inner axis, parity, repetition, frame
# _AXES2TUPLE = {
#     'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
#     'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
#     'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
#     'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
#     'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
#     'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
#     'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
#     'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}

# _TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())

# # For testing whether a number is close to zero
# _EPS4 = np.finfo(float).eps * 4.0

# _FLOAT_EPS = np.finfo(np.float64).eps


def flip_up(q):
    if q[0] < 0:
        q = -q
    return q


def multiply(r, a):

    b0 = r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]
    b1 = r[0] * a[1] + r[1] * a[0] + r[2] * a[3] - r[3] * a[2]
    b2 = r[0] * a[2] - r[1] * a[3] + r[2] * a[0] + r[3] * a[1]
    b3 = r[0] * a[3] + r[1] * a[2] - r[2] * a[1] + r[3] * a[0]
    return flip_up(np.array([b0, b1, b2, b3]))


def random_quaternions(n):

    x0 = np.random.uniform(0, 1, n)
    x1 = np.random.uniform(0, 2 * np.pi, n)
    x2 = np.random.uniform(0, 2 * np.pi, n)

    r1, r2 = np.sqrt(1 - x0), np.sqrt(x0)
    s1, c1 = np.sin(x1), np.cos(x1)
    s2, c2 = np.sin(x2), np.cos(x2)

    return np.array([s1 * r1, c1 * r1, s2 * r2, c2 * r2]).T


def multiply_first_component(r, a):
    return r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]


def rotate_into_fundamental_zone(q, generators):

    index = np.argmax([abs(multiply_first_component(q, g)) for g in generators])
    return multiply(q, generators[index])


def map_points_out(basis_points, basis_weights, superset, subset, map_indices):

    assert( len(map_indices) * len(subset) == len(superset) )

    superset = np.array(superset)
    subset = np.array(subset)

    mapped_points = []
    mapped_weights = []
    for g in superset[map_indices]:
        for b, w in zip(basis_points, basis_weights):
            r = multiply(b, g)
            r = rotate_into_fundamental_zone(r, subset)
            mapped_points += [r]
            mapped_weights += [w]

    return np.array(mapped_points), np.array(mapped_weights)


@njit
def quat2u(q):

    a, b, c, d = q

    u0 = a*a + b*b - c*c - d*d
    u1 = 2*b*c - 2*a*d
    u2 = 2*b*d + 2*a*c

    u3 = 2*b*c + 2*a*d
    u4 = a*a - b*b + c*c - d*d
    u5 = 2*c*d - 2*a*b

    u6 = 2*b*d - 2*a*c
    u7 = 2*c*d + 2*a*b
    u8 = a*a - b*b - c*c + d*d

    return np.array([[u0, u1, u2], [u3, u4, u5], [u6, u7, u8]])


@njit
def rod2quat(r):
	norm_sq = r[0]**2 + r[1]**2 + r[2]**2
	s = 1 / np.sqrt(1 + norm_sq)
	q = np.empty(4, dtype = np.float64)
	q[0] = s
	q[1] = s*r[0]
	q[2] = s*r[1]
	q[3] = s*r[2]
	return q


@njit
def quat2rod(q):
    
    q = np.asarray(q)
    quatmod = np.sqrt( np.sum(q*q) )
    
    qinn = q/quatmod # normalized quaternion
    thalf = np.arccos(qinn[0])
    R = np.zeros(3,)
    if thalf != 0:
        R = qinn[1:4] / np.cos(thalf)
    
    return R


@njit
def u2quat(u):
    
    r11, r12, r13 = u[0]
    r21, r22, r23 = u[1]
    r31, r32, r33 = u[2]

    q0 = (1.0 + r11 + r22 + r33) / 4.0
    q1 = (1.0 + r11 - r22 - r33) / 4.0
    q2 = (1.0 - r11 + r22 - r33) / 4.0
    q3 = (1.0 - r11 - r22 + r33) / 4.0
    q = np.array([q0, q1, q2, q3])
    q = np.sqrt(np.maximum(0, q))

    i = np.argmax(q)
    if i == 0:
        q[1] *= np.sign(r32 - r23)
        q[2] *= np.sign(r13 - r31)
        q[3] *= np.sign(r21 - r12)

    elif i == 1:
        q[0] *= np.sign(r32 - r23)
        q[2] *= np.sign(r21 + r12)
        q[3] *= np.sign(r13 + r31)

    elif i == 2:
        q[0] *= np.sign(r13 - r31)
        q[1] *= np.sign(r21 + r12)
        q[3] *= np.sign(r32 + r23)

    elif i == 3:
        q[0] *= np.sign(r21 - r12)
        q[1] *= np.sign(r31 + r13)
        q[2] *= np.sign(r32 + r23)

    return q / np.linalg.norm(q)


@njit
def quat2euler(q):
    q = np.asarray(q)
    u = quat2u(q)
    return u2euler(u,'rzxz')
#     the following tested to be not correct
#     euler angle conventions are taken from EMSoft
#     qq = q**2
#     q03 = qq[0] + qq[3]
#     q12 = qq[1] + qq[2]
#     chi = np.sqrt(q03 * q12)
#     if chi == 0:
#         if q12 == 0:
#             phi = 0
#             phi2 = 0
#             phi1 = np.arctan2(-2 * q[0] * q[3], qq[0] - qq[3])
#         else:
#             phi = np.pi
#             phi2 = 0
#             phi1 = np.arctan2(2 * q[1] * q[2], qq[1] - qq[2])
#     else:
#         phi = np.arctan2(2 * chi, q03 - q12)
#         phi1 = np.arctan2((-q[0] * q[2] + q[1] * q[3]) / chi, (-q[0] * q[1] - q[2] * q[3]) / chi)
#         phi2 = np.arctan2((+q[0] * q[2] + q[1] * q[3]) / chi, (-q[0] * q[1] + q[2] * q[3]) / chi)

#     res = np.array([phi1, phi, phi2]) % (2 * np.pi)
#     res[1] %= np.pi
#     return res


def euler2quat(e):
    #euler angle conventions are taken from EMSoft
    ee = np.array(e) / 2
    cphi = np.cos(ee[1])
    sphi = np.sin(ee[1])
    cm = np.cos(ee[0] - ee[2])
    sm = np.sin(ee[0] - ee[2])
    cp = np.cos(ee[0] + ee[2])
    sp = np.sin(ee[0] + ee[2])

    res = np.array([cphi * cp, -sphi * cm, -sphi * sm, -cphi * sp])    

    return flip_up(res)


def quat2axangle(quat, identity_thresh=None):
    ''' Convert quaternion to rotation of angle around axis

    Parameters
    ----------
    quat : 4 element sequence
       w, x, y, z forming quaternion.
    identity_thresh : None or scalar, optional
       Threshold below which the norm of the vector part of the quaternion (x,
       y, z) is deemed to be 0, leading to the identity rotation.  None (the
       default) leads to a threshold estimated based on the precision of the
       input.

    Returns
    -------
    theta : scalar
       angle of rotation.
    vector : array shape (3,)
       axis around which rotation occurs.

    Examples
    --------
    >>> vec, theta = quat2axangle([0, 1, 0, 0])
    >>> vec
    array([1., 0., 0.])
    >>> np.allclose(theta, np.pi)
    True

    If this is an identity rotation, we return a zero angle and an arbitrary
    vector:

    >>> quat2axangle([1, 0, 0, 0])
    (array([1., 0., 0.]), 0.0)

    If any of the quaternion values are not finite, we return a NaN in the
    angle, and an arbitrary vector:

    >>> quat2axangle([1, np.inf, 0, 0])
    (array([1., 0., 0.]), nan)

    Notes
    -----
    A quaternion for which x, y, z are all equal to 0, is an identity rotation.
    In this case we return a 0 angle and an arbitrary vector, here [1, 0, 0].

    The algorithm allows for quaternions that have not been normalized.
    '''
    _FLOAT_EPS = np.finfo(np.float64).eps
    
    quat = np.asarray(quat)
    Nq = np.sum(quat ** 2)
    if not np.isfinite(Nq):
        return np.array([1.0, 0, 0]), float('nan')
    if identity_thresh is None:
        try:
            identity_thresh = np.finfo(Nq.type).eps * 3
        except (AttributeError, ValueError): # Not a numpy type or not float
            identity_thresh = _FLOAT_EPS * 3
    if Nq < _FLOAT_EPS ** 2:  # Results unreliable after normalization
        return np.array([1.0, 0, 0]), 0.0
    if Nq != 1:  # Normalize if not normalized
        s = math.sqrt(Nq)
        quat = quat / s
    xyz = quat[1:]
    len2 = np.sum(xyz ** 2)
    if len2 < identity_thresh ** 2:
        # if vec is nearly 0,0,0, this is an identity rotation
        return np.array([1.0, 0, 0]), 0.0
    # Make sure w is not slightly above 1 or below -1
    theta = 2 * math.acos(max(min(quat[0], 1), -1))
    return  xyz / math.sqrt(len2), theta


@njit
def axangle2u(axis, angle, is_normalized=False):
    ''' Rotation matrix for rotation angle `angle` around `axis`

    Parameters
    ----------
    axis : 3 element sequence
       vector specifying axis for rotation.
    angle : scalar
       angle of rotation in radians.
    is_normalized : bool, optional
       True if `axis` is already normalized (has norm of 1).  Default False.

    Returns
    -------
    mat : array shape (3,3)
       rotation matrix for specified rotation

    Notes
    -----
    From: http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    '''
    x, y, z = axis[0], axis[1], axis[2]
    if not is_normalized:
        n = math.sqrt(x*x + y*y + z*z)
        x = x/n
        y = y/n
        z = z/n
    c = math.cos(angle); s = math.sin(angle); C = 1-c
    xs = x*s;   ys = y*s;   zs = z*s
    xC = x*C;   yC = y*C;   zC = z*C
    xyC = x*yC; yzC = y*zC; zxC = z*xC
    return np.array([
            [ x*xC+c,   xyC-zs,   zxC+ys ],
            [ xyC+zs,   y*yC+c,   yzC-xs ],
            [ zxC-ys,   yzC+xs,   z*zC+c ]])


@njit
def euler2u(e):
    # U matrix from Euler angles phi1, PHI, phi2 in radians.
    # INPUT: phi, PHI, and phi2 in radians
    # OUTPUT [U11 U12 U13; U21 U22 U23; U31 U32 U33]
    phi1 = e[0]
    PHI  = e[1]
    phi2 = e[2]
    U11 =  np.cos(phi1)*np.cos(phi2)-np.sin(phi1)*np.sin(phi2)*np.cos(PHI);
    U21 =  np.sin(phi1)*np.cos(phi2)+np.cos(phi1)*np.sin(phi2)*np.cos(PHI);
    U31 =  np.sin(phi2)*np.sin(PHI);
    U12 =  -np.cos(phi1)*np.sin(phi2)-np.sin(phi1)*np.cos(phi2)*np.cos(PHI);
    U22 =  -np.sin(phi1)*np.sin(phi2)+np.cos(phi1)*np.cos(phi2)*np.cos(PHI);
    U32 =  np.cos(phi2)*np.sin(PHI);
    U13 =  np.sin(phi1)*np.sin(PHI);   
    U23 =  -np.cos(phi1)*np.sin(PHI);
    U33 =  np.cos(PHI);
    U = np.array([[U11, U12, U13],
         [U21, U22, U23],
         [U31, U32, U33]]);
    return U


@njit
def u2euler(mat, axes='rzxz'):
    # github transforms3d/euler.py
    # by default using Bunge convention
    """Return Euler angles from rotation matrix for specified axis sequence in [0, 360] degrees.

    Note that many Euler angle triplets can describe one matrix.

    Parameters
    ----------
    mat : array-like shape (3, 3) or (4, 4)
        Rotation matrix or affine.
    axes : str, optional
        Axis specification; one of 24 axis sequences as string or encoded
        tuple - e.g. ``sxyz`` (the default).

    Returns
    -------
    ai : float
        First rotation angle (according to `axes`).
    aj : float
        Second rotation angle (according to `axes`).
    ak : float
        Third rotation angle (according to `axes`).

    Examples
    --------
    >>> R0 = euler2mat(1, 2, 3, 'syxz')
    >>> al, be, ga = mat2euler(R0, 'syxz')
    >>> R1 = euler2mat(al, be, ga, 'syxz')
    >>> np.allclose(R0, R1)
    True
    """
    # axis sequences for Euler angles
    _NEXT_AXIS = [1, 2, 0, 1]
    
    # For testing whether a number is close to zero
    # _EPS4 = np.finfo(float).eps * 4.0   # = 8.881784197001252e-16
    _EPS4 = 8.881784197001252e-16

#     # map axes strings to/from tuples of inner axis, parity, repetition, frame
#     _AXES2TUPLE = {
#         'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
#         'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
#         'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
#         'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
#         'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
#         'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
#         'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
#         'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}

#     _TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())
    _AXES2TUPLE = np.array([(0, 0, 0, 0), (0, 0, 1, 0), (0, 1, 0, 0),
                           (0, 1, 1, 0), (1, 0, 0, 0), (1, 0, 1, 0),
                           (1, 1, 0, 0), (1, 1, 1, 0), (2, 0, 0, 0),
                           (2, 0, 1, 0), (2, 1, 0, 0), (2, 1, 1, 0),
                           (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 0, 1),
                           (0, 1, 1, 1), (1, 0, 0, 1), (1, 0, 1, 1),
                           (1, 1, 0, 1), (1, 1, 1, 1), (2, 0, 0, 1),
                           (2, 0, 1, 1), (2, 1, 0, 1), (2, 1, 1, 1)], dtype='int32')
    axes_names = ['sxyz', 'sxyx', 'sxzy',
                 'sxzx', 'syzx', 'syzy',
                 'syxz', 'syxy', 'szxy',
                 'szxz', 'szyx', 'szyz',
                 'rzyx', 'rxyx', 'ryzx',
                 'rxzx', 'rxzy', 'ryzy',
                 'rzxy', 'ryxy', 'ryxz',
                 'rzxz', 'rxyz', 'rzyz']
    # Find the index of the given axis
    axes_indice = -1
    for i in range(len(axes_names)):
        if axes == axes_names[i]:
            axes_indice = i
            break
    if axes_indice == -1:
        raise ValueError("Invalid input of axes")

    # Parse axes
    firstaxis, parity, repetition, frame = _AXES2TUPLE[axes_indice]
               
    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    # M = np.array(mat, dtype=np.float64, copy=False)[:3, :3]
    M = np.asarray(mat, dtype=np.float64)[:3, :3]
    if repetition:
        sy = math.sqrt(M[i, j]*M[i, j] + M[i, k]*M[i, k])
        if sy > _EPS4:
            ax = math.atan2( M[i, j],  M[i, k])
            ay = math.atan2( sy,       M[i, i])
            az = math.atan2( M[j, i], -M[k, i])
        else:
            ax = math.atan2(-M[j, k],  M[j, j])
            ay = math.atan2( sy,       M[i, i])
            az = 0.0
    else:
        cy = math.sqrt(M[i, i]*M[i, i] + M[j, i]*M[j, i])
        if cy > _EPS4:
            ax = math.atan2( M[k, j],  M[k, k])
            ay = math.atan2(-M[k, i],  cy)
            az = math.atan2( M[j, i],  M[i, i])
        else:
            ax = math.atan2(-M[j, k],  M[j, j])
            ay = math.atan2(-M[k, i],  cy)
            az = 0.0

    if parity:
        ax, ay, az = -ax, -ay, -az
    if frame:
        ax, az = az, ax
    
    # move angles to the range of [0, 2*pi]
    if ax < 0:
        ax = ax + 2*np.pi
    if ay < 0:
        ay = ay + 2*np.pi
    if az < 0:
        az = az + 2*np.pi
    
    return np.rad2deg(ax), np.rad2deg(ay), np.rad2deg(az)
