exec(open('/data/id11/nanoscope/install_ImageD11_from_git.py').read())
# use my own python for this during development
PYTHONPATH = setup_ImageD11_from_git(os.path.join( os.environ['HOME'],'Code'), 'ImageD11' )

import os
import time
import numba
import numpy as np

nthreads = len(os.sched_getaffinity(os.getpid()))
numba.set_num_threads(nthreads)

# from ImageD11 import cImageD11
# # we do not want the processes to conflict
# os.environ['OMP_NUM_THREADS'] = '1'
# os.environ['OPENBLAS_NUM_THREADS'] = '1'
# os.environ['MKL_NUM_THREADS'] = '1'
# cImageD11.check_multiprocessing(patch=True) # forkserver
# cImageD11.cimaged11_omp_set_num_threads(1)

import ImageD11.sinograms.dataset
from ImageD11.sinograms import geometry
from ImageD11 import sym_u
from ImageD11 import transform as transformid11
from ImageD11.sinograms.point_by_point import PBPMap, PBPRefine
from ImageD11.sinograms.tensor_map import unitcell_to_b
import xfab


def get_point_data(params, xi0, yi0, peak_dict):
    # params: imaged11 pars
    # xi0, yi0: sx_grid[i, j], sy_grid[i, j]
    # peak_dict: the peaks
    
    # returns:
    # xpos: x coordinate of voxel in lab space
    # m: peak mask for this voxel
    # ydist: distance of all peaks from this voxel
    
    # does:
    # converts points in sample space to points in lab space
    # masks peaks to this voxel
    # returns lab space position, peak
    xpos = xi0*peak_dict['cosomega'] - yi0*peak_dict['sinomega']
    y = xi0*peak_dict['sinomega'] + yi0*peak_dict['cosomega']
    ydist = np.abs( y + peak_dict['dty'] )
    m = ydist <= params['TRANSLATION_STAGE_DY_STEP']
    # print('Voxel mask got', m.sum(), 'peaks')
    # print(np.arange(len(peak_dict['dty']))[m])
    return xpos, m, ydist

def get_gve(params, local_mask, xpos, peak_dict):
    # params: image11 pars
    # local_mask: peak mask for this voxel
    # xpos: x coordinate of voxel in lab space
    # peak_dict: the peaks
    
    # returns: corrected g-vectors
    D = params['distance']
    params['distance'] = params['distance'] - xpos[local_mask]
    tth, eta = ImageD11.transform.compute_tth_eta(
                        (peak_dict['sc'][local_mask], peak_dict['fc'][local_mask]),
                        **params)
    gve = ImageD11.transform.compute_g_vectors(tth,
                        eta,
                        peak_dict['omega'][local_mask],
                        params['wavelength'],
                        wedge=params['wedge'],
                        chi=params['chi'])
    params['distance'] = D
    return gve


def refine_point_axel(xi0, yi0, ubis, params, peak_dict, lattice):
    # xi0, yi0: sx_grid[i, j], sy_grid[i, j]
    # ubis: ubis at [i, j]
    # params: imaged11 pars
    # peak_dict: the peaks
    # lattice: ImageD11.sym_u.cubic()
    
    # get the masked peaks for this voxel
    _, local_mask, ydist = get_point_data(params, xi0, yi0, peak_dict)
    
    print('First local mask results...')
    print(peak_dict['sinomega'][local_mask].sum())
    
    # compute g-vectors at this voxel
    gve = get_gve(params, local_mask, peak_dict['xpos_refined'], peak_dict)
    
    # convert ubis to grains
    grains = [ImageD11.grain.grain(ubi) for ubi in ubis]

    tol = 0.1
    merge_tol = 0.05
    
    # iterate through grains
    for i,g in enumerate(grains):
        
        labels = np.zeros( gve.shape[1],'i')-1
        drlv2 = np.ones( gve.shape[1], 'd') 
        
        # score and assign the grain to the peaks at this voxel
        j = ImageD11.cImageD11.score_and_assign( g.ubi, gve.T, tol, drlv2, labels, i)
        li = np.round(labels).astype(int).copy()
        
        # which peaks at this voxel does this grain like?
        g.mask = g.m = (li == i)
        
        print('Assignment mask results...')
        print(peak_dict['sinomega'][local_mask][g.mask].sum())
        
        # print('Assigned', g.mask.sum(), 'peaks first time')
        # how many peaks does this grain have?
        g.npks = np.sum(g.mask)
        
        print('Assigned', g.npks, 'peaks first time out of', local_mask.sum())
        
        hkl_double = np.dot( g.ubi, gve[:, g.mask] ) 
        g.hkl = np.round ( hkl_double ).astype( int )
        g.etasign = np.sign( peak_dict['eta'][ local_mask ][ g.m ] )
        
        # 3D merge peaks in omega
        # only the peaks that this grain indexes
        merged = {'sc':[], 'fc':[], 'dty':[], 'omega':[], 'sum_intensity':[], 'xpos_refined':[]}

        peaktags = np.vstack( (g.hkl, g.etasign, peak_dict['iy'][ local_mask ][ g.m ]) )
        unitags, labels  = np.unique(peaktags, axis=1, return_inverse=True )
        print(labels[:10])
        print('Got', unitags.shape[1], 'unique peaks during merging')
        print(unitags[:, 0].sum())
        wI = peak_dict['sum_intensity'][ local_mask ][g.m]
        sI = np.bincount( labels, weights=wI )
        print('wI sum:', wI.sum())
        print('sI sum:', sI.sum())
        merged['sum_intensity']=sI
        sc = peak_dict['sc'][ local_mask ][g.m]
        fc = peak_dict['fc'][ local_mask ][g.m]
        om = peak_dict['omega'][ local_mask ][g.m]
        dty = peak_dict['dty'][ local_mask ][g.m]
        xpos_refined = peak_dict['xpos_refined'][ local_mask ][g.m]
        merged['sc'].extend( list( np.bincount( labels, weights=sc*wI  ) / sI ) )
        merged['fc'].extend( list( np.bincount( labels, weights=fc*wI  ) / sI ) )
        merged['omega'].extend( list( np.bincount( labels, weights=om*wI  ) / sI ) )
        merged['dty'].extend( list( np.bincount( labels, weights=dty*wI  ) / sI ) )
        merged['xpos_refined'].extend( list( np.bincount( labels, weights=xpos_refined*wI  ) / sI ) )

        for k in merged.keys():
            merged[k] = np.array(merged[k])
        
        print('merged_sc sum:',  merged['sc'].sum())
        
        merged['sinomega'] = np.sin( np.radians(merged['omega']) )
        merged['cosomega'] = np.cos( np.radians(merged['omega']) )
        
        print('merged_sc has shape', merged['sc'].shape)
        
        print('Merged', len(sc), 'peaks to', len(merged['sc']), 'peaks')
        
        print('Merged peak data for masking...')
        print('merged so sum', merged['sinomega'].sum())
        print('merged co sum', merged['cosomega'].sum())
        print('merged dty sum', merged['dty'].sum())
        
        # recompute the peaks this grain likes with merged data
        _, local_mask_of_grain, ydist = get_point_data(params, xi0, yi0, merged)
        print('Got', local_mask_of_grain.sum(), 'masked merged peaks')
        print('ydist sum before merged masking', ydist.sum())
        # recompute gvector of grain with merged data
        gve_grain = get_gve(params, local_mask_of_grain, merged['xpos_refined'], merged)
        
        print('gve_grain sum before merged assignment', gve_grain.sum())
        
        # reassign this grain with merged peaks data
        labels = np.zeros( gve_grain.shape[1],'i')-1
        drlv2 = np.ones( gve_grain.shape[1], 'd') 
        j = ImageD11.cImageD11.score_and_assign( g.ubi, gve_grain.T, merge_tol, drlv2, labels, i)
        li = np.round(labels).astype(int).copy()
        g.mask = g.m = (li == i)
        g.npks = np.sum(g.m)
        print('Assigned', g.npks, 'merged peaks out of', local_mask_of_grain.sum())
        g.gve = gve_grain[:,g.mask]
        hkl_double = np.dot( g.ubi, gve_grain[:, g.mask] ) 
        g.hkl = np.round ( hkl_double ).astype( int )
        g.ydist = ydist[local_mask_of_grain][g.mask]
        
        # print(g.npks, 'merged peaks assigned after merging')
        
        # if grain has more than 6 peaks assigned
        if np.sum(g.mask) > 6:
            
            
            print('Going into refinement...')
            print('ydist sum', g.ydist.sum())
            print('gve sum', g.gve.sum())
            print('hkl sum', g.hkl.sum())
            
            # a.T @ gve = h =>  gve.T @ a = h.T => a = np.linalg.pinv(gve.T) @ h.T, same for b and c 
            w = (1. / (g.ydist + 1) ).reshape(g.gve.shape[1], 1)
            ubifitT, residuals, rank, sing_vals = np.linalg.lstsq( w * g.gve.T, w * g.hkl.T, rcond=None )
            ubifit = ubifitT.T
            
            # TODO: fix bad cells....  and rank==3 and np.linalg.cond(ubifit)<1e16
            if ubifit is not None and rank==3 and np.linalg.cond(ubifit)<1e14 and np.linalg.det(ubifit) > 0 and np.linalg.matrix_rank(ubifit)==3:                    
                    
                _G = np.dot(ubifit, ubifit.T)
                _a, _b, _c = np.sqrt(np.diag(_G))

                vals = [(_G[1, 2] / _b / _c), (_G[0, 2] / _a / _c), (_G[0, 1] / _a / _b)]
                for v in vals:

                    if not isinstance(v, float):
                        print(ubifit, np.linalg.matrix_rank(ubifit), sing_vals, rank, np.linalg.det(ubifit), v, w)
                        raise

                    if np.imag(v)!=0:
                        print(ubifit, np.linalg.matrix_rank(ubifit), sing_vals, rank, np.linalg.det(ubifit), v, w)
                        raise

                    if v <=-1 or v >=1:
                        print(ubifit, np.linalg.matrix_rank(ubifit), sing_vals, rank, np.linalg.det(ubifit), v, w)
                        raise

                al = np.degrees(np.arccos(_G[1, 2] / _b / _c))
                be = np.degrees(np.arccos(_G[0, 2] / _a / _c))
                ga = np.degrees(np.arccos(_G[0, 1] / _a / _b))

                assert _a>0, str(ubifit)
                assert _b>0, str(ubifit)
                assert _c>0, str(ubifit)
                assert al>0 and al<180, str(ubifit)
                assert be>0 and be<180, str(ubifit)
                assert ga>0 and ga<180, str(ubifit)


                g.set_ubi( ubifit )

            labels = np.zeros( gve_grain.shape[1],'i')-1
            drlv2 = np.ones( gve_grain.shape[1], 'd') 

            j = ImageD11.cImageD11.score_and_assign( g.ubi, gve_grain.T, merge_tol, drlv2, labels, i)
            li = np.round(labels).astype(int).copy()
            g.mask = g.m = (li == i)
            g.npks = np.sum(g.m)

            g.eps_tensor  = None
            if ubifit is not None and rank==3 and np.linalg.cond(ubifit)<1e14 and np.linalg.det(ubifit) > 0 and np.linalg.matrix_rank(ubifit)==3:                    

                ####################################################################################################
                # And now fit a strain tensor decoupling strain and orientation for precision.
                
                print('Strain mask fit sum', g.mask.sum())
                
                g.gve = gve_grain[:,g.mask]
                hkl_double = np.dot( g.ubi, gve_grain[:, g.mask] ) 
                g.hkl = np.round ( hkl_double ).astype( int )
                g.ydist = ydist[local_mask_of_grain][g.mask]
                
                print('gve_grain_strainfit sum:',  g.gve.sum())
                print('ydist strainfit sum:', g.ydist.sum())

                dzero_cell = [params['cell_'+key] for key in ('_a','_b','_c','alpha','beta','gamma')]
                B0 = xfab.tools.form_b_mat(dzero_cell) / (np.pi*2)
                print('B0')
                print(B0)

                g.gve0 = g.u @ B0 @ g.hkl

                gTg0    = np.sum(g.gve*g.gve0, axis=0)
                gTg   = np.sum(g.gve*g.gve, axis=0)
                g.directional_strain  = (gTg0/gTg) - 1
                
                print('Directional strain')
                print(g.directional_strain)

                kappa = g.gve / np.linalg.norm(g.gve, axis=0)
                kx, ky, kz = kappa
                M = np.array( [ kx*kx, ky*ky, kz*kz, 2*kx*ky, 2*kx*kz, 2*ky*kz ] ).T

                w = (1. / (g.ydist + 1) ).reshape(g.gve.shape[1], 1)
                # The noise in the directional strain now propagates according to the linear transform
                gnoise_std = 1e-4
                a  = np.sum(g.gve0*(gnoise_std**2)*g.gve0, axis=0)
                strain_noise_std = np.sqrt( np.divide(a, gTg**2, out=np.ones_like(gTg), where=gTg!=0) )
                w = w * (1. / strain_noise_std.reshape(w.shape) )

                w[ g.directional_strain > np.mean(g.directional_strain) + np.std(g.directional_strain)*3.5 ] = 0 # outliers
                w[ g.directional_strain < np.mean(g.directional_strain) - np.std(g.directional_strain)*3.5 ] = 0 # outliers

                try:
                    w = w / np.max(w)
                    eps_vec = np.linalg.lstsq( w * M, w.flatten() * g.directional_strain, rcond=None )[0].T
                    sxx, syy, szz, sxy, sxz, syz = eps_vec
                    g.eps_tensor = np.array([[sxx, sxy, sxz],[sxy, syy, syz],[sxz, syz, szz]])
                    print(g.eps_tensor)
                except:
                    pass

                ###############################################################################################################

        else:
            g.eps_tensor  = None
    
    # sort grains array by npeaks
    grains = np.array(grains)[np.argsort([g.npks for g in grains])]
    # get symmetry-unique ubis of grains
    ubis = [ImageD11.sym_u.find_uniq_u(g.ubi, lattice) for g in grains]
    npks = [g.npks for g in grains]
    eps  = [g.eps_tensor for g in grains]
    
    return ubis, eps, npks



def call_axel(refine, points):
    print('Axel prep')
    # points are positions in the array
    all_pbpmap_ubis = refine.pbpmap.ubi
    all_npks = refine.pbpmap.ntotal
    
    # get params, peak_dict, lattice
    params = refine.icolf.parameters.get_parameters()
    params['TRANSLATION_STAGE_DY_STEP'] = refine.dset.ystep
    
    peak_dict = {
        'sc': refine.icolf.sc,
        'fc': refine.icolf.fc,
        'omega': refine.icolf.omega,
        'dty': refine.icolf.dty,
        'sinomega': refine.icolf.sinomega,
        'cosomega': refine.icolf.cosomega,
        'gx': refine.icolf.gx,
        'gy': refine.icolf.gy,
        'gz': refine.icolf.gz,
        'iy': refine.icolf.dtyi,
        'sum_intensity': refine.icolf.sum_intensity,
        'eta': refine.icolf.eta,
        'xpos_refined': refine.icolf.xpos_refined
    }
    lattice = sym_u.cubic()
    
    refine_inputs = []
    refine_results = []
    
    print('Iterating over points')
    for (si, sj) in points:
        # si, sj is in step space (could have negative values)
        # need to convert to reconstruction space to get sx_grid and sy_grid values
        ri, rj = geometry.step_to_recon(si, sj, refine.sx_grid.shape)
        # get pbpmap rows at this point
        pixel_mask = refine.pbpmap.get_pixel_mask(si, sj)
        # get ubis at this point
        ubis_in = list(np.rollaxis(all_pbpmap_ubis[:, :, pixel_mask], 2))
        # get npks at this point
        npks_in = list(all_npks[pixel_mask])
        
        # get sample positions at this point
        xi0 = refine.sx_grid[ri, rj]
        yi0 = refine.sy_grid[ri, rj]
        
        print('At ri rj', ri, rj)
        ubis_out, eps, npks_out = refine_point_axel(xi0, yi0, ubis_in, params, peak_dict, lattice)
        refine_inputs.append((ubis_in, npks_in))
        refine_results.append((ubis_out, eps, npks_out))
        
    return refine_inputs, refine_results

# stuff we need to compute g-vectors
# from ImageD11.transform
@numba.njit
def detector_rotation_matrix(tilt_x, tilt_y, tilt_z):
    """
    Return the tilt matrix to apply to peaks
    tilts are in radians
    typically applied to peaks rotating around beam center
    """
    r1 = np.array([[np.cos(tilt_z), -np.sin(tilt_z), 0.0],  # note this is r.h.
                   [np.sin(tilt_z), np.cos(tilt_z), 0.0],
                   [0.0,    0.0, 1.0]], np.float64)
    r2 = np.array([[np.cos(tilt_y), 0.0, np.sin(tilt_y)],
                   [0.0, 1.0,   0.0],
                   [-np.sin(tilt_y), 0.0, np.cos(tilt_y)]], np.float64)
    r3 = np.array([[1.0,          0.0,       0.0],
                   [0.0,  np.cos(tilt_x), -np.sin(tilt_x)],
                   [0.0,  np.sin(tilt_x), np.cos(tilt_x)]], np.float64)
    r2r1 = np.dot(np.dot(r3, r2), r1)
    return r2r1

@numba.njit
def compute_xyz_lab(sc, fc,
                    y_center=0., y_size=0., tilt_y=0.,
                    z_center=0., z_size=0., tilt_z=0.,
                    tilt_x=0.,
                    distance=0.,
                    o11=1.0, o12=0.0, o21=0.0, o22=-1.0):
    """
    Peaks is a 2 d array of x,y
    yc is the centre in y
    ys is the y pixel size
    ty is the tilt around y
    zc is the centre in z
    zs is the z pixel size
    tz is the tilt around z
    dist is the sample - detector distance
    detector_orientation is a matrix to apply to peaks arg to get
    ImageD11 convention
         (( 0, 1),( 1, 0)) for ( y, x)
         ((-1, 0),( 0, 1)) for (-x, y)
         (( 0,-1),(-1, 0)) for (-y,-x)
      etc...
      
    kwds are not used (but lets you pass in a dict with other things in it)
    """

    # Matrix for the tilt rotations
    r2r1 = detector_rotation_matrix(tilt_x, tilt_y, tilt_z)
    # Peak positions in 3D space
    #  - apply detector orientation
    peaks_on_detector = np.stack((sc, fc))
    peaks_on_detector[0, :] = (peaks_on_detector[0, :] - z_center) * z_size
    peaks_on_detector[1, :] = (peaks_on_detector[1, :] - y_center) * y_size
    #
    detector_orientation = [[o11, o12], [o21, o22]]
    # logging.debug("detector_orientation = "+str(detector_orientation))
    flipped = np.dot(np.array(detector_orientation, np.float64),
                     peaks_on_detector)
    #
    # vec = np.array([np.zeros(flipped.shape[1]),  # place detector at zero,
    #                 # sample at -dist
    #                 flipped[1, :],             # x in search, frelon +z
    #                 flipped[0, :]], np.float64)     # y in search, frelon -y
    
    vec = np.stack((np.zeros(flipped.shape[1]), flipped[1, :], flipped[0, :]))
    
    # Position of diffraction spots in 3d space after detector tilts about
    # the beam centre on the detector
    rotvec = np.dot(r2r1, vec)
    # Now add the distance (along x)
    rotvec[0, :] = rotvec[0, :] + distance
    return rotvec

@numba.njit
def compute_tth_eta_from_xyz(peaks_xyz,
                             t_x=0.0, t_y=0.0, t_z=0.0,
                             wedge=0.0,  # Wedge == theta on 4circ
                             chi=0.0):  # last line is for laziness -
    """
    Peaks is a 3 d array of x,y,z peak co-ordinates
    crystal_translation is the position of the grain giving rise to a diffraction spot
    in x,y,z ImageD11 co-ordinates
         x,y is with respect to the axis of rotation (usually also beam centre).
         z with respect to beam height, z centre
    omega data are needed if crystal translations are used
    
    computed via the arctan recipe.
    
    returns tth/eta in degrees
    """
    s1 = peaks_xyz
    
    

    # CHANGED to HFP convention 4-9-2007
    eta = np.degrees(np.arctan2(-s1[1, :], s1[2, :]))
    s1_perp_x = np.sqrt(s1[1, :] * s1[1, :] + s1[2, :] * s1[2, :])
    tth = np.degrees(np.arctan2(s1_perp_x, s1[0, :]))
    return tth, eta

@numba.njit
def compute_tth_eta(sc, fc,
                    y_center=0., y_size=0., tilt_y=0.,
                    z_center=0., z_size=0., tilt_z=0.,
                    tilt_x=0.,
                    distance=0.,
                    o11=1.0, o12=0.0, o21=0.0, o22=-1.0,
                    t_x=0.0, t_y=0.0, t_z=0.0,
                    wedge=0.0,
                    chi=0.0):
    """
    Finds x,y,z co-ordinates of peaks in the laboratory frame
    Computes tth/eta from these (in degrees)
    
    kwds are not used (left for convenience if you have a parameter dict)
    """
    peaks_xyz = compute_xyz_lab(
        sc, fc,
        y_center=y_center, y_size=y_size, tilt_y=tilt_y,
        z_center=z_center, z_size=z_size, tilt_z=tilt_z,
        tilt_x=tilt_x,
        distance=distance,
        o11=o11, o12=o12, o21=o21, o22=o22)

    tth, eta = compute_tth_eta_from_xyz(
        peaks_xyz,
        t_x=t_x, t_y=t_y, t_z=t_z,
        wedge=wedge,
        chi=chi)

    return tth, eta

@numba.njit
def compute_k_vectors(tth, eta, wvln):
    """
    generate k vectors - scattering vectors in laboratory frame
    """
    tth = np.radians(tth)
    eta = np.radians(eta)
    c = np.cos(tth / 2)  # cos theta
    s = np.sin(tth / 2)  # sin theta
    ds = 2 * s / wvln
    k = np.zeros((3, tth.shape[0]), np.float64)
    # x - along incident beam
    k[0, :] = -ds * s  # this is negative x
    # y - towards door
    k[1, :] = -ds * c * np.sin(eta)  # CHANGED eta to HFP convention 4-9-2007
    # z - towards roof
    k[2, :] = ds * c * np.cos(eta)
    return k

@numba.njit
def compute_g_from_k(k, omega, wedge=0, chi=0):
    """
    Compute g-vectors with cached k-vectors
    """
    om = np.radians(omega)
    # G-vectors - rotate k onto the crystal axes
    g = np.zeros((3, k.shape[1]), np.float64)
    t = np.zeros((3, k.shape[1]), np.float64)
    #
    # g =  R . W . k where:
    # R = ( cos(omega) , sin(omega), 0 )
    #     (-sin(omega) , cos(omega), 0 )
    #     (         0  ,         0 , 1 )
    #
    # W = ( cos(wedge) ,  0  ,  sin(wedge) )
    #     (         0  ,  1  ,          0  )
    #     (-sin(wedge) ,  0  ,  cos(wedge) )
    #
    # C = (         1  ,         0  ,      0     )
    #     (         0  ,  cos(chi)  , sin(chi)   )
    #     (         0  , -sin(chi)  , cos(chi)   )
    #
    if wedge != 0.0:
        c = np.cos(np.radians(wedge))
        s = np.sin(np.radians(wedge))
        t[0, :] = c * k[0, :] + s * k[2, :]
        t[1, :] = k[1, :]
        t[2, :] = -s * k[0, :] + c * k[2, :]
        k = t.copy()
    if chi != 0.0:
        c = np.cos(np.radians(chi))
        s = np.sin(np.radians(chi))
        t[0, :] = k[0, :]
        t[1, :] = c * k[1, :] + s * k[2, :]
        t[2, :] = -s * k[1, :] + c * k[2, :]
        k = t.copy()
    # This is the reverse rotation (left handed, k back to g)
    g[0, :] = np.cos(om) * k[0, :] + np.sin(om) * k[1, :]
    g[1, :] = -np.sin(om) * k[0, :] + np.cos(om) * k[1, :]
    g[2, :] = k[2, :]
    return g

@numba.njit
def compute_g_vectors(tth,
                      eta,
                      omega,
                      wvln,
                      wedge=0.0,
                      chi=0.0):
    """
    Generates spot positions in reciprocal space from
      twotheta, wavelength, omega and eta
    Assumes single axis vertical
    ... unless a wedge angle is specified
    """
    k = compute_k_vectors(tth, eta, wvln)
    return compute_g_from_k(k, omega, wedge, chi)


# from ImageD11/src/closest.c
@numba.njit
def score_and_assign(ubi, gvecs, tol, label):
    tolsq = tol ** 2

    labels = np.zeros(gvecs.shape[1],'i') - 1
    drlv2 = np.ones(gvecs.shape[1],'d')
    
    n = 0
    for gv_idx in np.arange(gvecs.shape[1]):
        gv = gvecs[:,gv_idx]
        hklf = ubi.dot(gv)
        hkle = np.round(hklf)
        err = hklf - hkle
        sumsq = (err ** 2).sum()

        if ((sumsq < tolsq) & (sumsq < drlv2[gv_idx])):
            labels[gv_idx] = label
            drlv2[gv_idx] = sumsq
            n += 1
        elif labels[gv_idx] == label:
            labels[gv_idx] = -1

    return n, labels, drlv2

@numba.njit
def count_unique_peaks(hkl, etasign, dtyi):

    # Combine the input arrays into one 2D array for easier sorting
    combined = np.empty((5, etasign.shape[0]), dtype=dtyi.dtype)
    combined[0:3, :] = hkl
    combined[3, :] = etasign
    combined[4, :] = dtyi
    
    # Create an array of indices to track the original order
    indices = np.arange(combined.shape[1])
    
    # Custom insertion sort to sort lexicographically
    def lexicographical_sort(arr, inds):
        for i in range(1, arr.shape[1]):
            key = arr[:, i].copy()  # Copy the ith column
            key_index = inds[i]  # Get the corresponding index
            j = i - 1
            # Lexicographical comparison: compare all elements in the column
            while j >= 0 and (arr[0, j] > key[0] or
                              (arr[0, j] == key[0] and arr[1, j] > key[1]) or
                              (arr[0, j] == key[0] and arr[1, j] == key[1] and arr[2, j] > key[2]) or
                              (arr[0, j] == key[0] and arr[1, j] == key[1] and arr[2, j] == key[2] and arr[3, j] > key[3]) or
                              (arr[0, j] == key[0] and arr[1, j] == key[1] and arr[2, j] == key[2] and arr[3, j] == key[3] and arr[4, j] > key[4])):
                arr[:, j + 1] = arr[:, j]
                inds[j + 1] = inds[j]  # Shift the indices accordingly
                j -= 1
            arr[:, j + 1] = key
            inds[j + 1] = key_index
    
    # Sort the combined array lexicographically
    lexicographical_sort(combined, indices)
    
    uniques = np.empty((5, etasign.shape[0]), dtype=dtyi.dtype)
    inverses = np.empty((etasign.shape[0]), dtype=np.int64)
    unique_ptr = 0
    
    # iterate through the sorted input array
    for idx_sorted in np.arange(combined.shape[1], dtype=np.int64):
        h, k, l, this_etasign, this_dtyi = combined[:, idx_sorted]
        # go from the index of the sorted array back to the index of the original array
        idx_combined = indices[idx_sorted]
        
        # did we see this row in the unique array yet?
        seen_array = False
        # iterate through the unique array
        for idx_unique in np.arange(unique_ptr, dtype=np.int64):
            # get row in the unique array
            h_uniq, k_uniq, l_uniq, etasign_uniq, dtyi_uniq = uniques[:, idx_unique]
            if (h == h_uniq) & (k == k_uniq) & (l == l_uniq) & (this_etasign == etasign_uniq) & (this_dtyi == dtyi_uniq):
                # saw it already!
                inverses[idx_combined] = idx_unique
                seen_array = True
                break
                
        if not seen_array:
            # got to the end of uniques and didn't see it
            uniques[:, unique_ptr] = [h, k, l, this_etasign, this_dtyi]
            inverses[idx_combined] = unique_ptr
            unique_ptr += 1
    
    return uniques[:,:unique_ptr], inverses


@numba.njit
def merge(hkl, etasign, dtyi, sum_intensity, sc, fc, omega, dty, xpos_refined):
    unitags, labels = count_unique_peaks(hkl, etasign, dtyi)
    print('Got', unitags.shape[1], 'unique peaks during merging')
    print(unitags[:, 0].sum())
    print(labels[:10])
    wI = sum_intensity
    sI = np.bincount( labels, weights=wI )
    print('wI sum:', wI.sum())
    print('sI sum:', sI.sum())
    
    merged_sum_intensity = sI
    merged_sc = np.bincount( labels, weights=sc*wI  )/sI
    merged_fc = np.bincount( labels, weights=fc*wI  )/sI
    merged_omega = np.bincount( labels, weights=omega*wI  )/sI
    merged_dty = np.bincount( labels, weights=dty*wI  )/sI
    merged_xpos_refined = np.bincount( labels, weights=xpos_refined*wI  )/sI
    
    print('merged_sc has shape', merged_sc.shape)
    print('merged_sc sum:', merged_sc.sum())

    return merged_sum_intensity, merged_sc, merged_fc, merged_omega, merged_dty, merged_xpos_refined

@numba.njit
def get_voxel_mask(y0, xi0, yi0, sinomega, cosomega, dty, ystep):
    # geometry.dtycalc_sincos
    dty_calc = y0 - xi0 * sinomega - yi0 * cosomega
    ydist = np.abs(dty_calc - dty)
    m = ydist <= ystep
    # print('Voxel mask got', m.sum(), 'peaks')
    return m, ydist


@numba.njit
def compute_gve(sc, fc, omega, xpos,
                distance, y_center, y_size, tilt_y, z_center, z_size, tilt_z, tilt_x,
                o11, o12, o21, o22,
                t_x, t_y, t_z, wedge, chi, wavelength):
        
        # compute g-vectors at this voxel
        this_distance = distance - xpos
        
        tth, eta = compute_tth_eta(sc, fc,
                                y_center=y_center,
                                y_size=y_size,
                                tilt_y=tilt_y,
                                z_center=z_center,
                                z_size=z_size,
                                tilt_z=tilt_z,
                                tilt_x=tilt_x,
                                distance=this_distance,
                                o11=o11,
                                o12=o12,
                                o21=o21,
                                o22=o22,
                                t_x=t_x,
                                t_y=t_y,
                                t_z=t_z,
                                wedge=wedge,
                                chi=chi
                               )
        
        gve = compute_g_vectors(tth, eta,
                                omega,
                                wavelength,
                                wedge=wedge,
                                chi=chi)
        
        return gve


@numba.njit
def weighted_lstsq_ubi_fit(ydist, gve, hkl):
    # run the weighted fit
    # a.T @ gve = h =>  gve.T @ a = h.T => a = np.linalg.pinv(gve.T) @ h.T, same for b and c 
    w = (1. / (ydist + 1) ).reshape(gve.shape[1], 1)
    a = w * gve.T
    b = w * hkl.T
    m, n = a.shape[-2:]
    rcond = np.finfo(b.dtype).eps * max(n, m)
    ubifitT, residuals, rank, sing_vals = np.linalg.lstsq(a, b, rcond=rcond)
    ubifit = ubifitT.T
    
    return w, ubifit, residuals, rank, sing_vals

@numba.njit
def gve_norm(gve):
    norms = np.zeros(gve.shape[1])
    for i in range(gve.shape[1]):
        gv = gve[:,i]
        norms[i] = np.sqrt(np.sum(np.power(gv, 2)))
    
    return norms


@numba.njit
def divide_where(arr1, arr2, out, wherearr):
    """Do arr1/arr2.
    In locations where wherearr == 0, return out instead"""
    div = np.divide(arr1, arr2)
    return np.where(wherearr != 0, div, out)

@numba.njit
def mine(refine_points, all_pbpmap_ubis, ri_col, rj_col, sx_grid, sy_grid,  # refinement stuff
         sc, fc, eta, sum_intensity, sinomega, cosomega, omega, dty, dtyi, xpos,  # icolf columns
         ystep, y0,
         B0,
         distance, y_center, y_size, tilt_y, z_center, z_size, tilt_z, tilt_x,
         o11, o12, o21, o22,
         t_x, t_y, t_z, wedge, chi, wavelength,
         tol=0.1, merge_tol=0.05
):
    # refine_points: list or Nx2 array of points in reconstruction space (ri, rj) to refine at
    # all_pbpmap_ubis: 3x3xM array of many-valued UBIs from point-by-point map
    # ri_col, rj_col: length-M array of ri and rj coordinates (reconstruction space) for each UBI in all_pbpmap_ubis
    # sx_grid, sy_grid: same as before (sample space sx, sy positions)
    # sinomega, cosomega, dtyi: same as before (icolf columns)
    # ystep, y0: same as before
    
    # iterate through point indices
    # they are in reconstruction space
    for refine_idx in numba.prange(len(refine_points)):
        ri, rj = refine_points[refine_idx]
        print('At ri rj', ri, rj)
        
        # mask all_ubis by the pbpmap points
        pbpmap_idx = (ri_col == ri) & (rj_col == rj)
        
        # get ubis at this point
        ubis_here = all_pbpmap_ubis[:, :, pbpmap_idx]

        # get xi0, xi0 at this point
        xi0 = sx_grid[ri, rj]
        yi0 = sy_grid[ri, rj]
        
        # get a mask to the peaks at this point
        # this is basically geometry.dtyimask_from_sincos
        # but we already have x, y
        
        local_mask, _ = get_voxel_mask(y0, xi0, yi0, sinomega, cosomega, dty, ystep)
        
        # does our local masking agree?
        
        print('First local mask results...')
        print(sinomega[local_mask].sum())
        
        gve_voxel = compute_gve(sc[local_mask], fc[local_mask], omega[local_mask], xpos[local_mask],
                          distance=distance, y_center=y_center, y_size=y_size, tilt_y=tilt_y, z_center=z_center, z_size=z_size, tilt_z=tilt_z, tilt_x=tilt_x,
                          o11=o11, o12=o12, o21=o21, o22=o22,
                          t_x=t_x, t_y=t_y, t_z=t_z, wedge=wedge, chi=chi, wavelength=wavelength)
        
        # iterate through the ubis at this voxel
        for ubi_idx in np.arange(ubis_here.shape[2]):
            ubi = ubis_here[:, :, ubi_idx]

            # assign this UBI to the peaks at this voxel
            j, labels, drlv2 = score_and_assign(ubi, gve_voxel, tol, ubi_idx)
            
            # print(labels.shape)
            # print(drlv2.shape)
            
            grain_peak_mask = labels == ubi_idx
            
            print('Assignment mask results...')
            print(sinomega[local_mask][grain_peak_mask].sum())
            
            grain_npks = np.sum(grain_peak_mask)
            
            if ubi_idx == 0:
                print('In')
                print('UBI Npks')
                print(ubi, grain_npks)
            
            print('Assigned', grain_npks, 'peaks first time out of', local_mask.sum())
            # compute hkl floats
            hkl_double = np.dot(ubi, gve_voxel[:, grain_peak_mask])
            # hkl ints
            hkl = np.round(hkl_double)
            etasign = np.sign(eta[local_mask][grain_peak_mask])
            
            # print('merging')
            # merge the assigned peaks in omega
            merged_sum_intensity, merged_sc, merged_fc, merged_omega, merged_dty, merged_xpos_refined = merge(hkl, etasign, dtyi[local_mask][grain_peak_mask], sum_intensity[local_mask][grain_peak_mask], sc[local_mask][grain_peak_mask], fc[local_mask][grain_peak_mask], omega[local_mask][grain_peak_mask], dty[local_mask][grain_peak_mask], xpos[local_mask][grain_peak_mask])
            
            merged_sinomega = np.sin(np.radians(merged_omega))
            merged_cosomega = np.cos(np.radians(merged_omega))
            
            print('Merged', len(xpos[local_mask][grain_peak_mask]), 'peaks to', len(merged_xpos_refined), 'peaks')
            
            print('Merged peak data for masking...')
            print('merged so sum', merged_sinomega.sum())
            print('merged co sum', merged_cosomega.sum())
            print('merged dty sum', merged_dty.sum())
            
            # print('did merge')
            # re-compute voxel peak mask with merged peaks
            local_mask_of_grain, ydist = get_voxel_mask(y0, xi0, yi0, merged_sinomega, merged_cosomega, merged_dty, ystep)
            print('Got', local_mask_of_grain.sum(), 'masked merged peaks')
            print('ydist sum before merged masking', ydist.sum())
            
            gve_voxel_merged = compute_gve(merged_sc[local_mask_of_grain], merged_fc[local_mask_of_grain], merged_omega[local_mask_of_grain], merged_xpos_refined[local_mask_of_grain],
                  distance=distance, y_center=y_center, y_size=y_size, tilt_y=tilt_y, z_center=z_center, z_size=z_size, tilt_z=tilt_z, tilt_x=tilt_x,
                  o11=o11, o12=o12, o21=o21, o22=o22,
                  t_x=t_x, t_y=t_y, t_z=t_z, wedge=wedge, chi=chi, wavelength=wavelength)
            
            print('gve_voxel_merged sum before merged assignment', gve_voxel_merged.sum())
            
            # print('reassigning')
            # reassign merged g-vectors with smaller merge_tol
            j, labels, drlv2 = score_and_assign(ubi, gve_voxel_merged, merge_tol, ubi_idx)
            grain_peak_mask = labels == ubi_idx
            grain_npks = np.sum(grain_peak_mask)
            print('Assigned', grain_npks, 'merged peaks out of', local_mask_of_grain.sum())
            
            # re-mask g-vectors to those assigned
            gve_grain = gve_voxel_merged[:, grain_peak_mask]
            hkl_double = np.dot(ubi, gve_grain)
            hkl = np.round(hkl_double)
            grain_ydist = ydist[local_mask_of_grain][grain_peak_mask]
            
            # print(grain_npks, 'merged peaks assigned after merging')
            
            # now we're ready to refine!
            
            if grain_npks > 6:
                
                print('Going into refinement...')
                print('ydist sum', grain_ydist.sum())
                print('gve sum', gve_grain.sum())
                print('hkl sum', hkl.sum())
                
                w, ubifit, residuals, rank, sing_vals = weighted_lstsq_ubi_fit(grain_ydist, gve_grain, hkl)
                
                # check the quality of the fit
                worth_fitting = (ubifit is not None) and (rank==3) and (np.linalg.cond(ubifit)<1e14) and (np.linalg.det(ubifit) > 0) and (np.linalg.matrix_rank(ubifit)==3)
                
                # check the unitcell that you create from the UBI
                if worth_fitting:
                    ucell = ubi_to_unitcell(ubifit)
                    _a, _b, _c, al, be, ga = ucell
                    
                    assert _a>0
                    assert _b>0
                    assert _c>0
                    assert al>0 and al<180
                    assert be>0 and be<180
                    assert ga>0 and ga<180
                
                # do we like the quality?
                if worth_fitting:
                    ubi_out = ubifit
                else:
                    ubi_out = ubi.copy()
                    
                if ubi_idx == 0:
                    print('Out')
                    print('UBI Npks')
                    print(ubi_out, grain_npks)
                
                # now reassign the (probably refined) UBI
                
                j, labels, drlv2 = score_and_assign(ubi_out, gve_voxel_merged, merge_tol, ubi_idx)
                grain_peak_mask = labels == ubi_idx
                grain_npks = np.sum(grain_peak_mask)
                
                # if we like the fit quality, redetermine the g-vectors from the new peak mask
                # then fit the strain tensor
                
                if worth_fitting:
                     # go for fancy eps fit
                        
                    print('Strain mask fit sum', grain_peak_mask.sum())
                    print('B0')
                    print(B0)
                    
                    gve_grain_strainfit = gve_voxel_merged[:, grain_peak_mask]
                    hkl_double = np.dot( ubi_out, gve_grain_strainfit ) 
                    hkl = np.round ( hkl_double )
                    ydist = ydist[local_mask_of_grain][grain_peak_mask]
                    
                    print('gve_grain_strainfit sum:', gve_grain_strainfit.sum())
                    print('ydist strainfit sum:', ydist.sum())
                    
                    # get U from UBI without using ImageD11 grain class
                    U = ubi_and_ucell_to_u(ubi_out, ucell)
                    print('U')
                    print(U)
                    gve0 = U @ B0 @ hkl
                    gTg0    = np.sum(gve_grain_strainfit*gve0, axis=0)
                    gTg   = np.sum(gve_grain_strainfit*gve_grain_strainfit, axis=0)
                    directional_strain  = (gTg0/gTg) - 1
                    print('Directional strain')
                    print(directional_strain)
                    
                    kappa = gve_grain_strainfit / gve_norm(gve_grain_strainfit)
                    kx, ky, kz = kappa
                    # M = np.array( [ kx*kx, ky*ky, kz*kz, 2*kx*ky, 2*kx*kz, 2*ky*kz ] ).T
                    M = np.column_stack((kx*kx, ky*ky, kz*kz, 2*kx*ky, 2*kx*kz, 2*ky*kz))

                    w = (1. / (ydist + 1) ).reshape(gve_grain_strainfit.shape[1], 1)
                    # The noise in the directional strain now propagates according to the linear transform
                    gnoise_std = 1e-4
                    a  = np.sum(gve0*(gnoise_std**2)*gve0, axis=0)
                    
                    strain_noise_std = np.sqrt(divide_where(a, gTg**2, out=np.ones_like(gTg), wherearr=gTg))
                    
                    # strain_noise_std = np.sqrt( np.divide(a, gTg**2, out=np.ones_like(gTg), where=gTg!=0) )
                    w = w * (1. / strain_noise_std.reshape(w.shape) )

                    w[directional_strain > np.mean(directional_strain) + np.std(directional_strain)*3.5 ] = 0 # outliers
                    w[directional_strain < np.mean(directional_strain) - np.std(directional_strain)*3.5 ] = 0 # outliers

                    try:
                        w = w / np.max(w)
                        
                        a = w * M
                        b = w.flatten() * directional_strain
                        m, n = a.shape[-2:]
                        rcond = np.finfo(b.dtype).eps * max(n, m)
                        
                        eps_vec = np.linalg.lstsq( a, b, rcond=rcond )[0].T
                        sxx, syy, szz, sxy, sxz, syz = eps_vec
                        eps_tensor = np.array([[sxx, sxy, sxz],[sxy, syy, syz],[sxz, syz, szz]])
                    except:
                        pass
               
                    print(eps_tensor)
                
@numba.njit
def ubi_to_unitcell(ubi):
    # fast numba version, can't use guvec version from tensor_map.py here unfortunately
    mt = np.dot(ubi, ubi.T)
    G = mt
    a, b, c = np.sqrt(np.diag(G))
    al = np.degrees(np.arccos(G[1, 2] / b / c))
    be = np.degrees(np.arccos(G[0, 2] / a / c))
    ga = np.degrees(np.arccos(G[0, 1] / a / b))
    return np.array([a, b, c, al, be, ga])

@numba.njit
def ubi_and_ucell_to_u(ubi, ucell):
    # compute B
    a, b, c = ucell[:3]
    ralpha, rbeta, rgamma = np.radians(ucell[3:])  # radians
    ca = np.cos(ralpha)
    cb = np.cos(rbeta)
    cg = np.cos(rgamma)
    g = np.full((3, 3), np.nan, float)
    g[0, 0] = a * a
    g[0, 1] = a * b * cg
    g[0, 2] = a * c * cb
    g[1, 0] = a * b * cg
    g[1, 1] = b * b
    g[1, 2] = b * c * ca
    g[2, 0] = a * c * cb
    g[2, 1] = b * c * ca
    g[2, 2] = c * c
    gi = np.linalg.inv(g)
    astar, bstar, cstar = np.sqrt(np.diag(gi))
    betas = np.degrees(np.arccos(gi[0, 2] / astar / cstar))
    gammas = np.degrees(np.arccos(gi[0, 1] / astar / bstar))
    
    B = np.zeros((3,3))
    
    B[0, 0] = astar
    B[0, 1] = bstar * np.cos(np.radians(gammas))
    B[0, 2] = cstar * np.cos(np.radians(betas))
    B[1, 0] = 0.0
    B[1, 1] = bstar * np.sin(np.radians(gammas))
    B[1, 2] = -cstar * np.sin(np.radians(betas)) * ca
    B[2, 0] = 0.0
    B[2, 1] = 0.0
    B[2, 2] = 1.0 / c
    
    u = np.dot(B, ubi).T
    return u
                

def call_mine(refine, points_step_space=None):
    # prepare simple numpy array objects
    # pass to Numba function
    
    # get columns for refine.pbpmap
    # which contain the ri, rj points
    # right now it's indexed in step space
    ri_col, rj_col = geometry.step_to_recon(refine.pbpmap.i, refine.pbpmap.j, refine.mask.shape)
    
    if points_step_space is None:
        # get the list of peak indices masks in reconstruction space
        points_recon_space = np.array(np.nonzero(refine.mask)).T
    else:
        # convert input points to reconstruction space, then pass to the function
        points_recon_space = [geometry.step_to_recon(si, sj, refine.mask.shape) for (si, sj) in points_step_space]
    
    # columnfile by [3, 3, (ri, rj)]
    all_pbpmap_ubis = refine.pbpmap.ubi
    
    pars = refine.icolf.parameters.get_parameters()
    
    dummy_var = np.eye(3)
    B0 = unitcell_to_b(refine.ref_ucell.lattice_parameters, dummy_var)
    
    mine(points_recon_space, all_pbpmap_ubis, ri_col, rj_col, refine.sx_grid, refine.sy_grid,
         refine.icolf.sc, refine.icolf.fc, refine.icolf.eta, refine.icolf.sum_intensity, refine.icolf.sinomega, refine.icolf.cosomega, refine.icolf.omega, refine.icolf.dty, refine.icolf.dtyi, refine.icolf.xpos_refined,
         refine.ystep, refine.y0,
         B0,
         pars['distance'], pars['y_center'], pars['y_size'], pars['tilt_y'], pars['z_center'], pars['z_size'], pars['tilt_z'], pars['tilt_x'],
         pars['o11'], pars['o12'], pars['o21'], pars['o22'],
         pars['t_x'], pars['t_y'], pars['t_z'], pars['wedge'], pars['chi'], pars['wavelength']
         )
    


def test_funcs(refine, points):
    # Reference
    
    print('Calling reference')
    start_time = time.perf_counter()
    inputs_axel, results_axel = call_axel(refine, points)
    # [([ubi0, ubi1], [eps0, eps1], [npks0, npks1]), ([ubi0, ubi1], [eps0, eps1], [npks0, npks1])]
    print('In')
    print('UBI Npks')
    print(inputs_axel[0][0][0], inputs_axel[0][1][0])
    print('Out')
    print('UBI Npks')
    print(results_axel[0][0][0], results_axel[0][2][0])
    end_time = time.perf_counter()
    time_taken = end_time - start_time
    print('Reference', time_taken)
    
    print('Calling mine')
    start_time = time.perf_counter()
    call_mine(refine, points)
    # # [([ubi0, ubi1], [eps0, eps1], [npks0, npks1]), ([ubi0, ubi1], [eps0, eps1], [npks0, npks1])]
    # print('In')
    # print('UBI Npks')
    # print(inputs_axel[0][0][0], inputs_axel[0][1][0])
    # print('Out')
    # print('UBI Npks')
    # print(results_axel[0][0][0], results_axel[0][2][0])
    end_time = time.perf_counter()
    time_taken = end_time - start_time
    print('Mine', time_taken)
    
    

def load_data():
    print('Loading data')
    # load dataset
    dset_file = '/data/visitor/ihma460/id11/20240621/PROCESSED_DATA/SS316/JADB/PROCESSED_DATA/SS316/SS316_s3DXRD_z10/SS316_s3DXRD_z10_dataset.h5'
    ds = ImageD11.sinograms.dataset.load(dset_file)
    ds.phases = ds.get_phases_from_disk()
    phase_str = 'Fe'
    
    # load pbp map
    pbpfile = os.path.join(ds.analysispath, ds.dsname + '_pbp.txt')
    pmap = PBPMap(pbpfile)
    
    # set single
    min_unique = 10
    pmap.choose_best(min_unique)
    
    # set up refinement
    refine = PBPRefine(dset=ds, y0=0.0, fpks=0.7, hkl_tol=0.0074, ds_tol=0.005, ifrac=1e-3)
    refine.setmap(pmap)
    refine.loadpeaks(refine.dset.refpeaksfile)
    refine.mask = np.load('mask.npy')
    refine.singlemap = np.load('singlemap.npy')
    ds.update_colfile_pars(refine.icolf, phase_name=phase_str)
    refine.ref_ucell = refine.dset.phases.unitcells[phase_str]
    
    print('Data loaded!')
    
    
    return refine


def main():
    refine = load_data()
    # points = [(344, -212), (323, 147), (-228, 303)]
    points = [(344, -212)]
    
    test_funcs(refine, points)
    # test_funcs(args, [10,50,100,500,1000,5000,10000])
    return

if __name__ == '__main__':
    main()