

import sys
import numpy as np
from ImageD11 import grain, parameters, blobcorrector, transform, \
    unitcell, indexing

def tth_ds_max( pars, detector_size ):
    # tth of the four corners.
    corners = [ (0, 0, detector_size[0],detector_size[0]),
                (0,detector_size[1],detector_size[1],0) ]
    tth_c, eta_c = transform.compute_tth_eta( corners, **pars.parameters)
    tthmax = tth_c.max()
    dsmax = 2*np.sin(np.radians(tthmax/2))/pars.get("wavelength")
    return tthmax, dsmax

def inscan( sc, fc, omega, detector_size):
    if sc <= 0 or sc > detector_size[0]:
        return False            
    if fc <= 0 or fc > detector_size[1]:
        return False
    return True

def forwards_project( gr,
                      pars,
                      detector_size,
                      spatial,
                      dsmax,
                      gid):
    latt = "P"
    for k in pars.parameters.keys():
        if k.startswith("cell_lattice"):
            latt = pars.get(k)
    cell = unitcell.unitcell( indexing.ubitocellpars( gr.ubi ), latt )
    ds_hkl = cell.gethkls( dsmax )
    hkls = np.array( [ x[1] for x in ds_hkl] )
    gvecs = np.dot( gr.UB, hkls.T )
    wvln, wedge, chi = [ pars.get(x) for x in ["wavelength", "wedge", "chi"] ]
    tth, (eta1, eta2), (omega1, omega2) = transform.uncompute_g_vectors(
        gvecs, wvln, wedge=wedge, chi=chi)
    osign = float(pars.get("omegasign"))
    if osign != 1:
        print("# Flipping omegasign %f"%(osign))
        np.multiply( omega1, osign, omega1 )
        np.multiply( omega2, osign, omega2 )
    pars.set( "t_x", gr.translation[0] )
    pars.set( "t_y", gr.translation[1] )
    pars.set( "t_z", gr.translation[2] )
    fc1, sc1 = transform.compute_xyz_from_tth_eta(tth, eta1, omega1,
                                                  **pars.parameters)
    fc2, sc2 = transform.compute_xyz_from_tth_eta(tth, eta2, omega2,
                                                  **pars.parameters)
    # Now make a single list of output stuff:
    alltth = np.concatenate( (tth, tth))
    alleta = np.concatenate( (eta1, eta2) )
    allsc =  np.concatenate( (sc1, sc2) )
    allfc =  np.concatenate( (fc1, fc2) )
    sraw = np.zeros(allsc.shape)
    fraw = np.zeros(allsc.shape)
    for i,(s,f) in enumerate(zip( allsc, allfc )):
        sraw[i], fraw[i] = spatial.distort( s, f )
    allomega = np.concatenate( (omega1, omega2) )
    allhkls =  np.concatenate( (hkls, hkls) )
    order = np.argsort( allomega )
    # output ?
    # omega in range
    # peak in detector
    #  etc
    fmt = " % -4d % -3d % -3d % -3d % 8.4f % 7.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f"
    strs = []
    for i in order:
        s, f, o = sraw[i], fraw[i], allomega[i]
        if alltth[i] == 0:
            continue
        if inscan(  s, f, o, detector_size):
            line = fmt%(gid,allhkls[i][0],allhkls[i][1],allhkls[i][2],
                        alltth[i], alleta[i], allomega[i], allsc[i], allfc[i],
                        sraw[i], fraw[i] )
            strs.append( (allomega[i], line) ) # to be omega sortable
    return strs


if __name__=="__main__":
    grains     = grain.read_grain_file( sys.argv[1] )
    pars = parameters.read_par_file( sys.argv[2] ) 
    detector_size = (2048, 2048)
    spline = sys.argv[3]
    if spline == "perfect":
        spatial = blobcorrector.perfect()
    else:
        spatial = blobcorrector.correctorclass( spline )
    outfile = sys.argv[4]
    tthmax, dsmax = tth_ds_max( pars, detector_size )
    print ("# id   h   k   l    tth       eta      omega     sc       fc     s_raw    f_raw")
    peaks = []
    for gid, gr in enumerate(grains):
        newpeaks = forwards_project( gr,
                                     pars,
                                     detector_size,
                                     spatial,
                                     dsmax,
                                     gid )
        peaks += newpeaks
        print ("# Grain",gid,"npks",len(newpeaks))
    peaks.sort()
    out = open(outfile, "w")
    out.write("# id   h   k   l    tth       eta      omega     sc       fc     s_raw    f_raw\n")
    for o, l in peaks:
        out.write(l+"\n")
    out.close()
        









