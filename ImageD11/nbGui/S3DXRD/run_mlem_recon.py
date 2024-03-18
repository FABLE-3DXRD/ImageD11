def read_hdf5(h5name, ginc):
    with h5py.File(h5name, 'r') as hin:
        grains_group = 'grains'
        gg = hin[grains_group][sorted(hin[grains_group].keys(), key=lambda x: int(x))[ginc]]
        # gg = hin[grains_group][str(gid)]
        
        ssino = gg['ssino'][:]
        sinoangles = gg['sinoangles'][:]
        circle_mask = gg['circle_mask'][:]
        y0 = float(gg.attrs['y0'][()])
        
    return ssino, sinoangles, circle_mask, y0


def main(h5name, ginc, dest, pad, niter, dohm, maskcen):
    # read the HDF5 file
    ssino, sinoangles, circle_mask, y0 = read_hdf5(h5name, ginc)
    # feed ssino, sinoangles, mask to run_iradon_mlem
    nthreads = len(os.sched_getaffinity(os.getpid()))
    # print(f"Running on {nthreads} threads")

    recon = ImageD11.sinograms.roi_iradon.run_mlem(ssino, sinoangles, mask=circle_mask, pad=pad, y0=y0, workers=nthreads, niter=niter, apply_halfmask=dohm, mask_central_zingers=maskcen)

    # write result to disk as Numpy array
    # save as TIFF instead?
    np.savetxt(dest, recon)


if __name__ == "__main__":
    # horrible workaround to include id11 code path
    import sys
    
    id11_code_path = sys.argv[1]
    
    sys.path.insert(0, id11_code_path)
    
    import os

    import numpy as np
    import h5py

    import ImageD11.sinograms.roi_iradon
    
    h5name = sys.argv[2]
    ginc = int(sys.argv[3])
    dest = sys.argv[4]
    pad = int(sys.argv[5])
    niter = int(sys.argv[6])
    dohm = sys.argv[7]
    maskcen = sys.argv[8]
    
    if dohm == "Yes":
        apply_halfmask = True
    else:
        apply_halfmask = False
        
    if maskcen == "Yes":
        mask_central_zingers = True
    else:
        mask_central_zingers = False
    
    main(h5name, ginc, dest, pad, niter, apply_halfmask, mask_central_zingers)
