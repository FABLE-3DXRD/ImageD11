import os

username = os.environ.get("USER")

id11_code_path = "/home/esrf/" + username + "/Code/ImageD11"

import sys

sys.path.insert(0, id11_code_path)

import numpy as np
import h5py

import ImageD11.sinograms.roi_iradon

def run_iradon_mlem(ssino, sinoangles, mask=None, pad=20, y0=0, workers=1, niter=20, apply_halfmask=False, mask_central_zingers=False):
    outsize = ssino.shape[0] + pad
    if apply_halfmask:
        halfmask = np.zeros_like(ssino)

        halfmask[:len(halfmask)//2-1, :] = 1
        halfmask[len(halfmask)//2-1, :] = 0.5
        
        ssino_to_recon = ssino * halfmask
    else:
        ssino_to_recon = ssino
    
    # Perform iradon transform of grain sinogram, store result (reconstructed grain shape) in g.recon
    recon =  ImageD11.sinograms.roi_iradon.mlem(ssino_to_recon, 
                       theta=sinoangles, 
                       mask=mask,
                       workers=workers,
                       output_size=outsize,
                       projection_shifts=np.full(ssino_to_recon.shape, -y0),
                       niter=niter)
    
    if mask_central_zingers:
        grs = recon.shape[0]
        xpr, ypr = -grs//2 + np.mgrid[:grs, :grs]
        inner_mask_radius = 25
        outer_mask_radius = inner_mask_radius + 2

        inner_circle_mask = (xpr ** 2 + ypr ** 2) < inner_mask_radius ** 2
        outer_circle_mask = (xpr ** 2 + ypr ** 2) < outer_mask_radius ** 2

        mask_ring = inner_circle_mask & outer_circle_mask
        # we now have a mask to apply
        fill_value = np.median(recon[mask_ring])
        recon[inner_circle_mask] = fill_value
    
    return recon

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
    recon = run_iradon_mlem(ssino, sinoangles, mask=circle_mask, pad=pad, y0=y0, workers=nthreads, niter=niter, apply_halfmask=dohm, mask_central_zingers=maskcen)
    # write result to disk as Numpy array
    # save as TIFF instead?
    np.savetxt(dest, recon)
    
if __name__ == "__main__":
    h5name = sys.argv[1]
    ginc = int(sys.argv[2])
    dest = sys.argv[3]
    # y0 = float(sys.argv[4])
    # pad = int(sys.argv[5])
    # niter = int(sys.argv[6])
    # dohm = sys.argv[7]
    # maskcen = sys.argv[8]
    
    pad = int(sys.argv[4])
    niter = int(sys.argv[5])
    dohm = sys.argv[6]
    maskcen = sys.argv[7]
    
    if dohm == "Yes":
        apply_halfmask = True
    else:
        apply_halfmask = False
        
    if maskcen == "Yes":
        mask_central_zingers = True
    else:
        mask_central_zingers = False
    
    main(h5name, ginc, dest, pad, niter, apply_halfmask, mask_central_zingers)
