def main(h5name, dsfile, group_name):
    # load the dataset
    ds = ImageD11.sinograms.dataset.load(dsfile)

    grainsinos = read_h5(h5name, ds, group_name=group_name)

    for gs in grainsinos:
        gs.recon(method="astra", astra_method="EM_CUDA")

    # mask recon after running
    for gs in grainsinos:
        gs.recons['astra'] = np.where(gs.recon_mask, gs.recons['astra'], 0)

    write_h5(h5name, grainsinos, overwrite_grains=True, group_name=group_name)


if __name__ == "__main__":
    # horrible workaround to include id11 code path
    import sys
    
    from ImageD11.sinograms.sinogram import read_h5, write_h5
    import ImageD11.sinograms.dataset
    import numpy as np
    
    h5name = sys.argv[1]
    dsfile = sys.argv[2]
    group_name = sys.argv[3]
    
    main(h5name, dsfile, group_name)
