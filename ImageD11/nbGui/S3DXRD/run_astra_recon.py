def main(h5name, dsfile, group_name):
    print('Loading dataset')
    # load the dataset
    ds = ImageD11.sinograms.dataset.load(dsfile)
    
    print('Loading grain sinograms from disk')
    grainsinos = read_h5(h5name, ds, group_name=group_name)
    
    print('Reconstructing grain sinograms')
    for inc, gs in enumerate(grainsinos):
        gs.recon(method="astra", astra_method="EM_CUDA")
        print('Reconstructed ' + str(inc+1) + '/' + str(len(grainsinos)))

    # mask recon after running
    print('Masking reconstructions')
    for gs in grainsinos:
        gs.recons['astra'] = np.where(gs.recon_mask, gs.recons['astra'], 0)
    
    print('Reconstructions finished')
    print('Writing grains to disk')
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
