def read_hdf5(h5name, ginc, ds):
    with h5py.File(h5name, "r") as hin:
        grains_group = hin['grains']
        grain_group = grains_group[str(ginc)]

        # import the grain object
        g = ImageD11.grain.grain.from_h5py_group(grain_group)

        # create the GrainSinogram object
        gs = GrainSinogram.from_h5py_group(grain_group, ds, g)

    return gs


def main(h5name, ginc, dsfile, dest, nthreads):
    # load the dataset
    ds = ImageD11.sinograms.dataset.load(dsfile)
    # read the GrainSinogram object from HDF5:

    grainsino = read_hdf5(h5name, ginc, ds)

    grainsino.recon(method="mlem", workers=nthreads)

    # write result to disk as Numpy array
    # save as TIFF instead?
    np.savetxt(dest, grainsino.recons["mlem"])


if __name__ == "__main__":
    # horrible workaround to include id11 code path
    import sys    

    import numpy as np
    import h5py

    import ImageD11.sinograms.dataset
    from ImageD11.sinograms.sinogram import GrainSinogram
    import ImageD11.grain
    
    h5name = sys.argv[1]
    ginc = int(sys.argv[2])
    dsfile = sys.argv[3]
    dest = sys.argv[4]
    nthreads = int(sys.argv[5])
    
    main(h5name, ginc, dsfile, dest, nthreads)
